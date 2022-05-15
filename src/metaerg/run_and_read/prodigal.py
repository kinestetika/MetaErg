from Bio.SeqFeature import FeatureLocation
from metaerg.run_and_read.data_model import MetaergSeqRecord
from metaerg.run_and_read import abc
from metaerg import utils


class Prodigal(abc.AbstractBaseClass):
    def __init__(self, genome, exec_env: abc.ExecutionEnvironment):
        super().__init__(genome, exec_env)
        self.prodigal_file = self.spawn_file('prodigal')

    def __repr__(self):
        return f'Prodigal({self.genome}, {self.exec})'

    def __purpose__(self) -> str:
        """Should return the purpose of the tool"""
        return 'coding sequence prediction with prodigal'

    def __programs__(self) -> tuple:
        """Should return a tuple with the programs needed"""
        return 'prodigal',

    def __result_files__(self) -> tuple:
        """Should return a tuple with the result files (Path objects) created by the programs"""
        return self.prodigal_file,

    def __run_programs__(self):
        """Should execute the helper programs to complete the analysis"""
        fasta_file = self.genome.make_masked_contig_fasta_file(self.spawn_file('masked'))
        utils.run_external(f'prodigal -g {self.genome.translation_table} -m -f gff -q -i {fasta_file} -o '
                           f'{self.prodigal_file}')

    def __read_results__(self) -> int:
        """Should parse the result files and return the # of positives."""
        cds_found = 0
        with open(self.prodigal_file) as prodigal_handle:
            for line in prodigal_handle:
                words = line.strip().split('\t')
                match(words):
                    case [str(word), *_] if word.startswith('#'):
                        continue
                    case [contig_name, _, _, start, end, _, strand, _, attributes]:
                        cds_found += 1
                        contig: MetaergSeqRecord = self.genome.contigs[contig_name]
                        location = FeatureLocation(int(start) - 1, int(end), strand=-1 if '+' == strand else 1)
                        feature = contig.spawn_feature('CDS', location, 'prodigal')
                        if 'partial=01' in attributes or 'partial=01' in attributes or 'partial=11' in attributes:
                            feature.note = 'partial'
        return cds_found
