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
                if line.startswith("#"):
                    continue
                words = line.split('\t')
                if len(words) < 9:
                    continue
                contig: MetaergSeqRecord = self.genome.contigs[words[0]]
                location = FeatureLocation(int(words[3]) - 1, int(words[4]), strand=-1 if '+' == words[6] else 1)
                feature = contig.spawn_feature('CDS', location, 'prodigal')
                if 'partial=01' in words[8] or 'partial=01' in words[8] or 'partial=11' in words[8]:
                    feature.note = 'partial'
                cds_found += 1
        return cds_found
