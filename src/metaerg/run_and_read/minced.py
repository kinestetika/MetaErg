from Bio.SeqFeature import FeatureLocation
from metaerg.run_and_read.data_model import MetaergSeqRecord
from metaerg.run_and_read import abc
from metaerg import utils


class Minced(abc.AbstractBaseClass):
    def __init__(self, genome, exec_env: abc.ExecutionEnvironment):
        super().__init__(genome, exec_env)
        self.minced_file = self.spawn_file('minced')

    def __repr__(self):
        return f'Minced({self.genome}, {self.exec})'

    def __purpose__(self) -> str:
        """Should return the purpose of the tool"""
        return 'CRISPR prediction with minced'

    def __programs__(self) -> tuple:
        """Should return a tuple with the programs needed"""
        return 'minced',

    def __result_files__(self) -> tuple:
        """Should return a tuple with the result files (Path objects) created by the programs"""
        return self.minced_file,

    def __run_programs__(self):
        """Should execute the helper programs to complete the analysis"""
        fasta_file = self.genome.make_masked_contig_fasta_file(self.spawn_file('masked'))
        utils.run_external(f'minced -gffFull {fasta_file} {self.minced_file}')

    def __read_results__(self) -> int:
        """Should parse the result files and return the # of positives"""
        crispr_region_count = 0
        with open(self.minced_file) as crispr_handle:
            for line in crispr_handle:
                if line.startswith('#'):
                    continue
                words = line.split('\t')
                if len(words) < 9:
                    continue
                if 'repeat_region' == words[2]:
                    crispr_region_count += 1
                elif 'repeat_unit' == words[2]:
                    contig: MetaergSeqRecord = self.genome.contigs[words[0]]
                    location = FeatureLocation(int(words[3]) - 1, int(words[4]), strand=-1 if '+' == words[6] else 1)
                    contig.spawn_feature('crispr_repeat', location, 'minced')
        return crispr_region_count
