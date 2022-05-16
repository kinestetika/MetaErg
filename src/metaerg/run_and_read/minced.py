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

    def _purpose(self) -> str:
        """Should return the purpose of the tool"""
        return 'CRISPR prediction with minced'

    def _programs(self) -> tuple:
        """Should return a tuple with the programs needed"""
        return 'minced',

    def _result_files(self) -> tuple:
        """Should return a tuple with the result files (Path objects) created by the programs"""
        return self.minced_file,

    def _run_programs(self):
        """Should execute the helper programs to complete the analysis"""
        fasta_file = self.genome.make_masked_contig_fasta_file(self.spawn_file('masked'))
        utils.run_external(f'minced -gffFull {fasta_file} {self.minced_file}')

    def _read_results(self) -> int:
        """Should parse the result files and return the # of positives"""
        crispr_region_count = 0
        with open(self.minced_file) as crispr_handle:  # this file is in gff format
            for line in crispr_handle:
                words = line.strip().split('\t')
                match words:
                    case [str(word), *_] if word.startswith('#'):
                        continue
                    case [_, _, 'repeat_region', _, _, _, _, _, _]:
                        crispr_region_count += 1
                    case [contig_name, _, 'repeat_unit', start, end, _, strand, _, _]:
                        contig: MetaergSeqRecord = self.genome.contigs[contig_name]
                        location = FeatureLocation(int(start) - 1, int(end), strand=-1 if '+' == strand else 1)
                        contig.spawn_feature('crispr_repeat', location, 'minced')
        return crispr_region_count
