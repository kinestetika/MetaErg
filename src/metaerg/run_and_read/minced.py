from metaerg.run_and_read.data_model import MetaergSeqRecord, FeatureType
from metaerg.run_and_read.abc import Annotator, register, ExecutionEnvironment
from metaerg import utils

@register
class Minced(Annotator):
    def __init__(self, genome, exec_env: ExecutionEnvironment):
        super().__init__(genome, exec_env)
        self.minced_file = self.spawn_file('minced')
        self._pipeline_position = 1
        self._purpose = 'CRISPR prediction with minced'
        self._programs = ('minced',)
        self._result_files = (self.minced_file,)

    def _run_programs(self):
        """Executes the helper programs to complete the analysis"""
        fasta_file = self.genome.write_fasta_files(self.spawn_file('masked'), masked=True)
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
                        contig.spawn_feature(int(start) - 1, int(end), 1 if '+' == strand else -1,
                                             FeatureType.crispr_repeat, inference = 'minced')
        return crispr_region_count
