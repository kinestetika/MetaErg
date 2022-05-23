from metaerg.run_and_read.abc import Annotator, ExecutionEnvironment, register
from metaerg.run_and_read.data_model import FeatureType


@register
class GeneWriter(Annotator):
    def __init__(self, genome, exec: ExecutionEnvironment):
        super().__init__(genome, exec)
        self.cds_aa_file = self.spawn_file('cds.faa')
        self.rna_nt_file = self.spawn_file('rna.nt')
        self._pipeline_position = 66
        self._purpose = 'Write files for proteins and rna genes.'
        self._result_files = (self.cds_aa_file, self.rna_nt_file)

    def _run_programs(self):
        """No helper progarms needed, simply writes genes to files."""
        self.genome.write_fasta_files(self.cds_aa_file, target=FeatureType.CDS)
        self.genome.write_fasta_files(self.rna_nt_file, target=(FeatureType.rRNA, FeatureType.ncRNA))

    def _read_results(self) -> int:
        return 0
