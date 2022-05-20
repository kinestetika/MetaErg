from metaerg.run_and_read.abc import Annotator, ExecutionEnvironment, register
from metaerg.run_and_read.data_model import FeatureType


@register
class GeneWriter(Annotator):
    def __init__(self, genome, exec: ExecutionEnvironment):
        super().__init__(genome, exec)
        self.cds_aa_file = self.spawn_file('cds.faa')
        self.rna_nt_file = self.spawn_file('rna.nt')
        self.pipeline_position = 66

    def __repr__(self):
        return f'GeneWriter({self.genome}, {self.exec})'

    def _purpose(self) -> str:
        """Should return the purpose of the tool"""
        return 'Write files for proteins and rna genes.'

    def _programs(self) -> tuple:
        """Should return a tuple with the programs needed"""
        return tuple()

    def _result_files(self) -> tuple:
        """Should return a tuple with the result files (Path objects) created by the programs"""
        return self.cds_aa_file, self.rna_nt_file

    def _run_programs(self):
        """No helper progarms needed, simply writes genes to files."""
        self.genome.write_fasta_files(self.cds_aa_file, target=FeatureType.CDS)
        self.genome.write_fasta_files(self.rna_nt_file, target=(FeatureType.rRNA, FeatureType.ncRNA))

    def _read_results(self) -> int:
        return 0
