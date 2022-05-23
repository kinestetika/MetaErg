import shutil
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor
from metaerg.run_and_read.data_model import MetaergSeqFeature, FeatureType
from metaerg.run_and_read.abc import Annotator, ExecutionEnvironment, register
from metaerg import utils


@register
class SignalP(Annotator):
    def __init__(self, genome, exec_env: ExecutionEnvironment):
        super().__init__(genome, exec_env)
        self.signalp_file = self.spawn_file('signalp')
        self._pipeline_position = 121
        self._purpose = 'signal peptide prediction with signalp'
        self._programs = ('signalp6',)
        self._result_files = (self.signalp_file,)

    def _run_programs(self):
        """Executes the helper programs to complete the analysis"""
        cds_aa_file = self.spawn_file('cds.faa')
        if self.exec.cpus_per_genome > 1:
            split_fasta_files = self.genome.write_fasta_files(cds_aa_file, self.exec.cpus_per_genome, target=FeatureType.CDS)
            split_signalp_files = [Path(self.signalp_file.parent, f'{self.signalp_file.name}.{i}')
                                   for i in range(len(split_fasta_files))]
            with ProcessPoolExecutor(max_workers=self.exec.cpus_per_genome) as executor:
                for split_input, split_output in zip(split_fasta_files, split_signalp_files):
                    executor.submit(utils.run_external, f'signalp6 --fastafile {split_input} --output_dir '
                                                        f'{split_output} --format none --organism other')

            self.signalp_file.mkdir()
            with open(Path(self.signalp_file, 'prediction_results.txt'), 'wb') as output:
                for split_cds_aa_file, split_signalp_dir in zip(split_fasta_files, split_signalp_files):
                    signalp_result_file = Path(split_signalp_dir, 'prediction_results.txt')
                    if signalp_result_file.exists():
                        with open(signalp_result_file, 'rb') as input:
                            shutil.copyfileobj(input, output)
                    else:
                        utils.log(f'({self.genome.id}) WARNING - missing part of signalp output!')
                    shutil.rmtree(split_signalp_dir)
                    split_cds_aa_file.unlink()
        else:
            utils.run_external(f'signalp6 --fastafile {cds_aa_file} --output_dir {self.signalp_file} --format none '
                               f'--organism other')

    def _read_results(self) -> int:
        """Should parse the result files and return the # of positives"""
        count = 0
        with open(Path(self.signalp_file, 'prediction_results.txt')) as signalp_handle:
            for line in signalp_handle:
                if line.startswith("#"):
                    continue
                words = line.split("\t")
                if "OTHER" == words[1]:
                    continue
                feature: MetaergSeqFeature = self.genome.get_feature(words[0].split()[0])
                feature.signal_peptide = words[1]
                count += 1
        return count
