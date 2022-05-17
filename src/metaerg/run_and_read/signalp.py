import shutil
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor
from metaerg.run_and_read.data_model import MetaergSeqFeature
from metaerg.run_and_read import abc
from metaerg import utils


class SignalP(abc.AbstractBaseClass):
    def __init__(self, genome, exec_env: abc.ExecutionEnvironment):
        super().__init__(genome, exec_env)
        self.signalp_file = self.spawn_file('signalp')

    def __repr__(self):
        return f'SignalP({self.genome}, {self.exec})'

    def _purpose(self) -> str:
        """Should return the purpose of the tool"""
        return 'signal peptide prediction with signalp'

    def _programs(self) -> tuple:
        """Should return a tuple with the programs needed"""
        return 'signalp6',

    def _result_files(self) -> tuple:
        """Should return a tuple with the result files (Path objects) created by the programs"""
        return self.signalp_file,

    def _run_programs(self):
        """Should execute the helper programs to complete the analysis"""
        cds_aa_file = self.spawn_file('cds.faa')
        if self.exec.threads > 1:
            split_fasta_files = self.genome.make_split_fasta_files(cds_aa_file, self.exec.threads)
            split_signalp_files = [Path(self.signalp_file.parent, f'{self.signalp_file.name}.{i}')
                                   for i in range(len(split_fasta_files))]
            with ProcessPoolExecutor(max_workers=self.exec.threads) as executor:
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
