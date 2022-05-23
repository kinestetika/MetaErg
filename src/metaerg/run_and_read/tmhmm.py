import shutil
from metaerg.run_and_read.abc import Annotator, ExecutionEnvironment, register
from metaerg import utils


@register
class TMHMM(Annotator):
    def __init__(self, genome, exec_env: ExecutionEnvironment):
        super().__init__(genome, exec_env)
        self.tmhmm_file = self.spawn_file('signalp')
        self._pipeline_position = 111
        self._purpose = 'transmembrane helix prediction with tmhmm'
        self._programs = ('tmhmm',)
        self._result_files = (self.tmhmm_file,)

    def _run_programs(self):
        """Executes the helper programs to complete the analysis"""
        cds_aa_file = self.spawn_file('cds.faa')
        with open(self.tmhmm_file, 'w') as output, open(cds_aa_file) as input:
            utils.run_external('tmhmm', stdin=input, stdout=output)
        # this is not thread-safe:
        for file in self.tmhmm_file.parent.glob(f'TMHMM_*'):
            if file.is_dir():
                shutil.rmtree(file)

    def _read_results(self) -> int:
        """Should parse the result files and return the # of positives"""
        count = 0
        current_feature = None
        current_txt = ""
        feature_tmh_count = 0
        with open(self.tmhmm_file) as tmhmm_handle:
            for line in tmhmm_handle:
                words = line.strip().split('\t')
                match words:
                    case[str(first_text), *_] if first_text.startswith('#'):
                        continue
                    case [_, _, 'TMhelix', start, end]:
                        feature_tmh_count += 1
                        current_txt += f'{start}-{end},'
                    case [feature_name, _, orientation, _, _] if orientation in ('inside', 'outside'):
                        if not current_feature or current_feature.id != feature_name:
                            if feature_tmh_count:
                                current_feature.transmembrane_helixes = " ".join((feature_tmh_count, current_txt[:-1]))
                                count += 1
                            new_feature = self.genome.get_feature(feature_name)
                            current_feature = new_feature
                            feature_tmh_count = 0
                            current_txt = 'i,' if 'inside' == orientation else 'o,'
            if feature_tmh_count:
                current_feature.transmembrane_helixes = " ".join((feature_tmh_count, current_txt[:-1]))
                count += 1
        return count
