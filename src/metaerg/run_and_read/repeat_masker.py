import shutil
from pathlib import Path
from metaerg.run_and_read.data_model import MetaergSeqFeature, MetaergSeqRecord, FeatureType
from metaerg.run_and_read.abc import Annotator, ExecutionEnvironment, register
from metaerg import utils


@register
class RepeatMasker(Annotator):
    def __init__(self, genome, exec_env: ExecutionEnvironment):
        super().__init__(genome, exec_env)
        self.repeatmasker_file = self.spawn_file('repeatmasker')
        self._pipeline_position = 51
        self._purpose = 'tandem repeat prediction with trf'
        self._programs = ('build_lmer_table', 'RepeatScout', 'filter-stage-1.prl', 'RepeatMasker')
        self._result_files = (self.repeatmasker_file,)

    def _run_programs(self):
        """Executes the helper programs to complete the analysis"""
        fasta_file, = self.genome.write_fasta_files(self.spawn_file('masked'), masked=True)
        lmer_table_file = self.spawn_file('lmer-table')
        repeatscout_file_raw = self.spawn_file('repeatscout-raw')
        repeatscout_file_filtered = self.spawn_file('repeatscout-filtered')

        utils.run_external(f'build_lmer_table -sequence {fasta_file} -freq {lmer_table_file}')
        utils.run_external(f'RepeatScout -sequence {fasta_file} -output {repeatscout_file_raw} -freq {lmer_table_file}')
        with open(repeatscout_file_filtered, 'w') as output, open(repeatscout_file_raw) as input:
            utils.run_external('filter-stage-1.prl', stdin=input, stdout=output)
        utils.run_external(f'RepeatMasker -pa {self.exec.cpus_per_genome} -lib {repeatscout_file_filtered} -dir . {fasta_file}')
        repeatmasker_output_file = Path(f'{fasta_file.name}.out')  # nothing we can do about that
        shutil.move(repeatmasker_output_file, self.repeatmasker_file)
        for file in Path.cwd().glob(f'{fasta_file.name}.*'):
            if file.is_dir():
                shutil.rmtree(file)
            else:
                file.unlink()

    def _read_results(self) -> int:
        """Should parse the result files and return the # of positives. Repeatmasker finds two types of repeats:
            (1) simple repeats, these are consecutive
            (2) unspecified repeats, these occur scattered and are identified by an id in words[9]. We only
            add those when they occur 10 or more times."""
        repeat_count = 0
        repeat_hash = dict()
        with open(self.repeatmasker_file) as repeatmasker_handle:
            for line in repeatmasker_handle:
                words = line.split()
                if len(words) < 11:
                    continue
                contig: MetaergSeqRecord = self.genome.contigs[words[4]]
                if 'Simple_repeat' == words[10]:
                    repeat_count += 1
                    feature = contig.spawn_feature(int(words[5]) - 1, int(words[6]), -1 if 'C' == words[8] else 1,
                                            FeatureType.repeat, inference='repeatmasker')
                    feature.notes.add(f'repeat {words[9]}')
                else:
                    repeat_list = repeat_hash.setdefault(words[9], list())
                    repeat_list.append({'start': int(words[5]) - 1,
                                        'end': int(words[6]),
                                        'strand': -1 if 'C' == words[8] else 1,
                                        'type': FeatureType.repeat,
                                        'inference': 'repeatmasker'})
        for repeat_list in repeat_hash.values():
            if len(repeat_list) >= 10:
                for f in repeat_list:
                    repeat_count += 1
                    feature = contig.spawn_feature(**f)
                    feature.notes.add(f' (occurs {len(repeat_list)}x)')
        return repeat_count
