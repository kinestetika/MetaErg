import shutil
from pathlib import Path
from Bio.SeqFeature import FeatureLocation
from metaerg.run_and_read.data_model import MetaergSeqFeature
from metaerg.run_and_read.data_model import MetaergSeqRecord
from metaerg.run_and_read import abc
from metaerg import utils


class RepeatMasker(abc.AbstractBaseClass):
    def __init__(self, genome, exec_env: abc.ExecutionEnvironment):
        super().__init__(genome, exec_env)
        self.repeatmasker_file = self.spawn_file('repeatmasker')

    def __repr__(self):
        return f'TandemRepeatFinder({self.genome}, {self.exec})'

    def _purpose(self) -> str:
        """Should return the purpose of the tool"""
        return 'tandem repeat prediction with trf'

    def _programs(self) -> tuple:
        """Should return a tuple with the programs needed"""
        return 'build_lmer_table', 'RepeatScout', 'filter-stage-1.prl', 'RepeatMasker'

    def _result_files(self) -> tuple:
        """Should return a tuple with the result files (Path objects) created by the programs"""
        return self.repeatmasker_file,

    def _run_programs(self):
        """Should execute the helper programs to complete the analysis"""
        fasta_file = self.genome.make_masked_contig_fasta_file(self.spawn_file('masked'))
        lmer_table_file = self.spawn_file('lmer-table')
        repeatscout_file_raw = self.spawn_file('repeatscout-raw')
        repeatscout_file_filtered = self.spawn_file('repeatscout-filtered')

        utils.run_external(f'build_lmer_table -sequence {fasta_file} -freq {lmer_table_file}')
        utils.run_external(f'RepeatScout -sequence {fasta_file} -output {repeatscout_file_raw} -freq {lmer_table_file}')
        with open(repeatscout_file_filtered, 'w') as output, open(repeatscout_file_raw) as input:
            utils.run_external('filter-stage-1.prl', stdin=input, stdout=output)
        utils.run_external(f'RepeatMasker -pa {self.exec.threads} -lib {repeatscout_file_filtered} -dir . {fasta_file}')
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
                location = FeatureLocation(int(words[5]) - 1, int(words[6]), strand=-1 if 'C' == words[8] else 1)
                contig: MetaergSeqRecord = self.genome.contigs[words[4]]
                feature = MetaergSeqFeature(contig, 'repeat_region', location, 'repeatmasker')
                if 'Simple_repeat' == words[10]:
                    repeat_count += 1
                    contig.features.append(feature)
                    feature.note = f'repeat {words[9]}'
                else:
                    try:
                        repeat_list = repeat_hash[words[9]]
                    except KeyError:
                        repeat_list = []
                        repeat_hash[words[9]] = repeat_list
                    repeat_list.append(feature)
        for repeat_list in repeat_hash.values():
            if len(repeat_list) >= 10:
                for f in repeat_list:
                    repeat_count += 1
                    contig.features.append(f)
                    f.note = f' (occurs {len(repeat_list)}x)'
        return repeat_count
