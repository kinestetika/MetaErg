from Bio.SeqFeature import FeatureLocation
from metaerg.run_and_read.data_model import MetaergSeqFeature
from metaerg.run_and_read.data_model import MetaergSeqRecord
from metaerg.run_and_read import abc
from metaerg import utils


class TandemRepeatFinder(abc.AbstractBaseClass):
    def __init__(self, genome, exec_env: abc.ExecutionEnvironment):
        super().__init__(genome, exec_env)
        self.trf_file = self.spawn_file('tandem-repeat-finder')

    def __repr__(self):
        return f'TandemRepeatFinder({self.genome}, {self.exec})'

    def _purpose(self) -> str:
        """Should return the purpose of the tool"""
        return 'tandem repeat prediction with trf'

    def _programs(self) -> tuple:
        """Should return a tuple with the programs needed"""
        return 'trf',

    def _result_files(self) -> tuple:
        """Should return a tuple with the result files (Path objects) created by the programs"""
        return self.trf_file,

    def _run_programs(self):
        """Should execute the helper programs to complete the analysis"""
        fasta_file = self.genome.make_masked_contig_fasta_file(self.spawn_file('masked'))
        with open(self.trf_file, 'w') as output:
            utils.run_external(f'trf {fasta_file} 2 7 7 80 10 50 500 -d -h -ngs', stdout=output)

    def _read_results(self) -> int:
        """Should parse the result files and return the # of positives"""
        tr_count = 0
        with open(self.trf_file) as trf_handle:
            for line in trf_handle:
                if line.startswith("@"):
                    contig: MetaergSeqRecord = self.genome.contigs[line[1:].strip()]
                    continue
                if not contig:
                    continue
                words = line.split()
                f = MetaergSeqFeature(contig)
                f.type = 'repeat_region'
                f.location = FeatureLocation(int(words[0]) - 1, int(words[1]))
                f.inference = 'tandem-repeat-finder'
                f.note = f'period size {words[2]}; copies {words[3]}'
                contig.features.append(f)
                tr_count += 1
        return tr_count
