import re
from metaerg.run_and_read.data_model import FeatureType
from metaerg.run_and_read.abc import Annotator, ExecutionEnvironment, register
from metaerg import utils


@register
class Aragorn(Annotator):
    def __init__(self, genome, exec: ExecutionEnvironment):
        super().__init__(genome, exec)
        self.aragorn_file = self.spawn_file("aragorn")
        self.pipeline_position = 11

    def __repr__(self):
        return f'Aragorn({self.genome}, {self.exec})'

    def _purpose(self) -> str:
        """Should return the purpose of the tool"""
        return 'tRNA prediction with aragorn'

    def _programs(self) -> tuple:
        """Should return a tuple with the programs needed"""
        return 'aragorn',

    def _result_files(self) -> tuple:
        """Should return a tuple with the result files (Path objects) created by the programs"""
        return self.aragorn_file,

    def _run_programs(self):
        """Should execute the helper programs to complete the analysis"""
        fasta_file = self.genome.write_fasta_files(self.spawn_file('masked'), masked=True)
        utils.run_external(f'aragorn -l -t -gc{self.genome.translation_table} {fasta_file} -w -o {self.aragorn_file}')

    def _read_results(self) -> int:
        """Should parse the result files and return the # of positives"""
        trna_count = 0
        current_contig = None
        coord_regexp = re.compile(r'(c*)\[(\d+),(\d+)]')
        with open(self.aragorn_file) as aragorn_handle:
            for line in aragorn_handle:
                words = line.strip().split()
                match words:
                    case ['>end']:
                        break
                    case [contig_name] if contig_name.startswith('>'):
                        current_contig = self.genome.contigs[contig_name[1:]]
                    case [_, trna, coordinates, _, codon]:
                        trna_count += 1
                        coord_match = coord_regexp.fullmatch(coordinates)
                        strand = -1 if 'c' == coord_match.group(1) else 1
                        start = max(0, int(coord_match.group(2)) - 1)
                        end = min(len(current_contig.sequence), int(coord_match.group(3)))
                        f = current_contig.spawn_feature(start, end, strand, FeatureType.tRNA,
                                                         inference='aragorn')
                        f.description = f'{trna}-{codon}'
        return trna_count
