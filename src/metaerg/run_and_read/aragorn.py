import re
from Bio.SeqFeature import FeatureLocation
from metaerg.run_and_read.data_model import MetaergSeqFeature
from metaerg.run_and_read import abc
from metaerg import utils

class Aragorn(abc.AbstractBaseClass):
    def __init__(self, genome, exec:abc.ExecutionEnvironment):
        super().__init__(genome, exec)
        self.aragorn_file = self.spawn_file("aragorn")

    def __repr__(self):
        return f'Aragorn({self.genome}, {self.exec})'

    def __purpose__(self) -> str:
        """Should return the purpose of the tool"""
        return 'tRNA prediction with aragorn'

    def __programs__(self) -> tuple:
        """Should return a tuple with the programs needed"""
        return 'aragorn',

    def __result_files__(self) -> tuple:
        """Should return a tuple with the result files (Path objects) created by the programs"""
        return self.aragorn_file,

    def __run_programs__(self):
        """Should execute the helper programs to complete the analysis"""
        fasta_file = self.genome.make_masked_contig_fasta_file(self.spawn_file('masked'))
        utils.run_external(f'aragorn -l -t -gc{self.genome.translation_table} {fasta_file} -w -o {self.aragorn_file}')

    def __read_results__(self) -> int:
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
                        end = min(len(current_contig.seq), int(coord_match.group(3)))
                        f = current_contig.spawn_feature('tRNA', FeatureLocation(start, end, strand=strand), 'aragorn')
                        f.description = f'{trna}-{codon}'
        return trna_count
