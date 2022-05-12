from Bio.SeqFeature import FeatureLocation
from metaerg.run_and_read.abc import AbstractBaseClass
from metaerg.run_and_read.data import MetaergSeqFeature
from metaerg import utils

class Aragorn(AbstractBaseClass):
    def __init__(self, genome, subsystem_hash, force=False, multi_mode=False):
        super().__init__(genome, subsystem_hash, force, multi_mode)

    def __purpose__(self) -> str:
        """Should return the purpose of the tool"""
        return 'tRNA prediction with aragorn'

    def __programs__(self) -> tuple:
        """Should return a tuple with the programs needed"""
        return tuple('aragorn')

    def __result_files__(self) -> tuple:
        """Should return a tuple with the result files (Path objects) created by the programs"""
        return tuple(self.spawn_file('aragorn'))

    def __run_programs__(self):
        """Should execute the helper programs to complete the analysis"""
        fasta_file = self.make_masked_contig_fasta_file()
        aragorn_file = self.spawn_file("aragorn")
        utils.run_external(f'aragorn -l -t -gc{self.genome.translation_table} {fasta_file} -w -o {aragorn_file}')

    def __read_results__(self) -> int:
        """Should parse the result files and return the # of positives"""
        trna_count = 0
        current_contig = None
        with open(self.spawn_file("aragorn")) as aragorn_handle:
            for line in aragorn_handle:
                if line.startswith('>end'):
                    break
                if line.startswith('>'):
                    current_contig = self.genome.contigs[line[1:].strip()]
                if not current_contig:
                    continue
                words = line.split()
                if len(words) < 5:
                    continue
                trna_count += 1
                strand = 1
                if words[2].startswith('c'):
                    strand = -1
                    words[2] = words[2][1:]
                pos_str = words[2][1:-1].split(',')
                start = max(0, int(pos_str[0]) - 1)
                end = min(len(current_contig.seq), int(pos_str[1]))
                f = MetaergSeqFeature(current_contig)
                f.location = FeatureLocation(start, end, strand=strand)
                f.type = 'tRNA'
                f.inference = 'aragorn'
                f.description = f'{words[1]}-{words[4]}'
                current_contig.features.append(f)
        return trna_count

    def __repr__(self):
        return f'aragorn({self.genome.name}, force={self.force}, multi_mode={self.multi_mode})'
