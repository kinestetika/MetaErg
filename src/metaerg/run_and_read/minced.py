from metaerg.run_and_read.abc import AbstractBaseClass
from metaerg.run_and_read.data import MetaergSeqRecord
from metaerg.run_and_read.data import MetaergSeqFeature
from metaerg import utils

class Minced(AbstractBaseClass):
    def __init__(self, genome, subsystem_hash, force=False, multi_mode=False):
        super().__init__(genome, subsystem_hash, force, multi_mode)

    def __purpose__(self) -> str:
        """Should return the purpose of the tool"""
        return 'CRISPR prediction with minced'

    def __programs__(self) -> tuple:
        """Should return a tuple with the programs needed"""
        return tuple('minced')

    def __result_files__(self) -> tuple:
        """Should return a tuple with the result files (Path objects) created by the programs"""
        return tuple(self.spawn_file('minced'))

    def __run_programs__(self):
        """Should execute the helper programs to complete the analysis"""
        fasta_file = self.make_masked_contig_fasta_file()
        utils.run_external(f'minced -gffFull {fasta_file} {self.spawn_file("minced")}')

    def __read_results__(self) -> int:
        """Should parse the result files and return the # of positives"""
        crispr_region_count = 0
        with open(self.spawn_file("minced")) as crispr_handle:
            for line in crispr_handle:
                if line.startswith('#'):
                    continue
                words = line.split('\t')
                if len(words) < 9:
                    continue
                if 'repeat_region' == words[2]:
                    crispr_region_count += 1
                elif 'repeat_unit' == words[2]:
                    contig:MetaergSeqRecord = self.contig_dictionary[words[0]]
                    feature = MetaergSeqFeature(parent=contig, gff_line=line)
                    feature.type = 'crispr_repeat'
                    contig.features.append(feature)
        return crispr_region_count

    def __repr__(self):
        return f'minced({self.genome_name}, force={self.force}, multi_mode={self.multi_mode})'


