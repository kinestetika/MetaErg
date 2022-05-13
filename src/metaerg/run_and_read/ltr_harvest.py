from Bio.SeqFeature import FeatureLocation
from metaerg.run_and_read.data import MetaergSeqFeature
from metaerg.run_and_read.data import MetaergSeqRecord
from metaerg.run_and_read import abc
from metaerg import utils

class LTRHarvest(abc.AbstractBaseClass):
    def __init__(self, genome, exec:abc.ExecutionEnvironment):
        super().__init__(genome, exec)

    def __repr__(self):
        return f'LTRHarvest({self.genome}, {self.exec})'

    def __purpose__(self) -> str:
        """Should return the purpose of the tool"""
        return 'retrotransposon prediction with ltrharvest'

    def __programs__(self) -> tuple:
        """Should return a tuple with the programs needed"""
        return 'ltr_index', 'ltr_harvest'

    def __result_files__(self) -> tuple:
        """Should return a tuple with the result files (Path objects) created by the programs"""
        return self.spawn_file('aragorn'),

    def __run_programs__(self):
        """Should execute the helper programs to complete the analysis"""
        fasta_file = self.make_masked_contig_fasta_file()
        ltr_index_file = self.spawn_file('ltr_index')
        ltr_harvest_file = self.spawn_file('ltr_harvest')

        utils.run_external(f'gt suffixerator -db {fasta_file} -indexname {ltr_index_file} -tis -suf -lcp -des -ssp -sds -dna')
        utils.run_external(f'gt ltrharvest -index {ltr_index_file} -gff3 {ltr_harvest_file} -seqids')
        # remove index files
        for file in ltr_harvest_file.parent.glob(f'{self.genome.name}.ltr_index*'):
            file.unlink()

    def __read_results__(self) -> int:
        """Should parse the result files and return the # of positives"""
        retrotransposon_count = 0
        with open(self.spawn_file('ltr_harvest')) as ltr_handle:
            for line in ltr_handle:
                if line.startswith('#'):
                    continue
                words = line.split('\t')
                if len(words) < 9:
                    continue
                if 'repeat_region' == words[2]:
                    retrotransposon_count += 1
                    contig: MetaergSeqRecord = self.genome.contigs[words[0]]
                    feature = MetaergSeqFeature(contig, gff_line=line)
                    feature.type = 'retrotransposon'
                    contig.features.append(feature)
        return retrotransposon_count