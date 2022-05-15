from Bio.SeqFeature import FeatureLocation
from metaerg.run_and_read.data_model import MetaergSeqRecord
from metaerg.run_and_read import abc
from metaerg import utils


class LTRHarvest(abc.AbstractBaseClass):
    def __init__(self, genome, exec: abc.ExecutionEnvironment):
        super().__init__(genome, exec)
        self.ltr_harvest_file = self.spawn_file('ltr_harvest')

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
        return self.ltr_harvest_file,

    def __run_programs__(self):
        """Should execute the helper programs to complete the analysis"""
        fasta_file = self.genome.make_masked_contig_fasta_file(self.spawn_file('masked'))
        ltr_index_file = self.spawn_file('ltr_index')

        utils.run_external(f'gt suffixerator -db {fasta_file} -indexname {ltr_index_file} -tis -suf -lcp '
                           f'-des -ssp -sds -dna')
        utils.run_external(f'gt ltrharvest -index {ltr_index_file} -gff3 {self.ltr_harvest_file} -seqids')
        # remove index files
        for file in ltr_index_file.parent.glob(f'{self.genome.name}.ltr_index*'):
            file.unlink()

    def __read_results__(self) -> int:
        """Should parse the result files and return the # of positives"""
        retrotransposon_count = 0
        with open(self.ltr_harvest_file) as ltr_handle:
            for line in ltr_handle:
                words = line.strip().split('\t')
                match words:
                    case [str(word), *_] if word.startswith('#'):
                        continue
                    case [contig_name, _, 'repeat_region', start, end, score, strand, frame, _]:
                        retrotransposon_count += 1
                        contig: MetaergSeqRecord = self.genome.contigs[contig_name]
                        location = FeatureLocation(int(start) - 1, int(end), strand=-1 if '+' == strand else 1)
                        contig.spawn_feature('retrotransposon', location, 'ltr_harvest')
        return retrotransposon_count
