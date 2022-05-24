from metaerg.run_and_read.data_model import MetaergGenome, MetaergSeqRecord, FeatureType
from metaerg.run_and_read.context import register, spawn_file, run_external


def _run_programs(genome:MetaergGenome, result_files):
    fasta_file = genome.write_fasta_files(spawn_file('masked', genome.id), masked=True)
    ltr_index_file = spawn_file('ltr_index', genome.id)

    run_external(f'gt suffixerator -db {fasta_file} -indexname {ltr_index_file} -tis -suf -lcp -des -ssp -sds -dna')
    run_external(f'gt ltrharvest -index {ltr_index_file} -gff3 {result_files[0]} -seqids')
    # remove index files
    for file in ltr_index_file.parent.glob(f'{genome.id}.ltr_index*'):
        file.unlink()


def _read_results(genome:MetaergGenome, result_files) -> int:
    retrotransposon_count = 0
    with open(result_files[0]) as ltr_handle:
        for line in ltr_handle:
            words = line.strip().split('\t')
            match words:
                case [str(word), *_] if word.startswith('#'):
                    continue
                case [contig_name, _, 'repeat_region', start, end, score, strand, frame, _]:
                    retrotransposon_count += 1
                    contig: MetaergSeqRecord = genome.contigs[contig_name]
                    contig.spawn_feature(int(start) - 1, int(end), 1 if '+' == strand else -1,
                                         FeatureType.retrotransposon, inference='ltr_harvest')
    return retrotransposon_count


@register
def run_and_read_ltr_harvest():
    return ({'pipeline_position': 31,
             'purpose': 'retrotransposon prediction with ltrharvest',
             'programs': ('ltr_index', 'ltr_harvest'),
             'result_files': ('ltr_harvest',),
             'run': _run_programs,
             'read': _read_results})
