from metaerg.run_and_read.data_model import MetaergGenome, MetaergSeqRecord, FeatureType
from metaerg.run_and_read.context import register, spawn_file, run_external


def _run_programs(genome:MetaergGenome, result_files):
    """Executes the helper programs to complete the analysis"""
    fasta_file = genome.write_fasta_files(spawn_file('masked', genome.id), masked=True)
    run_external(f'minced -gffFull {fasta_file} {result_files[0]}')


def _read_results(genome:MetaergGenome, result_files) -> int:
    """Should parse the result files and return the # of positives"""
    crispr_region_count = 0
    with open(result_files[0]) as crispr_handle:  # this file is in gff format
        for line in crispr_handle:
            words = line.strip().split('\t')
            match words:
                case [str(word), *_] if word.startswith('#'):
                    continue
                case [_, _, 'repeat_region', _, _, _, _, _, _]:
                    crispr_region_count += 1
                case [contig_name, _, 'repeat_unit', start, end, _, strand, _, _]:
                    contig: MetaergSeqRecord = genome.contigs[contig_name]
                    contig.spawn_feature(int(start) - 1, int(end), 1 if '+' == strand else -1,
                                         FeatureType.crispr_repeat, inference='minced')
    return crispr_region_count


@register
def run_and_read_minced():
    return ({'pipeline_position': 1,
             'purpose': 'CRISPR prediction with minced',
             'programs': ('minced',),
             'result_files': ("minced",),
             'run': _run_programs,
             'read': _read_results})

