from metaerg.data_model import Genome, FeatureType
from metaerg import context
from metaerg import bioparsers

def _run_programs(genome:Genome, result_files):
    fasta_file = context.spawn_file('masked', genome.id)
    bioparsers.write_genome_to_fasta_files(genome, fasta_file, mask=True)
    ltr_index_file = context.spawn_file('ltr_index', genome.id)

    context.run_external(f'gt suffixerator -db {fasta_file} -indexname {ltr_index_file} -tis -suf -lcp -des -ssp -sds -dna')
    context.run_external(f'gt ltrharvest -index {ltr_index_file} -gff3 {result_files[0]} -seqids')
    # remove index files
    for file in ltr_index_file.parent.glob(f'{genome.id}.ltr_index*'):
        file.unlink()


def _read_results(genome:Genome, result_files) -> int:
    retrotransposon_count = 0
    with bioparsers.GffParser(result_files[0], genome.contigs, inference='ltr_harvest',
                              target_feature_type_dict={'repeat_region': FeatureType.retrotransposon}) as parser:
        for c, f in parser:
            c.features.append(f)
            retrotransposon_count += 1
    return retrotransposon_count


@context.register_annotator
def run_and_read_ltr_harvest():
    return ({'pipeline_position': 31,
             'purpose': 'retrotransposon prediction with ltrharvest',
             'programs': ('gt',),
             'result_files': ('ltr_harvest',),
             'run': _run_programs,
             'read': _read_results})
