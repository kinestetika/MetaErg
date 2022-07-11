import pandas as pd
from metaerg import context
from metaerg.datatypes import fasta, gff


def _run_programs(genome_name, contig_dict, feature_data: pd.DataFrame, result_files):
    fasta_file = context.spawn_file('masked', genome_name)
    fasta.write_contigs_to_fasta(contig_dict, fasta_file, feature_data, genome_name,
                                 mask_targets=fasta.ALL_MASK_TARGETS)
    ltr_index_file = context.spawn_file('ltr_index', genome_name)

    context.run_external(f'gt suffixerator -db {fasta_file} -indexname {ltr_index_file} -tis -suf -lcp -des -ssp -sds -dna')
    context.run_external(f'gt ltrharvest -index {ltr_index_file} -gff3 {result_files[0]} -seqids')
    # remove index files
    for file in ltr_index_file.parent.glob(f'{genome_name}.ltr_index*'):
        file.unlink()


def _read_results(genome_name, contig_dict, feature_data: pd.DataFrame, result_files) -> tuple:
    new_features = []
    with gff.GffParser(result_files[0], contig_dict, inference='ltr_harvest',
                       target_feature_type_dict={'repeat_region': 'retrotransposon'}) as parser:
        for feature in parser:
            new_features.append(feature)
    feature_data = pd.concat([feature_data, pd.DataFrame(new_features)], ignore_index=True)
    return feature_data, len(new_features)


@context.register_annotator
def run_and_read_ltr_harvest():
    return ({'pipeline_position': 31,
             'purpose': 'retrotransposon prediction with ltrharvest',
             'programs': ('gt',),
             'result_files': ('ltr_harvest',),
             'run': _run_programs,
             'read': _read_results})
