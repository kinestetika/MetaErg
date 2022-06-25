import pandas as pd
from metaerg import context
from metaerg import bioparsers

def _run_programs(genome_name, contig_dict, feature_data: pd.DataFrame, result_files):
    """Executes the helper programs to complete the analysis"""
    fasta_file = context.spawn_file('masked', genome_name)
    bioparsers.write_contigs_to_fasta(genome_name, contig_dict, feature_data, fasta_file,
                                      mask_targets=bioparsers.ALL_MASK_TARGETS)
    context.run_external(f'minced -gffFull {fasta_file} {result_files[0]}')


def _read_results(genome_name, contig_dict, feature_data: pd.DataFrame, result_files) -> tuple:
    """Should parse the result files and return the # of positives"""
    new_features = []
    with bioparsers.GffParser(result_files[0], contig_dict, inference='minced',
                              target_feature_type_dict={'repeat_unit': 'crispr_repeat'}) as gff_parser:
        for contig, feature in gff_parser:
            new_features.append(feature)
    feature_data = pd.concat([feature_data, pd.DataFrame(new_features)], ignore_index=True)
    return feature_data, len(new_features)


@context.register_annotator
def run_and_read_minced():
    return ({'pipeline_position': 1,
             'purpose': 'CRISPR prediction with minced',
             'programs': ('minced',),
             'result_files': ("minced",),
             'run': _run_programs,
             'read': _read_results})

