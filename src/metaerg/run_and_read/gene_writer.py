import pandas as pd
from metaerg import context
from datatypes import fasta

RNA_TARGETS = set("rRNA tRNA tmRNA ncRNA retrotransposon".split())

def _run_programs(genome_name, contig_dict, feature_data: pd.DataFrame, result_files):
    feature_data = feature_data.sort_values(by='start')
    for i in range(len(feature_data.index)):
        feature_data.at[i]['id'] = context.DELIMITER.join((genome_name, feature_data.at[i]['contig'], f'{i:05d}'))
    feature_data.set_index('id')

    context.log(f'({genome_name}) Now writing CDS to fasta file...')
    bioparsers.write_features_to_fasta(feature_data, result_files[0], targets=set('CDS'))
    context.log(f'({genome_name}) Now writing RNA genes and features to fasta file...')
    bioparsers.write_features_to_fasta(feature_data, result_files[1], targets=RNA_TARGETS)

def _read_results(genome_name, contig_dict, feature_data: pd.DataFrame, result_files) -> tuple:
    return feature_data, len(feature_data.index)


@context.register_annotator
def run_and_read_trf():
    return ({'pipeline_position': 66,
             'purpose': 'feature ID generation',
             'programs': (),
             'result_files': ('cds.faa', 'rna.nt'),
             'run': _run_programs,
             'read': _read_results})
