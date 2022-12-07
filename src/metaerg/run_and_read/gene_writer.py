import pandas as pd
from metaerg import context
from metaerg.datatypes import fasta

def _run_programs(genome_name, contig_dict, feature_data: pd.DataFrame, result_files):
    pass


def _read_results(genome_name, contig_dict, feature_data: pd.DataFrame, result_files) -> tuple:
    feature_data = feature_data.sort_values(by=['contig', 'start'])
    j=0
    for i in feature_data.index:
        if context.RENAME_CONTIGS:
            # contigs already contain genome name
            feature_data.at[i, 'id'] = context.DELIMITER.join((feature_data.at[i, 'contig'], f'{j:05d}'))
        else:
            feature_data.at[i, 'id'] = context.DELIMITER.join((genome_name, feature_data.at[i, 'contig'], f'{j:05d}'))
        j += 1
    feature_data = feature_data.set_index('id', drop=False)
    feature_data = feature_data.fillna({'tmh': 0})
    feature_data = feature_data.fillna('')

    cds_count = len(feature_data[feature_data['type'] == 'CDS'])
    rna_count = len(feature_data[feature_data['type'].isin(context.RNA_TARGETS)])

    cds_file = context.spawn_file('cds.faa', genome_name)
    context.log(f'({genome_name}) Now writing {cds_count} proteins to fasta at {cds_file}...')
    fasta.write_features_to_fasta(feature_data, 'aa', cds_file, targets=('CDS',))
    rna_file = context.spawn_file('rna.fna', genome_name)
    context.log(f'({genome_name}) Now writing {rna_count} RNA genes and features to fasta at {rna_file}...')
    fasta.write_features_to_fasta(feature_data, 'nt', rna_file, targets=context.RNA_TARGETS)
    return feature_data, len(feature_data.index)


@context.register_annotator
def run_and_read_trf():
    return ({'pipeline_position': 66,
             'purpose': 'feature ID generation',
             'programs': (),
             'result_files': (),
             'run': _run_programs,
             'read': _read_results})
