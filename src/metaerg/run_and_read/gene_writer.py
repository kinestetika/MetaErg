import pandas as pd
from metaerg import context
from metaerg.datatypes import fasta

RNA_TARGETS = set("rRNA tRNA tmRNA ncRNA retrotransposon".split())


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

    context.log(f'({genome_name}) Now writing proteins to fasta file...')
    cds_file = context.spawn_file('cds.faa', genome_name)
    fasta.write_features_to_fasta(feature_data, cds_file, targets=('CDS',))
    context.log(f'({genome_name}) Now writing RNA genes and features to fasta file...')
    rna_file = context.spawn_file('rna.fna', genome_name)
    fasta.write_features_to_fasta(feature_data, rna_file, targets=RNA_TARGETS)
    return feature_data, len(feature_data.index)


@context.register_annotator
def run_and_read_trf():
    return ({'pipeline_position': 66,
             'purpose': 'feature ID generation',
             'programs': (),
             'result_files': (),
             'run': _run_programs,
             'read': _read_results})
