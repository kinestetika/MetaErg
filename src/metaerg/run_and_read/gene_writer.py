import pandas as pd

from metaerg import context
from metaerg.datatypes import fasta
from metaerg.datatypes import sqlite

def _run_programs(genome_name, contig_dict, feature_data: pd.DataFrame, result_files):
    pass


def _read_results(genome_name, contig_dict, db_connection, result_files) -> int:
    j = 0
    for feature in sqlite.read_all_features(db_connection):
        if context.RENAME_CONTIGS:
            # contigs already contain genome name
            feature.id = context.DELIMITER.join((feature.contig, f'{j:05d}'))
        else:
            feature.id = context.DELIMITER.join((genome_name, feature.contig, f'{j:05d}'))
        j += 1
        sqlite.update_feature_in_db(db_connection, feature)

    cds_count = sum(1 for f in sqlite.read_all_features(db_connection, type='CDS'))
    rna_count = sum(1 for f in sqlite.read_all_features(db_connection, type=sqlite.RNA_TARGETS))

    cds_file = context.spawn_file('cds.faa', genome_name)
    context.log(f'({genome_name}) Now writing {cds_count} proteins to fasta at {cds_file}...')
    fasta.write_features_to_fasta(db_connection, 'aa', cds_file, targets=('CDS',))
    rna_file = context.spawn_file('rna.fna', genome_name)
    context.log(f'({genome_name}) Now writing {rna_count} RNA genes and features to fasta at {rna_file}...')
    fasta.write_features_to_fasta(db_connection, 'nt', rna_file, targets=sqlite.RNA_TARGETS)
    return j


@context.register_annotator
def run_and_read_trf():
    return ({'pipeline_position': 66,
             'annotator_key': 'write_genes',
             'purpose': 'feature ID generation',
             'programs': (),
             'result_files': (),
             'run': _run_programs,
             'read': _read_results})
