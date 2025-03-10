from metaerg import context
from metaerg.datatypes import fasta, gff, sqlite

ANNOTATOR_KEY = 'minced'

def _run_programs(genome, contig_dict, db_connection, result_files):
    """Executes the helper programs to complete the analysis"""
    fasta_file = context.spawn_file('masked', genome.name)
    fasta.write_contigs_to_fasta(contig_dict, fasta_file, db_connection, genome.name,
                                 mask_targets=fasta.ALL_MASK_TARGETS)
    context.run_external(f'minced -gffFull {fasta_file} {result_files[0]}')


def _read_results(genome, contig_dict, db_connection, result_files) -> int:
    """Should parse the result files and return the # of positives"""
    with gff.GffParser(result_files[0], contig_dict, inference=ANNOTATOR_KEY,
                       target_feature_type_dict={'repeat_unit': 'CRISPR',
                                                 'repeat_region': 'repeat_region'}) as gff_parser:
        count = 0
        for feature in gff_parser:
            feature.genome = genome.name
            sqlite.add_new_feature_to_db(db_connection, feature)
            count += 1
    return count


@context.register_annotator  # (to enable minced, uncomment and add to __init__
def run_and_read_minced():
    return ({'pipeline_position': 1,
             'annotator_key': ANNOTATOR_KEY,
             'purpose': 'CRISPR prediction with minced',
             'programs': ('minced',),
             'result_files': ("minced",),
             'run': _run_programs,
             'read': _read_results})

