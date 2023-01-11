from metaerg import context
from metaerg.datatypes import fasta, gff, sqlite


def _run_programs(genome_name, contig_dict, db_connection, result_files):
    """Executes the helper programs to complete the analysis"""
    fasta_file = context.spawn_file('masked', genome_name)
    fasta.write_contigs_to_fasta(contig_dict, fasta_file, db_connection, genome_name,
                                 mask_targets=fasta.ALL_MASK_TARGETS)
    context.run_external(f'minced -gffFull {fasta_file} {result_files[0]}')


def _read_results(genome_name, contig_dict, db_connection, result_files) -> int:
    """Should parse the result files and return the # of positives"""
    with gff.GffParser(result_files[0], contig_dict, inference='minced',
                       target_feature_type_dict={'repeat_unit': 'crispr_repeat'}) as gff_parser:
        count = 0
        for feature in gff_parser:
            feature.genome = genome_name
            sqlite.add_new_feature_to_db(db_connection, feature)
            count += 1
    return count


@context.register_annotator
def run_and_read_minced():
    return ({'pipeline_position': 1,
             'annotator_key': 'minced',
             'purpose': 'CRISPR prediction with minced',
             'programs': ('minced',),
             'result_files': ("minced",),
             'run': _run_programs,
             'read': _read_results})

