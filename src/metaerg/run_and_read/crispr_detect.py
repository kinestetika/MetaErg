from metaerg import context
from metaerg.datatypes import fasta, gff, sqlite


def _run_programs(genome, contig_dict, db_connection, result_files):
    """Executes the helper programs to complete the analysis"""
    fasta_file = context.spawn_file('masked', genome.name)
    fasta.write_contigs_to_fasta(contig_dict, fasta_file, db_connection, genome.name,
                                 mask_targets=fasta.ALL_MASK_TARGETS)
    context.run_external(f'CRISPRDetect -f {fasta_file} -o {result_files[0][:-4]} -array_quality_score_cutoff 3 -minimum_word_repeatation 3')


def _read_results(genome, contig_dict, db_connection, result_files) -> int:
    """Should parse the result files and return the # of positives"""
    with gff.GffParser(result_files[0], contig_dict, inference='CRISPRDetect') as gff_parser:
        count = 0
        for feature in gff_parser:
            feature.genome = genome.name
            sqlite.add_new_feature_to_db(db_connection, feature)
            count += 1
    return count


@context.register_annotator
def run_and_read_crispr_detect():
    return ({'pipeline_position': 2,
             'annotator_key': 'crispr_detect',
             'purpose': 'CRISPR prediction with CRISPRDetect',
             'programs': ('CRISPRDetect',),
             'result_files': ("crispr_detect.gff",),
             'run': _run_programs,
             'read': _read_results})