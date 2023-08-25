from metaerg import context
from metaerg.datatypes import fasta, gff, sqlite


def _run_programs(genome, contig_dict, db_connection, result_files):
    """Executes the helper programs to complete the analysis"""
    fasta_file = context.spawn_file('masked', genome.name)
    fasta.write_contigs_to_fasta(contig_dict, fasta_file, db_connection, genome.name,
                                 mask_targets=fasta.ALL_MASK_TARGETS)
    result_file_base_name = result_files[0].name[:-4]
    context.run_external(f'CRISPRDetect.pl -f {fasta_file} -o {context.TEMP_DIR / result_file_base_name} -array_quality_score_cutoff 3 -minimum_word_repeatation 3')
    # remove unneeded output
    for file in context.TEMP_DIR.glob(f'{result_file_base_name}*'):
        if not file.name.endswith('gff'):
            # context.log(f'removing {file}')
            file.unlink()


def _read_results(genome, contig_dict, db_connection, result_files) -> int:
    """Should parse the result files and return the # of positives"""
    with gff.GffParser(result_files[0], contig_dict) as gff_parser:
        count = 0
        for feature in gff_parser:
            feature.genome = genome.name
            if 'direct_repeat' == feature.type:
                feature.type = 'CRISPR'
                count += 1
            sqlite.add_new_feature_to_db(db_connection, feature)
    return count


@context.register_annotator
def run_and_read_crispr_detect():
    return ({'pipeline_position': 2,
             'annotator_key': 'crispr_detect',
             'purpose': 'CRISPR prediction with CRISPRDetect',
             'programs': ('CRISPRDetect.pl',),
             'result_files': ("crispr_detect.gff",),
             'run': _run_programs,
             'read': _read_results})