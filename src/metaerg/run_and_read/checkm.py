import ast
from pathlib import Path

from metaerg import context


def _run_programs(genome, contig_dict, db_connection, result_files):
    pass

def _read_results(genome, contig_dict, db_connection, result_files) -> int:
    # checkm_results = {}
    if not context.CHECKM_DIR:
        context.log(f'({genome.name}) No dir with checkm results provided to estimate genome completeness.')
        return 0
    checkm_result_file = Path(context.CHECKM_DIR, 'storage', 'bin_stats_ext.tsv')
    # First try checkm:
    try:
        with open(checkm_result_file) as handle:
            for line in handle:
                words = line.split('\t')
                if words[0] == genome.input_fasta_file or words[0] == genome.name or words[0] == f'{genome.name}.{context.FILE_EXTENSION}':
                    checkm_result = ast.literal_eval(words[1])
                    genome.fraction_complete = float(checkm_result["Completeness"]) / 100
                    genome.fraction_contaminated = float(checkm_result["Contamination"]) / 100
                    return 1
        context.log(f'({genome.name}) No results found for {genome.name} ({genome.input_fasta_file}) in {checkm_result_file}.')
    except:
        # Now try checkm2:
        checkm_result_file = Path(context.CHECKM_DIR, 'quality_report.tsv')
        try:
            # First try checkm:
            with open(checkm_result_file) as handle:
                for line in handle:
                    w = line.split('\t')
                    if w[0] == genome.input_fasta_file or w[0] == genome.name or w[0] == f'{genome.name}.{context.FILE_EXTENSION}':
                        genome.fraction_complete = float(w[1]) / 100
                        genome.fraction_contaminated = float(w[2]) / 100
                        return 1
            context.log(f'({genome.name}) No results found for {genome.name} ({genome.input_fasta_file}) in {checkm_result_file}.')
        except:
            context.log(f'({genome.name}) No results found for checkm or checkm2 at {context.CHECKM_DIR}.')


@context.register_annotator
def run_and_read_trf():
    return ({'pipeline_position': 131,
             'annotator_key': 'checkm',
             'purpose': 'estimate genome completeness and contamination',
             'programs': (),
             'result_files': (),
             'run': _run_programs,
             'read': _read_results})

