
from metaerg import context
from metaerg.datatypes import fasta
from metaerg.datatypes import sqlite


ANNOTATOR_KEY = 'repeater'


def _run_programs(genome, contig_dict, db_connection, result_files):
    fasta_file = context.spawn_file('repeater', genome.name)
    fasta.write_contigs_to_fasta(contig_dict, fasta_file, db_connection, genome.name,
                                 mask_targets=fasta.ALL_MASK_TARGETS)
    context.run_external(f'repeater {fasta_file}')
    fasta_file.unlink(missing_ok=True)



def _read_results(genome, contig_dict, db_connection, result_files) -> int:
    # results look incomprehensible with genome-long repeats...
    pass


#@context.register_annotator
def run_and_read_repeatmasker():
    return ({'pipeline_position': 51,
             'annotator_key': ANNOTATOR_KEY,
             'purpose': 'repeat prediction with Repeater',
             'programs': ('repeater',),
             'result_files': ('repeater.gff',),
             'run': _run_programs,
             'read': _read_results})