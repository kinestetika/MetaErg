import re

from metaerg import context
from metaerg.datatypes import fasta
from metaerg.datatypes import sqlite


def _run_programs(genome_name, contig_dict, db_connection, result_files):
    fasta_file = context.spawn_file('masked', genome_name)
    fasta.write_contigs_to_fasta(contig_dict, fasta_file, db_connection, genome_name,
                                      mask_targets=fasta.ALL_MASK_TARGETS)
    context.run_external(f'aragorn -l -t -gc{context.TRANSLATION_TABLE} {fasta_file} -w -o {result_files[0]}')

def _read_results(genome_name, contig_dict, db_connection, result_files) -> int:
    count = 0
    current_contig = None
    coord_regexp = re.compile(r'(c*)\[(\d+),(\d+)]')
    with open(result_files[0]) as aragorn_handle:
        for line in aragorn_handle:
            words = line.strip().split()
            match words:
                case ['>end']:
                    break
                case [contig_name] if contig_name.startswith('>'):
                    current_contig = contig_dict[contig_name[1:]]
                case [_, trna, coordinates, _, codon]:
                    if coord_match := coord_regexp.fullmatch(coordinates):
                        strand = -1 if 'c' == coord_match.group(1) else 1
                        start = max(0, int(coord_match.group(2)) - 1)
                        end = min(len(current_contig['seq']), int(coord_match.group(3)))
                        seq = current_contig['seq'][start:end]
                        if strand < 0:
                            seq = fasta.reverse_complement(seq)
                        feature = sqlite.Feature(genome = genome_name,
                                   contig = current_contig['id'],
                                   start = start,
                                   end = end,
                                   strand = strand,
                                   type = 'tRNA',
                                   inference = 'aragorn',
                                   nt_seq = seq,
                                   descr = f'{trna}-{codon}')
                        sqlite.add_new_feature_to_db(db_connection, feature)
                        count += 1
    return count


@context.register_annotator
def run_and_read_aragorn():
    return ({'pipeline_position': 11,
             'annotator_key': 'aragorn',
             'purpose': 'tRNA prediction with aragorn',
             'programs': ('aragorn',),
             'result_files': ("aragorn",),
             'run': _run_programs,
             'read': _read_results})

