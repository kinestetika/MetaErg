import re

from metaerg import context
from metaerg.datatypes import fasta
from metaerg.datatypes import sqlite


def _run_programs(genome, contig_dict, db_connection, result_files):
    fasta_file = context.spawn_file('masked', genome.name)
    fasta.write_contigs_to_fasta(contig_dict, fasta_file, db_connection, genome.name,
                                      mask_targets=fasta.ALL_MASK_TARGETS)
    if not context.TRANSLATION_TABLE:
        set_translation_table(genome, 0, 0)
        context.run_external(f'aragorn -l -m -t -gcstd {fasta_file} -w -o {result_files[0]}')
    else:
        ## need to determine a suitable genetic code... for that we use prodigal
        mean_protein_length_of_default_code = 0
        for table_id in context.TRANSLATION_TABLE:
            context.run_external(f'prodigal -g {table_id} -m -f gff -q -i {fasta_file} -a {result_files[0]}')
            count = 0
            total_length = 0
            with fasta.FastaParser(result_files[0], cleanup_seq=False) as fasta_reader:
                for seq_rec in fasta_reader:
                    count += 1
                    total_length += len(seq_rec['seq'])
            if not count:
                context.log(f'({genome.name}) No proteins found with translation table {table_id}.')
                continue
            mean_protein_length = total_length / count
            if not mean_protein_length_of_default_code:
                mean_protein_length_of_default_code = mean_protein_length
            if mean_protein_length > 200:
                set_translation_table(genome, table_id, mean_protein_length)
                context.run_external(f'aragorn -l -t -gc{genome.genetic_code} {fasta_file} -w -o {result_files[0]}')
                return
            else:
                context.log(f'({genome.name}) Coding sequence prediction with translation table {table_id} failed, mean length {mean_protein_length:.1f} aa.')
        # if no suitable code was found, we default to the first value in context.TRANSLATION_TABLE
        set_translation_table(genome, context.TRANSLATION_TABLE[0], mean_protein_length_of_default_code)
        context.run_external(f'aragorn -l -t -gc{genome.genetic_code} {fasta_file} -w -o {result_files[0]}')


def _read_results(genome, contig_dict, db_connection, result_files) -> int:
    count = 0
    ln = 0
    current_contig = None
    coord_regexp = re.compile(r'(c*)\[(\d+),(\d+)]')
    with open(result_files[0]) as aragorn_handle:
        for line in aragorn_handle:
            words = line.strip().split()
            ln += 1
            match words:
                case ['>end']:
                    break
                case [contig_name] if contig_name.startswith('>'):
                    if space_index := contig_name.find(' ') > 0:
                        contig_name = contig_name[:space_index]
                    current_contig = contig_dict[contig_name[1:]]
                case [_, trna_type, coordinates, _, codon]:
                    if coord_match := coord_regexp.fullmatch(coordinates):
                        if current_contig:
                            if not trna_type.startswith('tRNA') and not trna_type.startswith('tmRNA'):
                                continue
                            strand = -1 if 'c' == coord_match.group(1) else 1
                            start = max(0, int(coord_match.group(2)) - 1)
                            end = min(len(current_contig['seq']), int(coord_match.group(3)))
                            seq = current_contig['seq'][start:end]
                            if strand < 0:
                                seq = fasta.reverse_complement(seq)
                            feature = sqlite.Feature(genome = genome.name,
                                       contig = current_contig['id'],
                                       start = start,
                                       end = end,
                                       strand = strand,
                                       type = 'tRNA' if trna_type.startswith('tRNA') else 'tmRNA',
                                       inference = 'aragorn',
                                       nt_seq = seq,
                                       descr = f'{trna_type}-{codon}')
                            sqlite.add_new_feature_to_db(db_connection, feature)
                            count += 1
                        else:
                            context.log(f'Warning: Unexpected tRNA coordinates in "{line}" '
                                    f'at line {ln - 1} of {result_files[0]}')
    genetic_code_file = context.spawn_file('genetic_code', genome.name)
    genome.genetic_code = int(genetic_code_file.read_text())
    return count


def set_translation_table(genome, table_id, mean_protein_length):
    context.log(
        f'({genome.name}) Coding sequence prediction with translation table {table_id}, mean length {mean_protein_length:.1f} aa.')
    genetic_code_file = context.spawn_file('genetic_code', genome.name)
    with open(genetic_code_file, 'w') as handle:
        handle.write(str(table_id))
    genome.genetic_code = table_id


@context.register_annotator
def run_and_read_aragorn():
    return ({'pipeline_position': 11,
             'annotator_key': 'aragorn',
             'purpose': 'tRNA prediction with aragorn',
             'programs': ('aragorn',),
             'result_files': ("aragorn",),
             'run': _run_programs,
             'read': _read_results})

