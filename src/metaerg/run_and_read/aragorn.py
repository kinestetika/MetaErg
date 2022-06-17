import re

import bioparsers
from metaerg.data_model import FeatureType, Genome, SeqFeature
from metaerg import context
from metaerg import bioparsers


def _run_programs(genome:Genome, result_files):
    fasta_file = genome, context.spawn_file('masked', genome.id)
    bioparsers.write_genome_to_fasta_files(genome, fasta_file, mask=True)
    context.run_external(f'aragorn -l -t -gc{genome.translation_table} {fasta_file} -w -o {result_files[0]}')

def _read_results(genome:Genome, result_files) -> int:
    trna_count = 0
    current_contig = None
    coord_regexp = re.compile(r'(c*)\[(\d+),(\d+)]')
    with open(result_files[0]) as aragorn_handle:
        for line in aragorn_handle:
            words = line.strip().split()
            match words:
                case ['>end']:
                    break
                case [contig_name] if contig_name.startswith('>'):
                    current_contig = genome.contigs[contig_name[1:]]
                case [_, trna, coordinates, _, codon]:
                    trna_count += 1
                    coord_match = coord_regexp.fullmatch(coordinates)
                    strand = -1 if 'c' == coord_match.group(1) else 1
                    start = max(0, int(coord_match.group(2)) - 1)
                    end = min(len(current_contig.seq), int(coord_match.group(3)))
                    seq = current_contig.seq[start:end]
                    if strand < 0:
                        seq = bioparsers.reverse_complement(seq)
                    feature = SeqFeature(start, end, strand, FeatureType.tRNA, 'aragorn', seq=seq,
                                         descr=f'{trna}-{codon}')
                    current_contig.features.append(feature)
    return trna_count


@context.register_annotator
def run_and_read_aragorn():
    return ({'pipeline_position': 11,
             'purpose': 'tRNA prediction with aragorn',
             'programs': ('aragorn',),
             'result_files': ("aragorn",),
             'run': _run_programs,
             'read': _read_results})

