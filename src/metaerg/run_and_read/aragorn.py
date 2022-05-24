import re
from metaerg.run_and_read.data_model import FeatureType, MetaergGenome
from metaerg.run_and_read.context import register_annotator, spawn_file, run_external


def _run_programs(genome:MetaergGenome, result_files):
    fasta_file = genome.write_fasta_files(spawn_file('masked', genome.id), masked=True)
    run_external(f'aragorn -l -t -gc{genome.translation_table} {fasta_file} -w -o {result_files[0]}')

def _read_results(genome:MetaergGenome, result_files) -> int:
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
                    end = min(len(current_contig.sequence), int(coord_match.group(3)))
                    f = current_contig.spawn_feature(start, end, strand, FeatureType.tRNA,
                                                     inference='aragorn')
                    f.description = f'{trna}-{codon}'
    return trna_count


@register_annotator
def run_and_read_aragorn():
    return ({'pipeline_position': 11,
             'purpose': 'tRNA prediction with aragorn',
             'programs': ('aragorn',),
             'result_files': ("aragorn",),
             'run': _run_programs,
             'read': _read_results})

