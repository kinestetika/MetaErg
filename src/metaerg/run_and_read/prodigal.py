from metaerg.data_model import Genome, SeqRecord, FeatureType
from metaerg import context
from metaerg import bioparsers


def _run_programs(genome:Genome, result_files):
    fasta_file, = bioparsers.write_genome_to_fasta_files(genome, context.spawn_file('masked', genome.id), mask=True)
    context.run_external(f'prodigal -g {genome.translation_table} -m -f gff -q -i {fasta_file} -o {result_files[0]}')


def _read_results(genome:Genome, result_files) -> int:
    cds_found = 0
    with open(result_files[0]) as prodigal_handle:
        for line in prodigal_handle:
            words = line.strip().split('\t')
            match words:
                case [str(word), *_] if word.startswith('#'):
                    continue
                case [contig_name, _, _, start, end, _, strand, _, attributes]:
                    cds_found += 1
                    contig: SeqRecord = genome.contigs[contig_name]
                    feature = contig.spawn_feature(int(start) - 1, int(end), 1 if '+' == strand else -1,
                                                   FeatureType.CDS, inference='prodigal')
                    if 'partial=01' in attributes or 'partial=01' in attributes or 'partial=11' in attributes:
                        feature.notes.add('partial protein')
    return cds_found


@context.register_annotator
def run_and_read_prodigal():
    return ({'pipeline_position': 61,
             'purpose': 'coding sequence prediction with prodigal',
             'programs': ('prodigal'),
             'result_files': ('prodigal',),
             'run': _run_programs,
             'read': _read_results})
