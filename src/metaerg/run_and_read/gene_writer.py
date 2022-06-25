from metaerg import context
from metaerg import bioparsers
from metaerg.data_model import FeatureType, Genome, RNA_FEATURES


def _run_programs(genome:Genome, result_files):
    genome.generate_feature_ids()
    context.log(f'({genome.id}) Now writing CDS to fasta file...')
    bioparsers.write_genome_to_fasta_files(genome, result_files[0], targets=(FeatureType.CDS,))
    context.log(f'({genome.id}) Now writing RNA genes and features to fasta file...')
    bioparsers.write_genome_to_fasta_files(genome, result_files[1], targets=RNA_FEATURES)

def _read_results(genome:Genome, result_files) -> int:
    return sum(len(c.features) for c in genome.contigs.values())


@context.register_annotator
def run_and_read_trf():
    return ({'pipeline_position': 66,
             'purpose': 'feature ID generation',
             'programs': (),
             'result_files': ('cds.faa', 'rna.nt'),
             'run': _run_programs,
             'read': _read_results})
