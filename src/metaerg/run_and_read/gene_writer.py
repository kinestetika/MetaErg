from metaerg import context
from metaerg import bioparsers
from metaerg.data_model import FeatureType, MetaergGenome, RNA_FEATURES


def _run_programs(genome:MetaergGenome, result_files):
    genome.generate_feature_ids()
    bioparsers.write_genome_fasta_files(genome, result_files[0], target=FeatureType.CDS)
    bioparsers.write_genome_fasta_files(genome, result_files[1], target=RNA_FEATURES)


def _read_results(genome:MetaergGenome, result_files) -> int:
    return 0


@context.register_annotator
def run_and_read_trf():
    return ({'pipeline_position': 66,
             'purpose': 'generate feature IDs, write files for proteins and rna genes',
             'programs': (),
             'result_files': ('cds.faa', 'rna.nt'),
             'run': _run_programs,
             'read': _read_results})
