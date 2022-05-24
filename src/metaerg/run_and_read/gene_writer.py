from metaerg.run_and_read.context import register_annotator
from metaerg.run_and_read.data_model import FeatureType, MetaergGenome


def _run_programs(genome:MetaergGenome, result_files):
    genome.write_fasta_files(result_files[0], target=FeatureType.CDS)
    genome.write_fasta_files(result_files[1], target=(FeatureType.rRNA, FeatureType.ncRNA))


def _read_results(genome:MetaergGenome, result_files) -> int:
    return 0


@register_annotator
def run_and_read_trf():
    return ({'pipeline_position': 66,
             'purpose': 'rrite files for proteins and rna genes',
             'programs': (),
             'result_files': ('cds.faa', 'rna.nt'),
             'run': _run_programs,
             'read': _read_results})
