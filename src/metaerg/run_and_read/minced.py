from metaerg.data_model import Genome, FeatureType
from metaerg import context
from metaerg import bioparsers

def _run_programs(genome:Genome, result_files):
    """Executes the helper programs to complete the analysis"""
    fasta_file = context.spawn_file('masked', genome.id)
    bioparsers.write_genome_to_fasta_files(genome, fasta_file, mask=True)
    context.run_external(f'minced -gffFull {fasta_file} {result_files[0]}')


def _read_results(genome:Genome, result_files) -> int:
    """Should parse the result files and return the # of positives"""
    crispr_count = 0
    with bioparsers.GffParser(result_files[0], genome.contigs, inference='minced',
                              target_feature_type_dict={'repeat_unit': FeatureType.crispr_repeat}) as gff_parser:
        for contig, feature in gff_parser:
            contig.features.append(feature)
            crispr_count += 1
    return crispr_count


@context.register_annotator
def run_and_read_minced():
    return ({'pipeline_position': 1,
             'purpose': 'CRISPR prediction with minced',
             'programs': ('minced',),
             'result_files': ("minced",),
             'run': _run_programs,
             'read': _read_results})

