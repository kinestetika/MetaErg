import shutil
from pathlib import Path

from metaerg import bioparsers
from metaerg import context
from metaerg.data_model import SeqFeature, Genome


def _run_programs(genome:Genome, result_files):
    """Should execute the helper programs to complete the analysis"""
    gbk_file = context.spawn_file('gbk', genome.id)
    with open(gbk_file, 'w') as handle:
        bioparsers.gbk_write_genome(handle, genome)
    context.run_external(f'antismash --genefinding-tool none --output-dir {result_files[0]} {gbk_file}')


def _read_results(genome:Genome, result_files) -> int:
    """Should parse the result files and return the # of positives."""
    antismash_hit_count = 0
    for f in sorted(result_files[0].glob('*region*.gbk')):
        with bioparsers.GbkFeatureParser(f) as reader:
            antismash_region_name = ''
            antismash_region_number = 0
            for f_as_dict in reader:
                antismash_hit_count += 1
                if 'region' == f_as_dict['type']:
                    antismash_region_name = '{} {}'.format(f_as_dict['product'], f_as_dict['rules'])
                    antismash_region_number = int(f_as_dict['region_number'])
                elif 'CDS' == f_as_dict['type']:
                    feature: SeqFeature = genome.get_feature(f_as_dict['locus_tag'])
                    antismash_gene_function = f_as_dict.get('gene_functions', '')
                    antismash_gene_category = f_as_dict.get('gene_kind', '')
                    if antismash_region_name:
                        feature.antismash = ' '.join((f'(region {antismash_region_number})',
                                                     antismash_region_name,
                                                     antismash_gene_function,
                                                     antismash_gene_category))
                        feature.subsystem.add('[secondary-metabolites]')
                        genome.subsystems.subsystems['[secondary-metabolites]'].add_hit(feature.id)
    if not antismash_hit_count:
        result_files[0].mkdir(exist_ok=True)  # to prevent re-doing fruitless searches
    return antismash_hit_count


@context.register_annotator
def run_and_read_antismash():
    return ({'pipeline_position': 91,
             'purpose': 'prediction of secondary metabolite genes with antismash',
             'programs': ('antismash',),
             'result_files': ('antismash',),
             'run': _run_programs,
             'read': _read_results})


@context.register_html_writer
def write_html(genome: Genome, dir):
    """need to copy the antismash result dir to the metaerg html dir."""
    dir.mkdir(exist_ok=True, parents=True)
    antismash_result_dir = context.spawn_file('antismash', genome.id)
    antismash_html_parent = Path(dir, genome.id, 'antismash')
    if antismash_html_parent.exists():
        shutil.rmtree(antismash_html_parent)
    shutil.copytree(antismash_result_dir, antismash_html_parent)

