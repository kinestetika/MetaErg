import shutil
import os
import pandas as pd
from pathlib import Path

from metaerg import context
from metaerg.datatypes import gbk


def _run_programs(genome_name, contig_dict, feature_data: pd.DataFrame, result_files):
    """Should execute the helper programs to complete the analysis"""
    gbk_file = context.spawn_file('gbk', genome_name)
    if context.MULTI_MODE:
        gbk_file = Path(gbk_file.parent, gbk_file.name + '.gbk')
    with open(gbk_file, 'w') as handle:
        gbk.gbk_write_genome(handle, contig_dict, feature_data)
    context.run_external(f'antismash --genefinding-tool none --output-dir {result_files[0]} {gbk_file}')


def _read_results(genome_name, contig_dict, feature_data: pd.DataFrame, result_files) -> tuple:
    """Should parse the result files and return the # of positives."""
    antismash_hit_count = 0
    for f in sorted(result_files[0].glob('*region*.gbk')):
        with gbk.GbkFeatureParser(f) as reader:
            antismash_region_name = ''
            antismash_region_number = 0
            for f_as_dict in reader:
                antismash_hit_count += 1
                if 'region' == f_as_dict['type']:
                    antismash_region_name = '{} {}'.format(f_as_dict['product'], f_as_dict['rules'])
                    antismash_region_number = int(f_as_dict['region_number'])
                elif 'CDS' == f_as_dict['type']:
                    antismash_gene_function = f_as_dict.get('gene_functions', '')
                    antismash_gene_category = f_as_dict.get('gene_kind', '')
                    if antismash_region_name:
                        try:
                            if not f_as_dict['locus_tag'] in feature_data.index:
                                # antismash sometimes predicts a new ORF, for example for a lassopeptide
                                # we do not add this to the feature table
                                continue
                        except KeyError:
                            continue
                        feature_data.at[f_as_dict['locus_tag'], 'antismash'] = \
                            ' '.join((f'(region {antismash_region_number})', antismash_region_name,
                                      antismash_gene_function, antismash_gene_category))
                        prev_subsystem_text = feature_data.at[f_as_dict['locus_tag'], 'subsystems']
                        if type(prev_subsystem_text) == str and '[secondary-metabolites]' not in prev_subsystem_text:
                                feature_data.at[f_as_dict['locus_tag'], 'subsystems'] += ' [secondary-metabolites]'
                        elif f_as_dict['locus_tag'] in feature_data.index:
                            print(f"starting: '{f_as_dict['locus_tag']}'")
                            feature_data.at[f_as_dict['locus_tag'], 'subsystems'] = '[secondary-metabolites]'
    if not antismash_hit_count:
        result_files[0].mkdir(exist_ok=True)  # to prevent re-doing fruitless searches
    return feature_data, antismash_hit_count


@context.register_annotator
def run_and_read_antismash():
    return ({'pipeline_position': 91,
             'purpose': 'prediction of secondary metabolite genes with antismash',
             'programs': ('antismash',),
             'result_files': ('antismash',),
             'run': _run_programs,
             'read': _read_results})


@context.register_html_writer
def write_html(genome_name, feature_data: pd.DataFrame, genome_properties:dict, dir):
    """need to copy the antismash result dir to the metaerg html dir."""
    dir.mkdir(exist_ok=True, parents=True)
    antismash_result_dir = context.spawn_file('antismash', genome_name)
    antismash_html_parent = Path(dir, genome_name, 'antismash')
    if antismash_html_parent.exists():
        shutil.rmtree(antismash_html_parent)
    if antismash_result_dir.exists():
        shutil.copytree(antismash_result_dir, antismash_html_parent)


@context.register_database_installer
def format_blast_databases():
    if 'A' not in context.TASKS:
        return
    antismash_database_path = context.DATABASE_DIR / 'antismash'
    context.log(f'Installing antismash database at {antismash_database_path}')
    antismash_database_path.mkdir(parents=True, exist_ok=True)
    os.system(f'download-antismash-databases --database-dir {antismash_database_path}')
