import shutil
import os
from pathlib import Path

from metaerg import context
from metaerg.datatypes import gbk
from metaerg.datatypes import sqlite
from metaerg.datatypes import functional_genes


def _run_programs(genome, contig_dict, db_connection, result_files):
    """Should execute the helper programs to complete the analysis"""
    gbk_file = context.spawn_file('gbk', genome.name, extension='gbk')
    with open(gbk_file, 'w') as handle:
        gbk.gbk_write_genome(handle, contig_dict, db_connection)
    if result_files[0].exists():
        shutil.rmtree(result_files[0])
    context.run_external(f'antismash --genefinding-tool none --output-dir {result_files[0]} {gbk_file}')


def _read_results(genome, contig_dict, db_connection, result_files) -> int:
    """Should parse the result files and return the # of positives."""
    antismash_hit_count = 0
    for f in sorted(result_files[0].glob('*region*.gbk')):
        contig = ''
        start = 0
        end = 0
        with open(f) as reader:
            for line in reader:
                line = line.strip()
                if line.startswith('LOCUS'):
                    contig = line.split()[1]
                elif line.startswith('Orig. start'):
                    start = line.split()[-1]
                elif line.startswith('Orig. end'):
                    end = line.split()[-1]
                    break
        if not contig or not end:
            context.log('Warning: could not locate antimash results.')
            continue
        region_id = f'antismash_region_{start}'
        region_feature = sqlite.Feature(genome=genome,
                                        contig=contig,
                                        start=max(int(start)-1, 0),
                                        end=int(end),
                                        strand=1,
                                        type='region',
                                        inference='antismash',
                                        id=region_id)
        region_feature.subsystems.append(functional_genes.SECONDARY_METABOLITE_GENE)
        antismash_genes = []
        with gbk.GbkFeatureParser(f) as reader:
            for f_as_dict in reader:
                if 'region' == f_as_dict['type']:
                    region_feature.descr = '{} {}'.format(f_as_dict['product'], f_as_dict['rules'])
                    sqlite.add_new_feature_to_db(db_connection, region_feature)
                elif 'CDS' == f_as_dict['type']:
                    antismash_genes.append(f_as_dict)
            for g in antismash_genes:
                antismash_gene_function = g.get('gene_functions', '')
                antismash_gene_category = g.get('gene_kind', '')
                if feature := sqlite.read_feature_by_id(db_connection, g['locus_tag']):
                    feature.parent = region_id
                    feature.notes += ' ' + ' '.join((region_feature.descr, antismash_gene_function,
                                                     antismash_gene_category))
                    feature.subsystems.append(functional_genes.SECONDARY_METABOLITE_GENE)
                    sqlite.update_feature_in_db(db_connection, feature)
                else:
                    # antismash sometimes predicts a new ORF, for example for a lassopeptide
                    new_feature = sqlite.Feature(id=f'antismash_{antismash_hit_count}',
                                                 genome=genome,
                                                 contig=contig,
                                                 start=g['start'],
                                                 end=g['end'],
                                                 strand=g['strand'],
                                                 type=g['type'],
                                                 inference='antismash',
                                                 parent = region_id,
                                                 descr= ' '.join((region_feature.descr, antismash_gene_function,
                                                                  antismash_gene_category)),
                                                 aa_seq=g['translation'])
                    #print('translation', new_feature.aa_seq)
                    new_feature.subsystems.append(functional_genes.SECONDARY_METABOLITE_GENE)
                    sqlite.add_new_feature_to_db(db_connection, new_feature)
                    # we do try to add this to the feature table
                antismash_hit_count += 1
    return antismash_hit_count


@context.register_annotator
def run_and_read_antismash():
    return ({'pipeline_position': 91,
             'annotator_key': 'antismash',
             'purpose': 'prediction of secondary metabolite genes with antismash',
             'programs': ('antismash',),
             'result_files': ('antismash',),
             'run': _run_programs,
             'read': _read_results})


@context.register_html_writer
def write_html(genome, db_connection, dir):
    """need to copy the antismash result dir to the metaerg html dir."""
    dir.mkdir(exist_ok=True, parents=True)
    antismash_result_dir = context.spawn_file('antismash', genome.name)
    antismash_html_parent = Path(dir, genome.name, 'antismash')
    if antismash_html_parent.exists():
        shutil.rmtree(antismash_html_parent)
    if antismash_result_dir.exists():
        shutil.copytree(antismash_result_dir, antismash_html_parent)


@context.register_database_installer
def format_antismash_databases():
    if 'A' not in context.DATABASE_TASKS:
        return
    antismash_database_path = context.DATABASE_DIR / 'antismash'
    context.log(f'Installing antismash database at {antismash_database_path}')
    antismash_database_path.mkdir(parents=True, exist_ok=True)
    os.system(f'download-antismash-databases --database-dir {antismash_database_path}')
