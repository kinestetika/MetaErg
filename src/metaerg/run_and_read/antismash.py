import shutil
import os
from pathlib import Path

from metaerg import context
from metaerg.datatypes import gbk
from metaerg.datatypes import sqlite
from metaerg.datatypes import functional_genes
from metaerg.run_and_read import gene_writer


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
            context.log(f'({genome.name}) Warning: could not locate antismash results in {f}.')
            continue
        region_feature = sqlite.Feature(genome=genome,
                                        contig=contig,
                                        start=max(int(start)-1, 0),
                                        end=int(end),
                                        strand=1,
                                        type='region',
                                        inference='antismash',
                                        id=f'antismash_region_{start}')
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
            # context.log(f'({genome.name}) Antismash results for {g["locus_tag"]}".')
            antismash_gene_function = g.get('gene_functions', '')
            antismash_gene_category = g.get('gene_kind', '')
            if feature := sqlite.read_feature_by_id(db_connection, g['locus_tag']):
                # context.log(f'({genome.name}) Located feature {g["locus_tag"]}".')
                feature.parent = region_feature.id
                feature.notes += ' ' + ' '.join((region_feature.descr, antismash_gene_function,
                                                 antismash_gene_category))
                feature.subsystems.append(functional_genes.SECONDARY_METABOLITE_GENE)
                sqlite.update_feature_in_db(db_connection, feature)
            else:
                # antismash sometimes predicts a new ORF, for example for a lassopeptide
                # find the closest gene as a reference for naming the new gene
                def create_new_id(existing_id):
                    if existing_id[-1].isalpha():
                        return existing_id[:-1] + chr(ord(existing_id[-1]) + 1)
                    else:
                        return existing_id + 'a'

                closest_gene = None
                closest_distance = 100000
                for f1 in sqlite.read_all_features(db_connection, contig=contig, start=max(g['start']-20000,0), end=g['end']+20000):
                    new_distance = f1.start - g['start']
                    if new_distance < closest_distance:
                        closest_gene = f1
                        closest_distance = new_distance
                if closest_gene:
                    new_id = create_new_id(closest_gene.id)
                else:
                    new_id = gene_writer.create_new_id(genome.name, contig, 0) + 'a'

                context.log(f'({genome.name}) New gene detected by antismash, created id {new_id}".')
                new_feature = sqlite.Feature(id=new_id,
                                             genome=genome.name,
                                             contig=contig,
                                             start=g['start'],
                                             end=g['end'],
                                             strand=g['strand'],
                                             type=g['type'],
                                             inference='antismash',
                                             descr= ' '.join((region_feature.descr, antismash_gene_function,
                                                              antismash_gene_category)),
                                             aa_seq=g['translation'])
                new_feature.parent = region_feature.id
                new_feature.subsystems.append(functional_genes.SECONDARY_METABOLITE_GENE)
                sqlite.add_new_feature_to_db(db_connection, new_feature)
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
    context.log(f'Now installing antismash database...')
    antismash_db = context.ANTISMASH_DATABASE if context.ANTISMASH_DATABASE else context.DATABASE_DIR / 'antismash'
    if not antismash_db.exists():
        context.log(f'Creating antismash database dir "{antismash_db}"')
        antismash_db.mkdir()
    else:
        context.log(f'Found existing antismash database dir "{antismash_db}"')
    context.run_external(f'download-antismash-databases')
    # the following is done while installing the program 'antismash' in installation.py:
    #os.system(f'ln -s {path_to_antismash_db} {antismash_database_python_dir}')
