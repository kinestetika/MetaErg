import shutil
import os
from pathlib import Path

from metaerg import context
from metaerg.datatypes import gff
from metaerg.datatypes import sqlite
from metaerg.datatypes.functional_genes import FunctionalGene

def _run_programs(genome, contig_dict, db_connection, result_files):
    gff_file = context.spawn_file('gff', genome.name, extension='gff')
    with open(gff_file, 'w') as handle:
        gff.gff_write_genome(handle, contig_dict, db_connection)
    cds_aa_file = context.spawn_file('cds.faa', genome.name)
    padloc_database_path = context.DATABASE_DIR / 'padloc'
    # crispr_detect_gff_path = context.spawn_file('crispr_detect.gff', genome.name)
    result_files[0].mkdir(exist_ok=True, parents=True)
    context.run_external(f'padloc --cpu {context.CPUS_PER_GENOME} --faa {cds_aa_file} --gff {gff_file} '
                         f'--outdir {result_files[0]} --force ')
                         #f'--crispr {crispr_detect_gff_path}')


def _read_results(genome, contig_dict, db_connection, result_files) -> int:
    # read functions from padloc database
    cds_aa_file = context.spawn_file('cds.faa', genome.name)
    padloc_result_file = result_files[0] / (cds_aa_file.stem + '_padloc.gff')
    if not padloc_result_file.exists() or not padloc_result_file.stat().st_size:
        context.log('({}) Missing expected result file {}; This happens when padloc '
                    'predicts no defenese genes.', (genome.name, padloc_result_file))
        return 0
    padloc_features = []
    padloc_feature_systems = {}
    with gff.GffParser(padloc_result_file) as gff_parser:
        for pf in gff_parser:
            # padloc does not report the origal gene id, so we find the genes by their start position
            found_feature = False
            for f in sqlite.read_all_features(db_connection, additional_sql = f'start = {pf.start}'):
                f.subsystems.append(FunctionalGene(f'Defense ({pf.type})', pf.id, 1))
                padloc_features.append(f)
                padloc_feature_systems[f] = pf.type  # we do not have easy access to the Functional Gene
                found_feature = True
                break
            if not found_feature:
                raise Exception(f'Found padloc result for unknown feature at position {pf.start}, '
                                f'may need to rerun metaerg with --force')

    # manage clusters
    region_features = []
    current_region_feature = None
    padloc_features.sort(key=lambda x: (x.contig, padloc_feature_systems[x], x.start))
    # first we simply create a region for cds that are close and have the same padloc "type"
    for i in range(1,len(padloc_features)):
        f1 = padloc_features[i-1]
        f2 = padloc_features[i]
        if f1.contig != f2.contig or abs(f2.start - f1.start) > 5000 \
                or padloc_feature_systems[f1] != padloc_feature_systems[f2]:
            current_region_feature = None
            continue
        if current_region_feature:
            f2.parent = current_region_feature.id
            current_region_feature.end = min(int(f2.end)+1, len(contig_dict[f1.contig]['seq']))
        else:
            current_region_feature = sqlite.Feature(genome=genome,
                                                    contig=f1.contig,
                                                    start=max(f1.start - 1, 0),
                                                    end=min(int(f2.end)+1, len(contig_dict[f1.contig]['seq'])),
                                                    strand=1,
                                                    type='region',
                                                    inference='padloc',
                                                    descr=padloc_feature_systems[f1],
                                                    id=f'padloc_region_{max(f1.start - 1, 0)}')
            f1.parent = current_region_feature.id
            f2.parent = current_region_feature.id
            region_features.append(current_region_feature)
        # then we merge identical regions, while updating their type and children
    removed_region_features = set()
    for i in range(len(region_features)):
        for j in range(i):
            r1 = region_features[i]
            r2 = region_features[j]
            if r1.start == r2.start and r1.end == r2.end:
                r1.descr = f'{r1.descr}; {r2.descr}'
                removed_region_features.add(r2.id)
                for f in padloc_features:
                    if r2.id == f.parent and f.start >= r1.start and f.end <= r1.end:
                        f.parent = r1.id
    # update features in db
    for f in padloc_features:
        sqlite.update_feature_in_db(db_connection, f)
    for r in region_features:
        if not r.id in removed_region_features:
            sqlite.update_feature_in_db(db_connection, r)
    return len(padloc_features)


@context.register_annotator
def run_and_read_antismash():
    return ({'pipeline_position': 92,
             'annotator_key': 'padloc',
             'purpose': 'prediction of antiviral defense mechanisms',
             'programs': ('padloc',),
             'result_files': ('padloc',),
             'run': _run_programs,
             'read': _read_results})


@context.register_database_installer
def format_padloc_databases():
    if 'D' not in context.DATABASE_TASKS:
        return
    context.log('Now installing padloc database...')
    padloc_db = context.PADLOC_DATABASE if context.PADLOC_DATABASE else context.DATABASE_DIR / 'padloc'
    if not padloc_db.exists():
        context.log(f'Creating padloc database dir "{padloc_db}"')
        padloc_db.mkdir()
    else:
        context.log(f'Found existing padloc database dir "{padloc_db}"')
    context.run_external(f'padloc --db-update')
    # the following is done while installing the program 'padloc' in installation.py:
    # ln -s /bio/data/databases/metaerg/padloc /bio/data/metaerg-test-install/padloc/data
