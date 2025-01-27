import os
import shutil
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor

from metaerg import context
from metaerg.datatypes import functional_genes
from metaerg.datatypes import fasta
from metaerg.datatypes import sqlite
from metaerg.datatypes.blast import DBentry, TabularBlastParser

CDD = {}
CDD_CLUSTERS = {}
CLUSTER_MAX_CDD_HITS = 5
CLUSTER_MIN_MATCH_SCORE = 0.2

ANNOTATOR_KEY = 'cdd'

def _run_programs(genome, contig_dict, db_connection, result_files):
    cds_aa_file = context.spawn_file('cds.faa', genome.name)
    cdd_database = Path(context.DATABASE_DIR, 'cdd', 'Cdd')
    context.run_external(f'rpsblast -db {cdd_database} -query {cds_aa_file} -out {result_files[0]} -outfmt 6 -evalue 1e-7'
                         f' -num_threads {context.CPUS_PER_GENOME}')


def _read_results(genome, contig_dict, db_connection, result_files) -> int:
    cdd_result_count = 0
    subsystem_result_count = 0
    cdd_cluster_count = 0

    def get_cdd_db_entry(id: str) -> DBentry:
        try:
            return CDD[int(id[4:])]
        except:
            raise Exception(f'While parsing cdd results, came across an entry that is not in cddid.tbl: {id}')
    with TabularBlastParser(result_files[0], 'BLAST', get_cdd_db_entry) as handle:
        for cdd_result in handle:
            feature = sqlite.read_feature_by_id(db_connection, cdd_result.query())
            if not feature:
                context.log(f'({genome.name}) FATAL ERROR: Found {ANNOTATOR_KEY} result for unknown feature {cdd_result.query()}, '
                                f'may need to rerun metaerg with --force')
                raise Exception(f'({genome.name}) Found {ANNOTATOR_KEY} result for unknown feature {cdd_result.query()}, '
                                f'may need to rerun metaerg with --force')
            # process the feature's cdd
            feature.cdd = cdd_result
            cdd_result_count += 1
            for new_match in functional_genes.match(cdd_result, number_of_hits_considered=5):
                if not new_match in feature.subsystems:
                    feature.subsystems.append(new_match)
                    subsystem_result_count += 1
            top_entry: DBentry = cdd_result.hits[0].hit
            feature.descr = f'{top_entry.accession}|{top_entry.gene} {top_entry.descr}'
            if len(feature.descr) > 35:
                feature.descr = feature.descr[:35] + '...'
            # write feature
            sqlite.update_feature_in_db(db_connection, feature)

    context.log(f'({genome.name}) Found {subsystem_result_count} matches to subsystems.')
    return cdd_result_count


def preload_db():
    cdd_index = Path(context.DATABASE_DIR, 'cdd', 'cddid.tbl')
    with open(cdd_index) as db_handle:
        for line in db_handle:
            words = line.strip().split("\t")
            CDD[int(words[0])] = DBentry(domain='cdd', accession=words[1], gene=words[2], descr=words[3],
                                         length=int(words[4]))
    context.log(f'Parsed {len(CDD)} entries from conserved domain database.')
    cdd_cluster_db = Path(context.DATABASE_DIR, 'cdd', 'cddid.clusters')
    if cdd_cluster_db.exists():
        with open(cdd_cluster_db) as db_handle:
            for line in db_handle:
                words = line.split('\t')
                CDD_CLUSTERS[words[0]] = float(words[1])
    context.log(f'Parsed {len(CDD_CLUSTERS)} gene-context-associations for profiles.')



@context.register_annotator
def run_and_read_cdd():
    return ({'pipeline_position': 71,
             'annotator_key': ANNOTATOR_KEY,
             'purpose': 'function prediction using RPSBlast and the conserved domain database',
             'programs': ('rpsblast',),
             'databases': (Path('cdd', 'cddid.tbl'), Path('cdd', 'Cdd.pal')),
             'result_files': ('cdd',),
             'run': _run_programs,
             'read': _read_results,
             'preload_db': preload_db})


@context.register_database_installer
def install_cdd_database():
    if 'C' not in context.DATABASE_TASKS:
       return
    cdd_dir = context.DATABASE_DIR / 'cdd'
    context.log(f'Installing the conserved domain database to {cdd_dir}...')
    if context.FORCE_INSTALLATION_OF_DB or not cdd_dir.exists():
        if cdd_dir.exists():
            shutil.rmtree(cdd_dir)
        cdd_dir.mkdir(exist_ok=True, parents=True)
        context.run_external(f'wget -P {cdd_dir} https://ftp.ncbi.nih.gov/pub/mmdb/cdd/cddid.tbl.gz')
        context.run_external(f'wget -P {cdd_dir} https://ftp.ncbi.nih.gov/pub/mmdb/cdd/little_endian/Cdd_LE.tar.gz')
        context.run_external(f'gunzip {cdd_dir / "cddid.tbl.gz"}')
        context.run_external(f'tar -C {cdd_dir} -xf {cdd_dir / "Cdd_LE.tar.gz"}')
        (cdd_dir / "Cdd_LE.tar.gz").unlink()
    else:
        context.log(f'{cdd_dir} already exists, use --force all to overwrite existing installation of cdd')
