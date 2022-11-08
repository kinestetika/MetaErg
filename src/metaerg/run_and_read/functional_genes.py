import os

import pandas as pd
from math import log10
from pathlib import Path
from concurrent import futures

from metaerg.datatypes.blast import DBentry, TabularBlastParser
from metaerg import context
from metaerg import functional_gene_configuration


def _run_programs(genome_name, contig_dict, feature_data: pd.DataFrame, result_files):
    cds_aa_file = context.spawn_file('cds.faa', genome_name)
    hmm_db = context.DATABASE_DIR / 'hmm' / 'functional_genes.hmm'
    # the domtbl format includes alignment starts and ends, so we can use that information...
    context.run_external(f'hmmscan -o /dev/null -E 1e-6 --domtblout {result_files[0]} {hmm_db} {cds_aa_file}')


def _read_results(genome_name, contig_dict, feature_data: pd.DataFrame, result_files) -> tuple:
    db_entry = None
    hmm_db = context.DATABASE_DIR / 'hmm' / 'functional_genes.hmm'
    functional_gene_db = {}
    with open(hmm_db) as handle:
        for line in handle:
            if line.startswith('NAME'):
                name = line.split()[1]
                db_entry = DBentry(domain='functional_genes', descr='', accession=name, gene=name)
                functional_gene_db[name] = db_entry
            elif line.startswith('ACC'):
                db_entry.accession = line[4:].strip()
            elif line.startswith('NC'):
                db_entry.min_score = float(line[4:].strip().split()[0])
            elif line.startswith('TC'):
                db_entry.min_t_score = float(line[4:].strip().split()[0])
            elif line.startswith('DESC'):
                db_entry.descr = line[4:].strip()
            elif line.startswith('LENG'):
                db_entry.length = int(line[4:].strip())
    context.log(f'({genome_name}) {len(functional_gene_db)} functional gene profiles in database.')

    def get_db_entry(db_id) -> DBentry:
        return functional_gene_db[db_id]

    with TabularBlastParser(result_files[0], 'HMMSCAN_DOM_TABLE', get_db_entry) as handle:
        hit_count = 0
        for blast_result in handle:
            if blast_result.query() not in feature_data.index:
                raise Exception(f'Found results for unknown feature {blast_result.query()}, '
                                f'may need to rerun metaerg with --force')
            for h in blast_result.hits:
                db_entry = h.hit
                confidence = 1.0
                if h.aligned_length / db_entry.length < 0.7:
                    continue
                if db_entry.min_score:
                    if h.score < db_entry.min_score:
                        continue
                    confidence = h.score / db_entry.min_t_score
                else:
                    if h.evalue > 0:
                        confidence = min(1.0, - log10(h.evalue) / 100)
                if new_subsystem := functional_gene_configuration.match_hit(h, confidence):
                    hit_count += 1
                    feature_data.at[blast_result.query(), 'subsystems'] = functional_gene_configuration.cleanup_subsystem_str(
                        feature_data.at[blast_result.query(), 'subsystems'].strip() + ' ' + new_subsystem)
                break
    return feature_data, hit_count


@context.register_annotator
def run_and_read_canthyd():
    return ({'pipeline_position': 101,
             'purpose': 'identification of functional genes with canthyd and metascan databases',
             'programs': ('hmmscan',),
             'databases': (Path('hmm', 'functional_genes.hmm'),),
             'result_files': ('hmm',),
             'run': _run_programs,
             'read': _read_results})

FUNCTIONAL_GENE_URLS = ('https://zenodo.org/record/6365663/files/carbon.cycle.sub.hmm',
                        'https://zenodo.org/record/6365663/files/c1.cycle.sub.hmm',
                        'https://zenodo.org/record/6365663/files/methane.cycle.sub.hmm',
                        'https://zenodo.org/record/6365663/files/misc.cycle.sub.hmm',
                        'https://zenodo.org/record/6365663/files/nitro.cycle.sub.hmm',
                        'https://zenodo.org/record/6365663/files/oxygen.cycle.sub.hmm',
                        'https://zenodo.org/record/6365663/files/sulphur.cycle.sub.hmm',
                        'https://zenodo.org/record/6365663/files/hydrogenase.hmm',
                        'https://github.com/dgittins/CANT-HYD-HydrocarbonBiodegradation/raw/main/HMMs/concatenated%20HMMs/CANT-HYD.hmm')


def add_hmm_file_to_db(hmm_file: Path, db_handle, names_done: set):
    with open(hmm_file) as hmm_reader:
        duplicates = []
        written = 0
        read = 0
        hmm_lines = []
        for line in hmm_reader:
            if line.startswith('HMMER3/f'):
                hmm_lines.clear()
                read += 1
                name = ''
            hmm_lines.append(line)
            if line.startswith('//'):
                if name in names_done:  # avoid duplicates
                    duplicates.append(name)
                else:
                    written += 1
                    for l in hmm_lines:
                        if l.startswith('ACC'):
                            pass
                        else:
                            db_handle.write(l)
                    names_done.add(name)
            elif line.startswith('NAME'):
                name = line[4:].strip()
    print(hmm_file.name, 'read, written, duplicates', read, written, len(duplicates))


@context.register_database_installer
def install_functional_gene_databases():
    if 'S' not in context.TASKS:
        return
    hmm_dir = context.DATABASE_DIR / 'hmm'
    hmm_dir.mkdir(exist_ok=True, parents=True)
    (hmm_dir / 'user_config').mkdir(exist_ok=True, parents=True)
    (hmm_dir / 'user_hmm').mkdir(exist_ok=True, parents=True)
    context.log(f'Installing functional gene hmm database at {hmm_dir}...')

    with futures.ThreadPoolExecutor() as executor:
        outcomes = []
        for url in FUNCTIONAL_GENE_URLS:
            destination_file = hmm_dir / Path(url).name
            if context.FORCE or not destination_file.exists() or not destination_file.stat().st_size:
                outcomes.append(executor.submit(context.download, url, destination_file))
        for future in futures.as_completed(outcomes):
            future.result()

    context.log('Checking data sanity...')
    all_hmms_file = hmm_dir / 'functional_genes.hmm'
    with open(all_hmms_file, 'w') as output:
        names_done = set()
        for f in hmm_dir.glob('*.hmm'):
            hmm_file = f.absolute()
            if f == all_hmms_file:
                continue
            add_hmm_file_to_db(hmm_file, output, names_done)
        for f in (hmm_dir / 'user_hmm').glob('*.hmm'):
            hmm_file = f.absolute()
            add_hmm_file_to_db(hmm_file, output, names_done)
    context.run_external(f'hmmpress -f {hmm_dir / "functional_genes.hmm"}')
