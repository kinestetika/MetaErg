import os
import gzip
import shutil
from math import log10
from pathlib import Path
from concurrent import futures

from metaerg.datatypes.blast import DBentry, TabularBlastParser
from metaerg.datatypes import sqlite
from metaerg import context
from metaerg import functional_gene_configuration


def _run_programs(genome_name, contig_dict, db_connection, result_files):
    cds_aa_file = context.spawn_file('cds.faa', genome_name)
    hmm_db = context.DATABASE_DIR / 'hmm' / 'functional_genes.hmm'
    # the domtbl format includes alignment starts and ends, so we can use that information...
    context.run_external(f'hmmscan -o /dev/null -E 1e-6 --domtblout {result_files[0]} {hmm_db} {cds_aa_file}')


def _read_results(genome_name, contig_dict, db_connection, result_files) -> int:
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
            feature = sqlite.read_feature_by_id(db_connection, blast_result.query())
            if not feature:
                raise Exception(f'Found results for unknown feature {blast_result.query()}, '
                                f'may need to rerun metaerg with --force')
            blast_result.hits = blast_result.hits[:10]
            feature.hmm = blast_result
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
                    feature.subsystems = functional_gene_configuration.cleanup_subsystem_str(
                        feature.subsystems.strip() + ' ' + new_subsystem)
                break
            sqlite.update_feature_in_db(db_connection, feature)
    return hit_count


@context.register_annotator
def run_and_read_canthyd():
    return ({'pipeline_position': 101,
             'annotator_key': 'hmm',
             'purpose': 'identification of functional genes with HMMs',
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


def add_hmm_file_to_db(hmm_file: Path, db_handle, names_done: set, targets: set):
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
                    continue
                names_done.add(name)
                if name in targets:
                    written += 1
                    for l in hmm_lines:
                        if l.startswith('ACC'):
                            pass
                        else:
                            db_handle.write(l)
                elif hmm_file.name == 'CANT-HYD.hmm':
                    print(''.join(hmm_lines[:20]))
            elif line.startswith('NAME'):
                name = line[4:].strip()
    print(f'{hmm_file.name}: read {read}, written {written}, duplicates {len(duplicates)}.')


@context.register_database_installer
def install_functional_gene_databases():
    if 'S' not in context.TASKS:
        return

    functional_gene_configuration.init_functional_gene_config()
    all_target_hmm_names = set()
    for t in functional_gene_configuration.SUBSYSTEM_DATA['profiles']:
        for t1 in t.split('|'):
            all_target_hmm_names.add(t1)

    hmm_dir = context.DATABASE_DIR / 'hmm'
    hmm_dir.mkdir(exist_ok=True, parents=True)
    user_hmm_dir = hmm_dir / 'user_hmm'
    (hmm_dir / 'user_config').mkdir(exist_ok=True, parents=True)  # used by functional_gene_configuration.py
    user_hmm_dir.mkdir(exist_ok=True, parents=True)  # used below
    context.log(f'Installing functional gene hmm database at {hmm_dir}...')

    current_working_dir = os.getcwd()
    # fetch bd type oxygen reductases
    bdor_dir = hmm_dir / 'bdor'
    if bdor_dir.exists():
        shutil.rmtree(bdor_dir)
    bdor_dir.mkdir()
    os.chdir(bdor_dir)
    os.system('git init')
    os.system('git remote add bdor https://github.com/ranjani-m/cytbd-superfamily.git')
    os.system('git fetch bdor')
    os.system('git checkout bdor/main -- Cytbd_HMM_2020_paper_version')
    with open(hmm_dir / 'bdor.hmm', 'w') as writer:
        for hmm_file in sorted((bdor_dir / 'Cytbd_HMM_2020_paper_version').glob('*.hmm')):
            os.system(f'sed -i "s|^NAME.*|NAME  bdor_{hmm_file.stem}|" {hmm_file}')
            context.run_external(f'/bio/bin/hmmer3/bin/hmmconvert {hmm_file}', stdout=writer)
    # fetch heme copper oxidase hmms
    hco_dir = hmm_dir / 'hco'
    if hco_dir.exists():
        shutil.rmtree(hco_dir)
    hco_dir.mkdir()
    os.chdir(hco_dir)
    os.system('git init')
    os.system('git remote add hco https://github.com/ranjani-m/HCO.git')
    os.system('git fetch hco')
    os.system('git checkout hco/main -- HMM')
    with open(hmm_dir / 'hco.hmm', 'w') as writer:
        for hmm_file in sorted((hco_dir / 'HMM').glob('*.hmm')):
            os.system(f'sed -i "s|^NAME.*|NAME  hco_{hmm_file.stem}|" {hmm_file}')
            context.run_external(f'hmmconvert {hmm_file}', stdout=writer)

    # unzip hmms from data
    hmm_data_dir = Path(__file__).parent / 'data'
    print (hmm_data_dir)
    for hmm_file in hmm_data_dir.glob('*.hmm.gz'):
        with gzip.open(hmm_file, 'rb') as f_in:
            dest_hmm_file = hmm_dir / hmm_file.stem
            print ('unzipping to ', dest_hmm_file)
            with open(dest_hmm_file, 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)

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
            add_hmm_file_to_db(hmm_file, output, names_done, all_target_hmm_names)
        for f in user_hmm_dir.glob('*.hmm'):
            hmm_file = f.absolute()
            add_hmm_file_to_db(hmm_file, output, names_done, all_target_hmm_names)
    context.run_external(f'hmmpress -f {hmm_dir / "functional_genes.hmm"}')
    os.chdir(current_working_dir)
