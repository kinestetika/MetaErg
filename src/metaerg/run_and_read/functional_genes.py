import os
import gzip
import shutil
from pathlib import Path
from concurrent import futures
from zipfile import ZipFile

from metaerg.datatypes.blast import DBentry, TabularBlastParser
from metaerg.datatypes import sqlite, fasta
from metaerg import context
from metaerg.datatypes import functional_genes
from metaerg.calculations.codon_usage_bias import compute_codon_bias_estimate_doubling_time

ANNOTATOR_KEY = 'hmm'

def _run_programs(genome, contig_dict, db_connection, result_files):
    cds_aa_file = context.spawn_file('cds.faa', genome.name)
    hmm_db = context.DATABASE_DIR / 'hmm' / 'functional_genes.hmm'
    # the domtbl format includes alignment starts and ends, so we can use that information...
    context.run_external(f'hmmscan -o /dev/null -E 1e-6 --domtblout {result_files[0]} {hmm_db} {cds_aa_file}')


def _read_results(genome, contig_dict, db_connection, result_files) -> int:
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
    context.log(f'({genome.name}) {len(functional_gene_db)} functional gene profiles in database.')

    def get_db_entry(db_id) -> DBentry:
        try:
            return functional_gene_db[db_id]
        except KeyError:
            context.log(f'({genome.name}) Unkown hmm database id "{db_id}". Did the hmm database change after this prediction was done?')
            raise Exception(f'({genome.name}) Unkown hmm database id "{db_id}". Did the hmm database change after this prediction was done?')


    with TabularBlastParser(result_files[0], 'HMMSCAN_DOM_TABLE', get_db_entry) as handle:
        hit_count = 0
        for blast_result in handle:
            feature = sqlite.read_feature_by_id(db_connection, blast_result.query())
            if not feature:
                context.log(
                    f'({genome.name}) FATAL ERROR: Found {ANNOTATOR_KEY} result for unknown feature {blast_result.query()}, '
                    f'may need to rerun metaerg with --force')
                raise Exception(f'Found {ANNOTATOR_KEY} result for unknown feature {blast_result.query()}, '
                                f'may need to rerun metaerg with --force all')
            if new_matches := functional_genes.match(blast_result, number_of_hits_considered=1):
                hit_count += 1
                for new_match in new_matches:
                    if not new_match in feature.subsystems:
                        feature.subsystems.append(new_match)
            blast_result.hits = blast_result.hits[:5]
            feature.hmm = blast_result
            sqlite.update_feature_in_db(db_connection, feature)
    # with all subsystems annotated, we are ready to update the genome's subsystems:
    genome.codon_usage_bias, genome.doubling_time = compute_codon_bias_estimate_doubling_time(genome.name, db_connection)
    genome.subsystems = functional_genes.aggregate(db_connection)
    for subsystem, subsystem_genes in genome.subsystems.items():
        subsystem_completeness = functional_genes.get_subsystem_completeness(subsystem, subsystem_genes)
        genome.subsystem_summary[subsystem] = subsystem_completeness

    return hit_count


@context.register_annotator
def run_and_read_canthyd():
    return ({'pipeline_position': 101,
             'annotator_key': ANNOTATOR_KEY,
             'purpose': 'identification of functional genes with HMMs',
             'programs': ('hmmscan',),
             'databases': (Path('hmm', 'functional_genes.hmm'),),
             'result_files': ('hmm',),
             'run': _run_programs,
             'read': _read_results,
             'preload_db': functional_genes.init_functional_gene_config})


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
    if 'S' not in context.DATABASE_TASKS:
        return

    functional_genes.init_functional_gene_config()
    all_target_hmm_names = set()
    for gene_def in functional_genes.GENES:
        all_target_hmm_names.update(gene_def.cues)
        all_target_hmm_names.update(gene_def.anti_cues)

    hmm_dir = context.DATABASE_DIR / 'hmm'
    hmm_dir.mkdir(exist_ok=True, parents=True)
    user_hmm_dir = hmm_dir / 'user_hmm'
    (hmm_dir / 'user_config').mkdir(exist_ok=True, parents=True)  # used by functional_genes.py
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
            context.run_external(f'hmmconvert {hmm_file}', stdout=writer)

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

    # fetch nitrite reductase data from Pold et al. (2024) ISME Reports...
    nir_hmm_file = hmm_dir / 'nir.hmm'
    hmm_nir_install_dir = hmm_dir / 'nir'
    hmm_nir_install_dir.mkdir(exist_ok=True, parents=True)
    os.chdir(hmm_nir_install_dir)
    # see https://api.figshare.com/v2/articles/23913078
    os.system('wget -q https://ndownloader.figshare.com/files/42604285')
    nir_archive_file = hmm_nir_install_dir / '42604285'
    with ZipFile(nir_archive_file, 'r') as zipped_file:
        zipped_file.extractall(path=hmm_nir_install_dir)
    #nir_hmm_archive_file = hmm_nir_install_dir / 'clade_specific_HMMs'
    #with ZipFile(nir_hmm_archive_file, 'r') as zipped_file:
    #    zipped_file.extractall(path=hmm_nir_install_dir)
    nirk_count = 0
    nirs_count = 0
    with open(nir_hmm_file, 'w') as nir_hmm_writer:
        for f in sorted((hmm_nir_install_dir / 'clade_specific_HMMs' / 'NirK').glob('*.hmm')):
            with open(f, 'r') as reader:
                for l in reader:
                    if l.startswith('NAME'):
                        l = 'NAME  NirK_' + l[6:]
                    nir_hmm_writer.write(l)
                    nirk_count += 1

        for f in sorted((hmm_nir_install_dir / 'clade_specific_HMMs' / 'NirS').glob('*.hmm')):
            with open(f, 'r') as reader:
                for l in reader:
                    if l.startswith('NAME  Clade') or l.startswith('NAME  Arch'):
                        l = 'NAME  NirS_' + l[6:]
                    nir_hmm_writer.write(l)
                    nirs_count += 1
    context.log(f'Downloaded and added {nirs_count} and {nirk_count} profiles for NirS and Nir respectively.')
    shutil.rmtree(hmm_nir_install_dir)

    # unzip hmms from metaerg python module
    hmm_data_dir = Path(__file__).parent / 'data'
    context.log(f'Now unpacking hmm data from {hmm_data_dir}')
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
            if context.FORCE_INSTALLATION_OF_DB or not destination_file.exists() or not destination_file.stat().st_size:
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


def collect_nir_hmms():
    context.log(f'Installing Nir HMMs from Pold et al. (2024) ISME Reports...')
    hmm_nir_install_dir = context.DATABASE_DIR / 'hmm' / 'nir'
    hmm_nir_install_dir.mkdir(exist_ok=True, parents=True)
    os.chdir(hmm_nir_install_dir)
    os.system('wget -q https://figshare.com/ndownloader/articles/23913078/versions/1')
    nir_archive_file = hmm_nir_install_dir / '1'
    with ZipFile(nir_archive_file, 'r') as zipped_file:
        zipped_file.extractall(path=hmm_nir_install_dir)
    nir_hmm_archive_file = hmm_nir_install_dir / 'clade_specific_HMMs.zip'
    with ZipFile(nir_hmm_archive_file, 'r') as zipped_file:
        zipped_file.extractall(path=hmm_nir_install_dir)
    nir_hmm_archive_dir = hmm_nir_install_dir / 'clade_specific_HMMs'
    nirk_hmm_archive_dir = nir_hmm_archive_dir / 'NirK'
    nirk_fasta_dir = nirk_hmm_archive_dir / 'fasta'
    nirk_hmm_dir = nirk_hmm_archive_dir / 'hmm'
    nirk_fasta_dir.mkdir(exist_ok=True, parents=True)
    nirk_hmm_dir.mkdir(exist_ok=True, parents=True)
    for f in nirk_hmm_archive_dir.glob('*.hmm'):
        with open(f, 'r') as reader, open(nirk_hmm_dir / f.name, 'w') as writer:
            for l in reader:
                if l.startswith('NAME'):
                    l = 'NAME  NirK_' + l[6:]
                writer.write(l)
    for f in nirk_hmm_archive_dir.glob('*.fasta'):
        with fasta.FastaParser(f, cleanup_seq=True) as reader, open(nirk_fasta_dir / f.name, 'w') as writer:
            for fe in reader:
                fe['seq'] = fe['seq'].replace('X', '')
                fasta.write_fasta(writer, fe)
    nirs_hmm_archive_dir = nir_hmm_archive_dir / 'NirS'
    nirs_fasta_dir = nirs_hmm_archive_dir / 'fasta'
    nirs_hmm_dir = nirs_hmm_archive_dir / 'hmm'
    nirs_fasta_dir.mkdir(exist_ok=True, parents=True)
    nirs_hmm_dir.mkdir(exist_ok=True, parents=True)
    for f in nirs_hmm_archive_dir.glob('*.hmm'):
        with open(f, 'r') as reader, open(nirs_hmm_dir / f.name, 'w') as writer:
            for l in reader:
                if l.startswith('NAME  Clade') or l.startswith('NAME  Arch'):
                    l = 'NAME  NirS_' + l[6:]
                writer.write(l)
    for f in nirs_hmm_archive_dir.glob('*.fasta'):
        with fasta.FastaParser(f, cleanup_seq=True) as reader, open(nirs_fasta_dir / f.name, 'w') as writer:
            for fe in reader:
                fe['seq'] = fe['seq'].replace('X', '')
                fasta.write_fasta(writer, fe)

def main():
    context.DATABASE_DIR = Path('/bio/data/databases/metaerg')
    collect_nir_hmms()

if __name__ == "__main__":
    main()

