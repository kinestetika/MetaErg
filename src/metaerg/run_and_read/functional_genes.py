import pandas as pd
from pathlib import Path
from concurrent import futures

from metaerg.datatypes.blast import DBentry, TabularBlastParser
from metaerg import context
from metaerg import subsystems


def _run_programs(genome_name, contig_dict, feature_data: pd.DataFrame, result_files):
    cds_aa_file = context.spawn_file('cds.faa', genome_name)
    hmm_db = context.DATABASE_DIR / 'hmm' / 'functional_genes.hmm'
    context.run_external(f'hmmscan -E 1e-6 --tblout {result_files[0]} {hmm_db} {cds_aa_file}')


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
                db_entry.min_score = int(line[4:].strip().split()[0])
            elif line.startswith('TC'):
                db_entry.min_t_score = int(line[4:].strip().split()[0])
            elif line.startswith('DESC'):
                db_entry.descr = line[4:].strip()
            elif line.startswith('LENG'):
                db_entry.length = int(line[4:].strip())
    context.log(f'({genome_name}) {len(functional_gene_db)} functional gene profiles in database.')

    def get_db_entry(db_id) -> DBentry:
        return functional_gene_db[db_id]

    with TabularBlastParser(result_files[0], 'HMMSCAN', get_db_entry) as handle:
        hit_count = 0
        for blast_result in handle:
            if blast_result.query() not in feature_data.index:
                raise Exception(f'Found results for unknown feature {blast_result.query()}, '
                                f'may need to rerun metaerg with --force')
            for h in blast_result.hits:
                db_entry = h.hit
                confidence = 'similar to '
                if db_entry.min_score:
                    if h.score < db_entry.min_score:
                        continue
                    if h.score > db_entry.min_t_score:
                        confidence = ''
                else:
                    if h.evalue > 1e-25:
                        continue
                    if h.evalue < 1e-100:
                        confidence = ''
                if descr := h.hit.descr:
                    hit_count += 1
                    # feature_data.at[blast_result.query(), 'descr'] = f'{confidence}{descr}'
                    feature_data.at[blast_result.query(), 'subsystems'] = subsystems.match_hit(h)
                else:
                    context.log(f'Warning, missing description for hmm {h.hit}...')
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

@context.register_database_installer
def install_functional_gene_databases():
    if 'S' not in context.CREATE_DB_TASKS:
        return
    hmm_dir = context.DATABASE_DIR / 'hmm'
    if context.FORCE or not hmm_dir.exists():
        context.log(f'Installing functional gene databases at {hmm_dir}...')
        hmm_dir.mkdir(exist_ok=True, parents=True)
        with futures.ThreadPoolExecutor() as executor:
            outcomes = []
            for url in FUNCTIONAL_GENE_URLS:
                destination_file = hmm_dir / Path(url).name
                outcomes.append(executor.submit(context.download, url, destination_file))
            for future in futures.as_completed(outcomes):
                future.result()
        context.log('Checking data sanity...')

        all_hmms_file = hmm_dir / 'functional_genes.hmm'
        with open(all_hmms_file, 'w') as output:
            names_done = set()
            for url in FUNCTIONAL_GENE_URLS:
                hmm_file = hmm_dir / Path(url).name
                print(hmm_file)
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
                                        output.write(l)
                                names_done.add(name)
                        elif line.startswith('NAME'):
                            name = line[4:].strip()
                print(hmm_file.name, 'read, written, duplicates', read, written, len(duplicates))
        context.run_external(f'hmmpress -f {hmm_dir / "functional_genes.hmm"}')
    else:
        context.log(f'Keeping existing functional gene database in {hmm_dir}, use --force to overwrite.')

