import pandas as pd
from pathlib import Path

from metaerg.datatypes.blast import BlastResult, DBentry, TabularBlastParser
from metaerg import context
from metaerg import subsystems


def _run_programs(genome_name, contig_dict, feature_data: pd.DataFrame, result_files):
    cds_aa_file = context.spawn_file('cds.faa', genome_name)
    hmm_db = context.DATABASE_DIR / 'hmm' / 'functional_genes.hmm'
    context.run_external(f'hmmscan -E e-6 --tblout {result_files[0]} {hmm_db} {cds_aa_file}')


def _read_results(genome_name, contig_dict, feature_data: pd.DataFrame, result_files) -> tuple:
    db_entry = None
    hmm_db = context.DATABASE_DIR / 'hmm' / 'functional_genes.hmm'
    functional_gene_db = {}
    with open(hmm_db) as handle:
        for line in handle:
            if line.startswith('NAME'):
                name = line.split()[1]
                db_entry = DBentry(domain='functional_genes', accession=name, gene=name)
                functional_gene_db[name] = db_entry
            elif line.startswith('ACC'):
                db_entry.accession = line[4:].strip()
            elif line.startswith('NC'):
                db_entry.min_score = int(line[4:].strip())
            elif line.startswith('TC'):
                db_entry.min_t_score = int(line[4:].strip())
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
            for h in blast_result.hits:
                db_entry = h.hit
                confidence = 'similar to '
                if db_entry.min_score:
                    if h.score < db_entry.min_score:
                        continue
                    if h.score > db_entry.min_t_score:
                        confidence = ''
                else:
                    if h.evalue > 1e-80:
                        continue
                    if h.evalue < 1e-100:
                        confidence = ''
                if descr := h.hit.descr:
                    hit_count += 1
                    feature_data.at[blast_result.query, 'descr'] = \
                        f'{confidence}{descr}'
                    feature_data.at[blast_result.query(), 'subsystems'] = subsystems.match_hit(h)
                else:
                    context.log(f'Warning, missing description for hmm {h.hit}...')
                break
        return feature_data, hit_count


@context.register_annotator
def run_and_read_canthyd():
    return ({'pipeline_position': 101,
             'purpose': 'prediction of hydrocarbon degradation genes with canthyd',
             'programs': ('hmmscan',),
             'databases': (Path('canthyd', 'CANT-HYD.hmm'),),
             'result_files': ('canthyd',),
             'run': _run_programs,
             'read': _read_results})


@context.register_database_installer
def install_canthyd_database():
    if 'S' not in context.CREATE_DB_TASKS:
        return
    canthyd_dir = context.DATABASE_DIR / 'canthyd'
    if context.FORCE or not canthyd_dir.exists():
        context.log(f'Installing the conserved domain database to {canthyd_dir}...')
        canthyd_dir.mkdir(exist_ok=True, parents=True)
        context.run_external(f'wget -P {canthyd_dir} https://github.com/dgittins/CANT-HYD-HydrocarbonBiodegradation/raw/'
                             f'main/HMMs/concatenated%20HMMs/CANT-HYD.hmm')
        context.run_external(f'hmmpress -f {canthyd_dir / "CANT-HYD.hmm"}')
    else:
        context.log(f'Keeping existing cangthyd database in {canthyd_dir}, use --force to overwrite.')

https://zenodo.org/record/6365663/files/carbon.cycle.sub.hmm