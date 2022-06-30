import pandas as pd
from pathlib import Path

from metaerg.datatypes.blast import BlastResult, DBentry, TabularBlastParser
from metaerg import context


def _run_programs(genome_name, contig_dict, feature_data: pd.DataFrame, result_files):
    cds_aa_file = context.spawn_file('cds.faa', genome_name)
    canthyd_db = Path(context.DATABASE_DIR, 'canthyd', 'CANT-HYD.hmm')
    context.run_external(f'hmmscan --cut_nc --tblout {result_files[0]} {canthyd_db} {cds_aa_file}')


def _read_results(genome_name, contig_dict, feature_data: pd.DataFrame, result_files) -> tuple:
    canthyd_trusted_cutoffs = {}
    canthyd_descr = {'AlkB': 'alkane hydrolase',
                     'AlmA_GroupI': 'flavin-binding alkane monooxygenase',
                     'AlmA_GroupIII': 'flavin-binding alkane monooxygenase',
                     'CYP153': 'alkane oxidizing cytochrome P450',
                     'LadA_alpha': 'long-chain alkane hydrolase',
                     'LadA_beta': 'long-chain alkane hydrolase',
                     'LadB': 'long-chain alkane hydrolase',
                     'pBmoA': 'membrane-bound alkane monooxygenase subunit A',
                     'pBmoB': 'membrane-bound alkane monooxygenase subunit B',
                     'pBmoC': 'membrane-bound alkane monooxygenase subunit C',
                     'PrmA': 'propane 2-monooxygenase large subunit',
                     'PrmC': 'propane 2-monooxygenase small subunit',
                     'sBmoX': 'soluble alkane monooxygenase subunit A',
                     'sBmoY': 'soluble alkane monooxygenase subunit B',
                     'DmpO': 'phenol/toluene 2-monooxygenase (NADH dependent)',
                     'DszC': 'dibenzothiophene monooxygenase',
                     'MAH_alpha': 'benzene/toluene/naphtalene dioxygenase subunit alpha',
                     'MAH_beta': 'benzene/toluene/naphtalene dioxygenase subunit beta',
                     'NdoB': 'benzene/toluene/naphtalene dioxygenase subunit alpha',
                     'non_NdoB_type': 'similar to benzene/toluene/naphtalene dioxygenase subunit alpha',
                     'NdoC': 'benzene/toluene/naphtalene dioxygenase subunit beta',
                     'TmoA_BmoA': 'toluene monooxygenase subunit A',
                     'TmoB_BmoB': 'toluene monooxygenase subunit B',
                     'TmoE': 'toluene monooxygenase system protein E',
                     'TomA1': 'phenol/toluene monooxygenase/hydroxylase (NADH dependent)',
                     'TomA3': 'phenol/toluene monooxygenase/hydroxylase (NADH dependent)',
                     'TomA4': 'phenol/toluene monooxygenase/hydroxylase (NADH dependent)',
                     'ahyA': 'molybdopterin-family alkane C2 methylene hydroxylase',
                     'AssA': 'alkylsuccinate synthase',
                     'AbcA_1': 'benzene carboxylase',
                     'BssA': 'benzylsuccinate synthase',
                     'CmdA': 'molybdopterin-family ethylbenzene dehydrogenase subunit alpha',
                     'EbdA': 'molybdopterin-family ethylbenzene dehydrogenase subunit alpha',
                     'K27540': 'naphtalene carboxylase',
                     'NmsA': 'naphtylmethyl succinate synthase'}
    current_name = None
    canthyd_db = Path(context.DATABASE_DIR, 'canthyd', 'CANT-HYD.hmm')
    with open(canthyd_db) as handle:
        for line in handle:
            if line.startswith('NAME'):
                current_name = line.split()[1]
            elif line.startswith('TC'):
                canthyd_trusted_cutoffs[current_name] = int(line.split()[1])
    context.log(f'({genome_name}) Parsed {len(canthyd_trusted_cutoffs)} entries from CantHyd database.')

    def get_db_entry(db_id) -> DBentry:
        return DBentry(domain='canthyd', gene=db_id, descr=canthyd_descr.get(db_id, ''),
                       pos=canthyd_trusted_cutoffs[db_id])

    with TabularBlastParser(result_files[0], 'HMMSCAN', get_db_entry) as handle:
        canthyd_hit_count = 0
        for blast_result in handle:
            for h in blast_result.hits:
                if descr := h.hit.descr:
                    canthyd_hit_count += 1
                    confidence = 'high' if h.score > h.hit.pos else 'low'  # cutoff is stored in 'pos'
                    feature_data.at[blast_result.query, 'descr'] = \
                        f'{descr}, {h.hit.id} (CantHyd DB, {confidence} confidence)'
                    if len(feature_data.at[blast_result.query(), 'subsystems']):
                        if '[hydrocarbon degradation]' not in feature_data.at[blast_result.query(), 'subsystems']:
                            feature_data.at[blast_result.query(), 'subsystems'] += ' [hydrocarbon degradation]'
                    else:
                        feature_data.at[blast_result.query(), 'subsystems'] = '[hydrocarbon degradation]'
                else:
                    context.log(f'Warning, missing description for cant-hyd hmm {h.hit}...')
        return feature_data, canthyd_hit_count


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
    canthyd_dir = Path(context.DATABASE_DIR, 'canthyd')
    if context.FORCE or not canthyd_dir.exists():
        context.log(f'Installing the conserved domain database to {canthyd_dir}...')
        canthyd_dir.mkdir(exist_ok=True, parents=True)
        context.run_external(f'wget -P {canthyd_dir} https://github.com/dgittins/CANT-HYD-HydrocarbonBiodegradation/raw/'
                             f'main/HMMs/concatenated%20HMMs/CANT-HYD.hmm')
        context.run_external(f'hmmpress -f {Path(canthyd_dir, "CANT-HYD.hmm")}')
    else:
        context.log(f'Keeping existing cangthyd database in {canthyd_dir}, use --force to overwrite.')
