from pathlib import Path

from metaerg.run_and_read.data_model import MetaergSeqFeature, TabularBlastParser, DBentry, MetaergGenome
from metaerg.run_and_read.context import register, spawn_file, run_external, DATABASE_DIR, log


def _run_programs(genome:MetaergGenome, result_files):
    cds_aa_file = spawn_file('cds.faa', genome.id)
    canthyd_db = Path(DATABASE_DIR, 'CANT-HYD.hmm')
    run_external(f'hmmscan --cut_nc --tblout {result_files[0]} {canthyd_db} {cds_aa_file}')


def _read_results(genome:MetaergGenome, result_files) -> int:
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
    canthyd_db = Path(DATABASE_DIR, 'CANT-HYD.hmm')
    with open(canthyd_db) as handle:
        for line in handle:
            if line.startswith('NAME'):
                current_name = line.split()[1]
            elif line.startswith('TC'):
                canthyd_trusted_cutoffs[current_name] = int(line.split()[1])
    log(f'Parsed {len(canthyd_trusted_cutoffs)} entries from CantHyd database.')

    def get_db_entry(db_id) -> DBentry:
        return DBentry(db_id, '', canthyd_descr.get(db_id, ''), '', 0, canthyd_trusted_cutoffs[db_id])

    with TabularBlastParser(result_files[0], 'HMMSCAN', get_db_entry) as handle:
        canthyd_hit_count = 0
        for blast_result in handle:
            feature: MetaergSeqFeature = genome.get_feature(blast_result.query)
            for h in blast_result.hits:
                if descr := h.hit.descr:
                    canthyd_hit_count += 1
                    confidence = 'high' if h.score > h.hit.pos else 'low'  # cutoff is stored in 'pos'
                    feature.product = f'{descr}, {h.hit.id} (CantHyd DB, {confidence} confidence)'
                    feature.subsystem.add('[hydrocarbon degradation]')
                    genome.subsystems.subsystems['[hydrocarbon degradation]'].add_hit(feature.id)
                else:
                    log(f'Warning, missing description for cant-hyd hmm {h.hit}...')
        return canthyd_hit_count


@register
def run_and_read_canthyd():
    return ({'pipeline_position': 101,
             'purpose': 'prediction of hydrocarbon degradation genes with canthyd',
             'programs': ('hmmscan'),
             'databses': ('CANT-HYD.hmm'),
             'result_files': ('canthyd'),
             'run': _run_programs,
             'read': _read_results})
