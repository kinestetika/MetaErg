from pathlib import Path
from metaerg.run_and_read.data_model import MetaergSeqFeature, TabularBlastParser, DBentry
from metaerg.run_and_read.abc import Annotator, ExecutionEnvironment, register
from metaerg import utils


@register
class CantHyd(Annotator):
    def __init__(self, genome, exec_env: ExecutionEnvironment):
        super().__init__(genome, exec_env)
        self.canthyd_file = self.spawn_file('canthyd')
        self.db_canthyd = Path(self.exec.database_dir, "canthyd", "CANT-HYD.hmm")
        self._pipeline_position = 101
        self._purpose = 'prediction of hydrocarbon degradation genes with canthyd'
        self._programs = ('hmmscan',)
        self._databases = (self.db_canthyd,)
        self._result_files = (self.canthyd_file,)
        self.canthyd_trusted_cutoffs = {}
        self.canthyd_descr = {'AlkB': 'alkane hydrolase',
                 'AlmA_GroupI':	'flavin-binding alkane monooxygenase',
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

    def _run_programs(self):
        """Executes the helper programs to complete the analysis"""
        cds_aa_file = self.spawn_file('cds.faa')
        utils.run_external(f'hmmscan --cut_nc --tblout {self.canthyd_file} {self.db_canthyd} {cds_aa_file}')

    def get_db_entry(self, db_id) -> DBentry:
        return DBentry(db_id, '', self.canthyd_descr.get(db_id, ''), '', 0, self.canthyd_trusted_cutoffs[db_id])

    def _read_results(self) -> int:
        """Should parse the result files and return the # of positives."""
        if self.db_canthyd.exists():
            current_name = None
            with open(self.db_canthyd) as handle:
                for line in handle:
                    if line.startswith('NAME'):
                        current_name = line.split()[1]
                    elif line.startswith('TC'):
                        self.canthyd_trusted_cutoffs[current_name] = int(line.split()[1])
            utils.log(f'Parsed {len(self.canthyd_trusted_cutoffs)} entries from CantHyd database.')
        with TabularBlastParser(self.canthyd_file, 'HMMSCAN', self.get_db_entry) as handle:
            canthyd_hit_count = 0
            for blast_result in handle:
                feature: MetaergSeqFeature = self.genome.get_feature(blast_result.query)
                for h in blast_result.hits:
                    if descr := h.hit.descr:
                        canthyd_hit_count += 1
                        confidence = 'high' if h.score > h.hit.pos else 'low'  # cutoff is stored in 'pos'
                        feature.product = f'{descr}, {h.hit.id} (CantHyd DB, {confidence} confidence)'
                        feature.subsystem.add('[hydrocarbon degradation]')
                        self.genome.subsystems.subsystems['[hydrocarbon degradation]'].add_hit(feature.id)
                    else:
                        utils.log(f'Warning, missing description for cant-hyd hmm {h.hit}...')
            return canthyd_hit_count
