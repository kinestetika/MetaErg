from pathlib import Path
from metaerg.run_and_read.data_model import MetaergSeqFeature
from metaerg.run_and_read import abc
from metaerg import utils


class CantHyd(abc.AbstractBaseClass):
    def __init__(self, genome, exec_env: abc.ExecutionEnvironment):
        super().__init__(genome, exec_env)
        self.canthyd_file = self.spawn_file('canthyd')
        self.db_canthyd = Path(self.exec.database_dir, "canthyd", "CANT-HYD.hmm")
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

    def __repr__(self):
        return f'CantHyd({self.genome}, {self.exec})'

    def __purpose__(self) -> str:
        """Should return the purpose of the tool"""
        return 'prediction of hydrocarbon degradation genes with canthyd'

    def __programs__(self) -> tuple:
        """Should return a tuple with the programs needed"""
        return 'hmmscan',

    def __databases__(self) -> tuple:
        """Should return a tuple with database files needed"""
        return self.db_canthyd,

    def __result_files__(self) -> tuple:
        """Should return a tuple with the result files (Path objects) created by the programs"""
        return self.canthyd_file,

    def __run_programs__(self):
        """Should execute the helper programs to complete the analysis"""
        cds_aa_file = self.spawn_file('cds.faa')
        utils.run_external(f'hmmscan --cut_nc --tblout {self.canthyd_file} {self.db_canthyd} {cds_aa_file}')

    def __read_results__(self) -> int:
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
        with utils.TabularBlastParser(self.canthyd_file, 'HMMSCAN') as handle:
            canthyd_hit_count = 0
            for blast_result in handle:
                feature: MetaergSeqFeature = self.genome.get_feature(blast_result.query)
                for h in blast_result.hits:
                    descr = self.canthyd_descr.get(h.hit, None)
                    if descr:
                        canthyd_hit_count += 1
                        confidence = 'high' if h.score > self.canthyd_trusted_cutoffs[h.hit] else 'low'
                        feature.description = f'{descr}, {h.hit} (CantHyd DB, {confidence} confidence)'
                        feature.subsystem.append('[hydrocarbon degradation]')
                    else:
                        utils.log(f'Warning, missing description for cant-hyd hmm {h.hit}...')
            return canthyd_hit_count
