from pathlib import Path
import shutil
from collections import namedtuple
from concurrent.futures import ProcessPoolExecutor
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
        return 'hmmsearch',

    def __databases__(self) -> tuple:
        """Should return a tuple with database files needed"""
        return self.db_canthyd,

    def __result_files__(self) -> tuple:
        """Should return a tuple with the result files (Path objects) created by the programs"""
        return self.canthyd_file,

    def __run_programs__(self):
        """Should execute the helper programs to complete the analysis"""
        cds_aa_file = self.spawn_file('cds.faa')
        utils.run_external(f'hmmsearch --cut_nc --tblout {self.canthyd_file} {self.db_canthyd} {cds_aa_file}')

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

        with utils.TabularHMMSearchParser(self.canthyd_file) as handle:
            for blast_result in handle:
                feature: MetaergSeqFeature = self.genome.get_feature(blast_result.query)


                for h in hits.values():
                    if h["hmm_id"] not in databases.CANTHYD_DESCR.keys():
                        utils.log(f'Warning, missing description for cant-hyd hmm {h["hmm_id"]}...')
                        continue
                    coord = utils.decipher_metaerg_id(h['hit_id'])
                    feature = contig_dict[coord["contig_id"]].features[coord["gene_number"]]
                    prev_canthyd_qualifier = utils.get_feature_qualifier(feature, 'canthyd')
                    if not prev_canthyd_qualifier:
                        if h["score"] > databases.CANTHYD_TRUSTED_CUTOFFS[h["hmm_id"]]:
                            confidence = 'high confidence'
                        else:
                            confidence = 'low confidence'
                        utils.set_feature_qualifier(feature, 'canthyd', f'{databases.CANTHYD_DESCR[h["hmm_id"]]}'
                                                                        f' ({h["hmm_id"]}) [{confidence}]')
                        subsystems.add_subsystem_to_feature(feature, '[hydrocarbon degradation]',
                                                            phrase=None, assignments=subsystem_hash)
