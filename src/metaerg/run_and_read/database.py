from pathlib import Path
from metaerg import utils
from collections import namedtuple

CDD_INDEX_FILENAME = "cddid.tbl"

CDDEntry = namedtuple('CDDEnry', ['name', 'gene', 'descr', 'length'])

class MetaergDatabase:
    def __init__(self, database_dir: Path):
        self.database_dir = database_dir
        self.descriptions = {}
        self.taxonomy = {}
        self.cdd = {}
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
        self._load()

    def __repr__(self):
        return 'Database(database_dir={})'.format(self.database_dir)

    def get_diamond_entry(self, id):
        words = id.split('~')  # org_acc gene_acc [pev] gene# decr# taxon#
        return DiamondEntry(self.taxonomy[words[2]][int(words[5])], self.descriptions[words[2]][int(words[4])],
                            int(words[6]), int(words[3]))

    def is_dabase_ready(self):
        return len(self.cdd) and len(self.descriptions) and len(self.taxonomy)

    def write_database(self, folder:Path):
        with open(Path(folder, DB_DESCR_FILENAME), 'w') as file_handle:
            for kingdom in self.descriptions:
                for db_id in self.descriptions[kingdom]:
                    file_handle.write(f'{kingdom}\t{db_id}\t{self.descriptions[kingdom][db_id]}\n')
        with open(Path(dir, DB_TAXON_FILENAME), 'w') as file_handle:
            for kingdom in TAXONOMY_CACHE:
                for db_id in TAXONOMY_CACHE[kingdom]:
                    file_handle.write(f'{kingdom}\t{db_id}\t{TAXONOMY[kingdom][db_id]}\n')
        with open(Path(dir, CDD_INDEX_FILENAME), 'w') as file_handle:
            for cdd_id in CDD_CACHE:
                cdd_item = CDD[cdd_id]
                file_handle.write(f'{cdd_id}\t{cdd_item[0]}\t{cdd_item[1]}\t{cdd_item[2]}\t{cdd_item[3]}\n')

    def _load(self):
        # load cdd
        cdd_file = Path(self.database_dir, CDD_INDEX_FILENAME)
        if cdd_file.exists() and cdd_file.stat().st_size:
            with open(cdd_file) as cdd_descr_handle:
                for line in cdd_descr_handle:
                    words = line.split("\t")
                    # id, id_name, gene_name, descr, length
                    self.cdd[int(words[0])] = CDDEntry(words[1], words[2], words[3], int(words[4]))
        # load canthyd
        canthyd_file = Path(self.database_dir, 'canthyd', 'CANT-HYD.hmm')
        if canthyd_file.exists() and canthyd_file.stat().st_size:
            current_name = None
            with open(canthyd_file) as handle:
                for line in handle:
                    if line.startswith('NAME'):
                        current_name = line.split()[1]
                    elif line.startswith('TC'):
                        self.canthyd_trusted_cutoffs[current_name] = int(line.split()[1])
