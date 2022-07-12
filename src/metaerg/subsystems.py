import re
import pandas as pd

from metaerg import subsystems_data
from metaerg.datatypes.blast import BlastResult, BlastHit

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

current_subsystem = None
subsystems = []
index_data = []
for line in subsystems_data.subsystem_data().split('\n'):
    line = line.strip()
    if line.startswith("#") or not len(line):
        continue
    elif line.startswith(">"):
        current_subsystem = line[1:]
    elif current_subsystem is not None:
        si = line.find(' ')
        subsystems.append({'profiles': line[:si], 'genes': ''})
        index_data.append((current_subsystem, line[si+1:]))
DATAFRAME_INDEX = pd.MultiIndex.from_tuples(index_data, names=['subsystem', 'function'])
SUBSYSTEM_DATA = pd.DataFrame(subsystems, dtype=str, index=DATAFRAME_INDEX)


def match(blast_result:BlastResult) -> str:
    for h in blast_result.hits:
        if subsystems := match_hit(h):
            return subsystems
    return ''


def match_hit(blast_hit:BlastHit) -> str:
    matching_subsystems = set()
    for sf in SUBSYSTEM_DATA.itertuples():
        for profile in sf.profiles.split('|'):
            if profile == blast_hit.hit.accession and blast_hit.aligned_length >= 0.7 * blast_hit.hit.length:
                matching_subsystems.add(f'[{sf.Index[0]}|{sf.Index[1]}]')
    return ' '.join(matching_subsystems).strip()


def aggregate(feature_data: pd.DataFrame):
    aggregated_subsystem_data = SUBSYSTEM_DATA.copy()
    unstructured_subsystems = {}
    for feature in feature_data.itertuples():
        if not len(feature.subsystems):
            continue
        for s in re.split(r'\s(?=\[)', feature.subsystems):  # split at space if followed b y [
            s = s[1:-1]
            if '|' in s:
                gene_subsystem, gene_function = s.split('|')
                try:
                    aggregated_subsystem_data.at[(gene_subsystem, gene_function), 'genes'] += f' {feature.id}'
                except TypeError:
                    aggregated_subsystem_data.at[(gene_subsystem, gene_function), 'genes'] = feature.id
            else:  # this is a subsystem with undefined structure, add a row to the dataframe
                try:
                    unstructured_subsystems[s] += f' {feature.id}'
                except KeyError:
                    unstructured_subsystems[s] = feature.id
    new_dataframe_rows = [{'profiles': '', 'genes': v} for k, v in unstructured_subsystems.items()]
    new_dataframe_index = [(k, k) for  k, v in unstructured_subsystems.items()]
    new_dataframe_index = pd.MultiIndex.from_tuples(new_dataframe_index, names=['subsystem', 'function'])
    new_dataframe_data = pd.DataFrame(new_dataframe_rows, dtype=str, index=new_dataframe_index)
    aggregated_subsystem_data = pd.concat([aggregated_subsystem_data, new_dataframe_data])
    # print (aggregated_subsystem_data)
    return aggregated_subsystem_data

