import re
import pandas as pd

from metaerg import subsystems_data
from metaerg.datatypes.blast import BlastResult, BlastHit

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

