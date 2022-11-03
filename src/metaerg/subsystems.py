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
    subsystems = ''
    i = 0
    for h in blast_result.hits:
        i += 1
        if i > 5:
            break
        if new_subsystem := match_hit(h):
            subsystems += ' ' + new_subsystem
    return subsystems


def match_hit(blast_hit:BlastHit, confidence=1) -> str:
    matching_subsystems = set()
    for sf in SUBSYSTEM_DATA.itertuples():
        for profile in sf.profiles.split('|'):
            if profile == blast_hit.hit.accession and blast_hit.aligned_length >= 0.7 * blast_hit.hit.length:
                matching_subsystems.add(f'[{sf.Index[0]}|{sf.Index[1]}|{confidence:.1f}]')
    return ' '.join(matching_subsystems).strip()


def cleanup_subsystem_str(subsystem_str: str) -> str:
    subsystem_str = subsystem_str.strip()
    collected_subsystems = {}
    for s in re.split(r'\s(?=\[)', subsystem_str):  # split at space if followed by [
        s = s.strip()[1:-1]
        if '|' in s:
            gene_subsystem, gene_function, confidence = s.split('|')
            confidence = float(confidence)
            hash_str = f'{gene_subsystem}|{gene_function}'
            try:
                prev_confidence = collected_subsystems[hash_str]
                if confidence > prev_confidence:
                    collected_subsystems[hash_str] = confidence
            except KeyError:
                collected_subsystems[hash_str] = confidence
        else:
            collected_subsystems[s] = 0
    subsystem_str = ' '.join((f'[{val[0]}|{val[1]:.1f}]' for val in sorted(collected_subsystems.items(), key=lambda x: (-x[1], x[0]))))
    return subsystem_str.strip()


def aggregate(feature_data: pd.DataFrame):
    aggregated_subsystem_data = SUBSYSTEM_DATA.copy().sort_index()
    unstructured_subsystems = {}
    for feature in feature_data.itertuples():
        if not len(feature.subsystems):
            continue
        for s in re.split(r'\s(?=\[)', feature.subsystems):  # split at space if followed b y [
            s = s[1:-1]
            try:
                gene_subsystem, gene_function, confidence = s.split('|')
                try:
                    aggregated_subsystem_data.loc[(gene_subsystem, gene_function), 'genes'] += \
                        f' {feature.id}@{confidence}'
                except TypeError:
                    aggregated_subsystem_data.loc[(gene_subsystem, gene_function), 'genes'] = feature.id
            except ValueError:
                gene_subsystems = s.split('|')
                try:
                    unstructured_subsystems[gene_subsystems[0]] += f' {feature.id}'
                except KeyError:
                    unstructured_subsystems[gene_subsystems[0]] = feature.id
    new_dataframe_rows = [{'profiles': '', 'genes': v} for k, v in unstructured_subsystems.items()]
    new_dataframe_index = [(k, k) for  k, v in unstructured_subsystems.items()]
    new_dataframe_index = pd.MultiIndex.from_tuples(new_dataframe_index, names=['subsystem', 'function'])
    new_dataframe_data = pd.DataFrame(new_dataframe_rows, dtype=str, index=new_dataframe_index)
    aggregated_subsystem_data = pd.concat([aggregated_subsystem_data, new_dataframe_data])
    # print (aggregated_subsystem_data)
    return aggregated_subsystem_data
