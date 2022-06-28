import re
import pandas as pd

import subsystems_data

DATAFRAME_COLUMNS = 'function subsystem genes'.split()
SUBSYSTEM_DATA = pd.DataFrame(columns=DATAFRAME_COLUMNS, dtype=str)

current_subsystem = None
for line in subsystems_data.subsystem_data().split('\n'):
    subsystems = []
    line = line.strip()
    if line.startswith("#") or not len(line):
        continue
    elif line.startswith(">"):
        current_subsystem = line[1:]
    elif current_subsystem is not None:
        subsystems.append({'subsystem': current_subsystem, 'function': line})
    SUBSYSTEM_DATA = pd.concat([SUBSYSTEM_DATA, pd.DataFrame(subsystems)], ignore_index=True)
    SUBSYSTEM_DATA.set_index('function', drop=False, inplace=True)


def match(descriptions) -> str:
    for d in descriptions:
        for sf in SUBSYSTEM_DATA.itertuples():
            if len(d) > len(sf.function) + 20:
                continue
            match = re.search(r'\b' + sf.function + r'\b', d)
            if match and match.start() < 10:
                return f'[{sf.subsystem}:{sf.function}]'
    return ''


def aggregate(feature_data: pd.DataFrame):
    aggregated_subsystem_data = SUBSYSTEM_DATA.copy()
    unstructured_subsystems = {}
    for feature in feature_data.itertuples():
        if not len(feature.subsystems):
            continue
        for s in re.split(r'\s(?=\[)', feature.subsystems):  # split at space if followed b y [
            s = s[1:-1]
            if ':' in s:
                gene_subsystem, gene_function = s.split(':')
                try:
                    aggregated_subsystem_data.at[gene_function, 'genes'] += f' {feature.id}'
                except TypeError:
                    aggregated_subsystem_data.at[gene_function, 'genes'] = feature.id
            else:  # this is a subsystem with undefined structure, add a row to the dataframe
                try:
                    unstructured_subsystems[s] += f' {feature.id}'
                except KeyError:
                    unstructured_subsystems[s] = feature.id
    new_dataframe_rows = [{'function': k, 'subsystem': k, 'genes': v} for k, v in unstructured_subsystems.items()]
    aggregated_subsystem_data = pd.concat([aggregated_subsystem_data, pd.DataFrame(new_dataframe_rows)])
    aggregated_subsystem_data.set_index('function', drop=False, inplace=True)
    return aggregated_subsystem_data

