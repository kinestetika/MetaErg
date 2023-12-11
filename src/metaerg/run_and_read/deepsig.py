from pathlib import Path

from metaerg import context
from metaerg.datatypes import sqlite


def _run_programs(genome, contig_dict, db_connection, result_files):
    cds_aa_file = context.spawn_file('cds.faa', genome.name)
    context.run_external(f'deepsig -f {cds_aa_file} -o {result_files[1]} -k gramn')
    context.run_external(f'deepsig -f {cds_aa_file} -o {result_files[0]} -k gramp')


def _read_results(genome, contig_dict, db_connection, result_files) -> int:
    count = 0
    for file, source in zip(result_files, ('gram-positive', 'gram-negative')):
        with open(file) as reader:
            for line in reader:
                words = line.strip().split('\t')
                match words:
                    case [feature_id, _, 'Chain', start, end, score, _, _, _]:
                        pass
                    case [feature_id, _, 'Signal peptide', start, end, score, _, _, _]:
                        feature = sqlite.read_feature_by_id(db_connection, feature_id)
                        if not feature:
                            raise Exception(f'({genome.name}) Found deepsig result for unknown feature {feature_id}, '
                                            f'may need to rerun metaerg with --force')
                        if feature.signal_peptide:
                            feature.signal_peptide += ','
                        else:
                            count += 1
                        feature.signal_peptide += f'{source}({end})'
                        sqlite.update_feature_in_db(db_connection, feature)
                    case [feature_id, _, x, start, end, score, _, _, _]:
                        context.log(f'({genome.name}) Warning: Deepsig result file "{file}" contains unknown/new signal type "{x}"')
    return count


@context.register_annotator
def run_and_read_pureseqtm():
    return ({'pipeline_position': 122,
             'annotator_key': 'deepsig',
             'purpose': 'signal peptide prediction with DeepSig',
             'programs': ('deepsig',),
             'result_files': ('deepsig_gram-negative','deepsig_gram-positive'),
             'run': _run_programs,
             'read': _read_results})



# What follows next was used to test deepsig and compare results to signalp
# There are quite some differences, results only about 75% identical...
# This code can be deleted in the future


def load_data(file:Path, db_connection):
    count = 0
    data = {}
    with open(file) as reader:
        for line in reader:
            words = line.strip().split('\t')
            match words:
                case [feature_id, _, 'Chain', start, end, score, _, _, _]:
                    if feature_id not in data.keys():
                        data[feature_id] = '-'
                    continue
                case [feature_id, _, 'Signal peptide', start, end, score, _, _, _]:
                    # feature = sqlite.read_feature_by_id(db_connection, feature_id)
                    # if not feature:
                    #     raise Exception(f'Found results for unknown feature {feature_id}, '
                    #                     f'may need to rerun metaerg with --force')
                    # feature.signal_peptide = f'SP-{end}'
                    # sqlite.update_feature_in_db(db_connection, feature)
                    count += 1
                    data[feature_id] = 'SP'
                case [feature_id, _, x, start, end, score, _, _, _]:
                    print('found feature type ', x)
    print(file, count, len(data))
    return data




def load_signalp(file:Path):
    count = 0
    data = {}
    with open(file) as signalp_handle:
        for line in signalp_handle:
            if line.startswith("#"):
                continue
            words = line.split("\t")
            if "OTHER" == words[1]:
                continue
            feature_id = words[0].split()[0]
            data[feature_id] = words[1]
            count += 1
    print(file, count)
    return(data)


def main():
    deepsig_data_gramn = load_data(Path('/bio/fast/phormidium/temp/pha.gramn.deepsig'), 0)
    deepsig_data_gramp = load_data(Path('/bio/fast/phormidium/temp/pha.gramp.deepsig'), 0)
    signalp_data = load_signalp(Path('/bio/fast/phormidium/temp/pha.signalp/prediction_results.txt'))
    gramp_gramn = 0
    gramp_signalp = 0
    gramn_signalp = 0
    signalp = 0
    gramn = 0
    gramp = 0
    all = 0
    for i in deepsig_data_gramn:
        if deepsig_data_gramn[i] =='SP' or deepsig_data_gramp[i] =='SP' or (i in signalp_data.keys()):
            print(i, deepsig_data_gramn[i], deepsig_data_gramp[i], signalp_data.get(i, '-'))
        if deepsig_data_gramn[i] =='SP' and deepsig_data_gramp[i] =='SP' and (i in signalp_data.keys()):
            all += 1
        elif deepsig_data_gramn[i] =='SP' and deepsig_data_gramp[i] =='SP':
            gramp_gramn += 1
        elif deepsig_data_gramp[i] =='SP' and (i in signalp_data.keys()):
            gramp_signalp += 1
        elif deepsig_data_gramn[i] =='SP' and (i in signalp_data.keys()):
            gramn_signalp += 1
        elif deepsig_data_gramn[i] =='SP':
            gramn += 1
        elif deepsig_data_gramp[i] =='SP':
            gramp += 1
        elif i in signalp_data.keys():
            signalp += 1
    print('gramn', gramn)
    print('gramp', gramp)
    print('signalp', signalp)
    print('gramn + gramp', gramp_gramn)
    print('gramp + signalp', gramp_signalp)
    print('gramn + signalp', gramn_signalp)
    print('all', all)


if __name__ == "__main__":
    main()
