from pathlib import Path

from metaerg import context
from metaerg.datatypes import sqlite
from metaerg.datatypes.fasta import FastaParser

def _run_programs(genome, contig_dict, db_connection, result_files):
    cds_aa_file = context.spawn_file('cds.faa', genome.name)
    context.run_external(f'PureseqTM_proteome.sh -i {cds_aa_file} -o {result_files[0]} -c {context.CPUS_PER_GENOME}')


def _read_results(genome, contig_dict, db_connection, result_files) -> int:
    count = 0
    with FastaParser(result_files[0], cleanup_seq=False) as fasta_reader:
        for orf in fasta_reader:
            tmh_list = []
            orf['seq'] = orf['seq'][len(orf['seq'])//2:]  # the first half is the aa seq
            if orf['seq'][0] == '1':
                current_tmh = {'start': 0}
            else:
                current_tmh = None
            for pos in range(1,len(orf['seq'])):
                if '0' == orf['seq'][pos-1] and '1' == orf['seq'][pos]:
                    current_tmh = {'start': pos}
                if '1' == orf['seq'][pos-1] and '0' == orf['seq'][pos]:
                    current_tmh['end'] = pos
                    tmh_list.append(current_tmh)
                    current_tmh = None
            if current_tmh:
                current_tmh['end'] = len(orf['seq']) - 1
                tmh_list.append(current_tmh)
            if len(tmh_list):
                feature = sqlite.read_feature_by_id(db_connection, orf['id'])
                if not feature:
                    raise Exception(f'({genome.name}) Found pureseqtm result for unknown feature {orf["id"]}, '
                                    f'may need to rerun metaerg with --force')
                feature.tmh = len(tmh_list)
                feature.tmh_topology = ','.join(f'{tmh["start"]}-{tmh["end"]}' for tmh in tmh_list)
                sqlite.update_feature_in_db(db_connection, feature)
                count += 1
    return count


@context.register_annotator
def run_and_read_pureseqtm():
    return ({'pipeline_position': 112,
             'annotator_key': 'pureseqtm',
             'purpose': 'transmembrane helix prediction with PureseqTM',
             'programs': ('PureseqTM_proteome.sh',),
             'result_files': ('pureseqtm',),
             'run': _run_programs,
             'read': _read_results})


# What follows next was used to test pureseqtm and compare results to TMHMM
# The results were 95% identical...
# This code can be deleted in the future

def load_tm_helixe_info_from_file(file:Path):
    count = 0
    with FastaParser(file, cleanup_seq=False) as fasta_reader:
        orfs = {}
        for orf in fasta_reader:
            tmh_list = []
            orf['seq'] = orf['seq'][len(orf['seq'])//2:]
            if orf['seq'][0] == '1':
                current_tmh = {'start': 0}
            else:
                current_tmh = None
            for pos in range(1,len(orf['seq'])):
                if '0' == orf['seq'][pos-1] and '1' == orf['seq'][pos]:
                    current_tmh = {'start': pos}
                if '1' == orf['seq'][pos-1] and '0' == orf['seq'][pos]:
                    current_tmh['end'] = pos
                    #print(f'{orf["id"]}: {current_tmh["start"]}-{current_tmh["end"]}' )
                    tmh_list.append(current_tmh)
                    current_tmh = None
                    count += 1
            if current_tmh:
                current_tmh['end'] = len(orf['seq']) - 1
                tmh_list.append(current_tmh)
                count += 1
            orfs[orf['id']] = tmh_list
            #if len(tmh_list) > 6:
            #    print(tmh_list)
        print(f'{file}: parsed {count} TMH.')
    return orfs

def determine_overlap(results1, results2):
    total_overlap = 0
    total_1 = 0
    total_2 = 0
    for orf_id, tmh_list1 in results1.items():
        #print(tmh_list1)
        tmh_list2 = results2.get(orf_id, [])
#        if len(tmh_list1) + len(tmh_list2):
#            print(orf_id)
#            print('1: ' + ', '.join(f'{tmh["start"]}-{tmh["end"]}' for tmh in tmh_list1))
#            print('2: ' + ', '.join(f'{tmh["start"]}-{tmh["end"]}' for tmh in tmh_list2))
        overlap = 0
        for tmh1 in tmh_list1:
            for tmh2 in tmh_list2:
                if tmh1['end'] > tmh2['start'] and tmh1['start'] < tmh2['end']:
                    overlap += 1
                    break
        total_overlap += overlap
        total_1 += len(tmh_list1)
        total_2 += len(tmh_list2)
    print(total_overlap, total_1, total_2, f'{total_overlap / min(total_1, total_2) * 100: .1f}')

def load_tmhmm_results(file:Path):
    count = 0
    current_orf_name = ""
    orfs = {}
    tmh_list = []
    with open(file) as tmhmm_handle:
        for line in tmhmm_handle:
            words = line.strip().split()
            match words:
                case[first_text, *_] if first_text.startswith('#'):
                    continue
                case [_, _, 'TMhelix', start, end]:
                    tmh_list.append({'start': int(start),
                                     'end': int(end)})
                    count += 1
                case [next_orf_name, _, orientation, _, _] if orientation in ('inside', 'outside'):
                    if current_orf_name and next_orf_name != current_orf_name:
                        #print(current_orf_name, next_orf_name)
                        orfs[current_orf_name] = tmh_list
                        #print(tmh_list)
                        tmh_list = []
                    current_orf_name = next_orf_name
        if current_orf_name:
            orfs[current_orf_name] = tmh_list
    print(f'{file}: parsed {count} TMH.')
    return orfs


def main():
    #results_fast_sp = load_tm_helixe_info_from_file(Path('/bio/bin/Pureseq_test_pha_fast_sp'))
    results_fast = load_tm_helixe_info_from_file(Path('/bio/bin/Pureseq_test_pha_fast'))
    results_slow_sp = load_tm_helixe_info_from_file(Path('/bio/bin/Pureseq_test_pha_sp'))
    results_tmhmm = load_tmhmm_results(Path('/bio/fast/phormidium/temp/pha.tmhmm'))

    determine_overlap(results_fast, results_tmhmm)
    determine_overlap(results_slow_sp, results_tmhmm)

if __name__ == "__main__":
    main()

