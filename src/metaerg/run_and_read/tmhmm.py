import shutil

from metaerg import context
from metaerg.datatypes import sqlite

def _run_programs(genome_name, contig_dict, db_connection, result_files):
    cds_aa_file = context.spawn_file('cds.faa', genome_name)
    with open(result_files[0], 'w') as output, open(cds_aa_file) as input:
        context.run_external('tmhmm', stdin=input, stdout=output)


def _read_results(genome_name, contig_dict, db_connection, result_files) -> int:
    count = 0
    current_feature = None
    current_txt = ""
    feature_tmh_count = 0
    with open(result_files[0]) as tmhmm_handle:
        for line in tmhmm_handle:
            words = line.strip().split()
            match words:
                case[first_text, *_] if first_text.startswith('#'):
                    continue
                case [_, _, 'TMhelix', start, end]:
                    feature_tmh_count += 1
                    current_txt += f'{start}-{end},'
                case [next_feature_name, _, orientation, _, _] if orientation in ('inside', 'outside'):
                    next_feature = sqlite.read_feature_by_id(db_connection, next_feature_name)
                    if not next_feature:
                        raise Exception(f'Found results for unknown feature {next_feature_name}, '
                                        f'may need to rerun metaerg with --force')
                    if not current_feature or next_feature.id != current_feature.id:
                        if feature_tmh_count:
                            current_feature.tmh = feature_tmh_count
                            current_feature.tmh_topology = current_txt[:-1]
                            sqlite.update_feature_in_db(db_connection, current_feature)
                            count += 1
                        current_feature = next_feature
                        feature_tmh_count = 0
                        current_txt = 'i,' if 'inside' == orientation else 'o,'
        if feature_tmh_count:
            current_feature.tmh = feature_tmh_count
            current_feature.tmh_topology = current_txt[:-1]
            sqlite.update_feature_in_db(db_connection, current_feature)
            count += 1
    return count


@context.register_annotator
def run_and_read_tmhmm():
    return ({'pipeline_position': 111,
             'annotator_key': 'tmhmm',
             'purpose': 'transmembrane helix prediction with tmhmm',
             'programs': ('tmhmm',),
             'result_files': ('tmhmm',),
             'run': _run_programs,
             'read': _read_results})


def cleanup(dir):
    for file in dir.glob(f'TMHMM_*'):
        if file.is_dir():
            shutil.rmtree(file)
