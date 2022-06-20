import shutil
from metaerg import context
from metaerg.data_model import Genome


def _run_programs(genome:Genome, result_files):
    cds_aa_file = context.spawn_file('cds.faa', genome.id)
    with open(result_files[0], 'w') as output, open(cds_aa_file) as input:
        context.run_external('tmhmm', stdin=input, stdout=output)
    # this is not thread-safe:
    for file in result_files[0].parent.glob(f'TMHMM_*'):
        if file.is_dir():
            shutil.rmtree(file)


def _read_results(genome:Genome, result_files) -> int:
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
                case [feature_name, _, orientation, _, _] if orientation in ('inside', 'outside'):
                    if not current_feature or current_feature.id != feature_name:
                        if feature_tmh_count:
                            current_feature.transmembrane_helixes = " ".join((str(feature_tmh_count), current_txt[:-1]))
                            count += 1
                        new_feature = genome.get_feature(feature_name)
                        current_feature = new_feature
                        feature_tmh_count = 0
                        current_txt = 'i,' if 'inside' == orientation else 'o,'
        if feature_tmh_count:
            current_feature.transmembrane_helixes = " ".join((feature_tmh_count, current_txt[:-1]))
            count += 1
    return count


@context.register_annotator
def run_and_read_tmhmm():
    return ({'pipeline_position': 111,
             'purpose': 'transmembrane helix prediction with tmhmm',
             'programs': ('tmhmm',),
             'result_files': ('tmhmm',),
             'run': _run_programs,
             'read': _read_results})
