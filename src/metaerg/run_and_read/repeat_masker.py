import shutil
import pandas as pd
from pathlib import Path

from metaerg import context
from metaerg.datatypes import fasta


def _run_programs(genome_name, contig_dict, feature_data: pd.DataFrame, result_files):
    fasta_file = context.spawn_file('masked', genome_name)
    fasta.write_contigs_to_fasta(contig_dict, fasta_file, feature_data, genome_name,
                                 mask_targets=fasta.ALL_MASK_TARGETS)
    lmer_table_file = context.spawn_file('lmer-table', genome_name)
    repeatscout_file_raw = context.spawn_file('repeatscout-raw', genome_name)
    repeatscout_file_filtered = context.spawn_file('repeatscout-filtered', genome_name)

    context.run_external(f'build_lmer_table -sequence {fasta_file} -freq {lmer_table_file}')
    context.run_external(f'RepeatScout -sequence {fasta_file} -output {repeatscout_file_raw} -freq {lmer_table_file}')
    with open(repeatscout_file_filtered, 'w') as output, open(repeatscout_file_raw) as input:
        context.run_external('filter-stage-1.prl', stdin=input, stdout=output)
    repeatmasker_output_file = Path(f'{fasta_file.name}.out')  # nothing we can do about that
    if repeatscout_file_filtered.stat().st_size > 0:
        context.run_external(f'RepeatMasker -pa {context.CPUS_PER_GENOME} -lib {repeatscout_file_filtered} -dir . '
                             f'{fasta_file}')
    else:
        context.log(f'({genome_name}) No repeats detected by repeatmasker.')
        with open(repeatmasker_output_file, 'w') as handle:
            handle.write('#No repeats detected by repeatmasker')
    shutil.move(repeatmasker_output_file, result_files[0])
    for file in Path.cwd().glob(f'{fasta_file.name}.*'):
        if file.is_dir():
            shutil.rmtree(file)
        else:
            file.unlink()


def words2feature(words: list[str], contig, genome_name:str):
    start = int(words[5]) - 1
    end = int(words[6])
    strand = -1 if 'C' == words[8] else 1
    seq = contig['seq'][start:end]
    if strand < 0:
        seq = fasta.reverse_complement(seq)
    return {'genome': genome_name,
            'contig': contig['id'],
            'start': start,
            'end': end,
            'strand': strand,
            'type': 'repeat',
            'inference': 'repeatmasker',
            'seq': seq}


def _read_results(genome_name, contig_dict, feature_data: pd.DataFrame, result_files) -> tuple:
    """(1) simple repeats, these are consecutive
       (2) unspecified repeats, these occur scattered and are identified by an id in words[9]. We only
           add those when they occur 10 or more times."""
    new_features = []
    repeat_hash = dict()
    with open(result_files[0]) as repeatmasker_handle:
        for line in repeatmasker_handle:
            words = line.split()
            if len(words) < 11 or words[0] in ('SW', 'score'):
                continue
            try:
                contig = contig_dict[words[4]]
            except KeyError:
                context.log(f'({genome_name}) Warning: Unknown contig id "{words[4]}"')
                continue
            if 'Simple_repeat' == words[10]:
                feature = words2feature(words, contig, genome_name)
                new_features.append(feature)
                feature['notes'] = f'repeat {words[9]}'
            else:
                repeat_list = repeat_hash.setdefault(words[9], list())
                repeat_list.append(words2feature(words, contig, genome_name))
    for repeat_list in repeat_hash.values():
        if len(repeat_list) >= 10:
            for feature in repeat_list:
                new_features.append(feature)
                feature['notes'] = f' (occurs {len(repeat_list)}x)'
    feature_data = pd.concat([feature_data, pd.DataFrame(new_features)], ignore_index=True)
    return feature_data, len(new_features)


@context.register_annotator
def run_and_read_repeatmasker():
    return ({'pipeline_position': 51,
             'purpose': 'repeat prediction with repeatmasker',
             'programs': ('build_lmer_table', 'RepeatScout', 'filter-stage-1.prl', 'RepeatMasker'),
             'result_files': ('repeatmasker',),
             'run': _run_programs,
             'read': _read_results})
