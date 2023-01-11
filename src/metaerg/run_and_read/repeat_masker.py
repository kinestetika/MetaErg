import shutil
from pathlib import Path

from metaerg import context
from metaerg.datatypes import fasta
from metaerg.datatypes import sqlite

def _run_programs(genome_name, contig_dict, db_connection, result_files):
    fasta_file = context.spawn_file('masked', genome_name)
    fasta.write_contigs_to_fasta(contig_dict, fasta_file, db_connection, genome_name,
                                 mask_targets=fasta.ALL_MASK_TARGETS)
    lmer_table_file = context.spawn_file('lmer-table', genome_name)
    repeatscout_file_raw = context.spawn_file('repeatscout-raw', genome_name)
    repeatscout_file_filtered = context.spawn_file('repeatscout-filtered', genome_name)
    repeatscout_file_filtered2 = context.spawn_file('repeatscout-filtered2', genome_name)

    context.run_external(f'build_lmer_table -sequence {fasta_file} -freq {lmer_table_file}')
    context.run_external(f'RepeatScout -sequence {fasta_file} -output {repeatscout_file_raw} -freq {lmer_table_file}')
    with open(repeatscout_file_filtered, 'w') as output, open(repeatscout_file_raw) as input:
        context.run_external('filter-stage-1.prl', stdin=input, stdout=output)

    repeatmasker_output_file = fasta_file.parent / f'{fasta_file.name}.out'  # nothing we can do about that

    if not repeatscout_file_filtered.stat().st_size:
        context.log(f'({genome_name}) No repeats remaining after filter stage 1.')
        finalize(False, genome_name, fasta_file, repeatmasker_output_file, result_files[0])
        return

    context.run_external(f'RepeatMasker -pa {context.CPUS_PER_GENOME} -lib {repeatscout_file_filtered} '
                         f'-dir {fasta_file.parent} {fasta_file}')
    with open(repeatscout_file_filtered2, 'w') as output, open(repeatscout_file_filtered) as input:
        context.run_external(f'filter-stage-2.prl --cat={repeatmasker_output_file} --thresh=10', stdin=input, stdout=output)

    if not repeatscout_file_filtered2.stat().st_size:
        context.log(f'({genome_name}) No repeats remaining after filter stage 2.')
        finalize(False, genome_name, fasta_file, repeatmasker_output_file, result_files[0])
        return

    context.run_external(f'RepeatMasker -pa {context.CPUS_PER_GENOME} -lib {repeatscout_file_filtered2} '
                         f'-dir {fasta_file.parent} {fasta_file}')
    finalize(True, genome_name, fasta_file, repeatmasker_output_file, result_files[0])


def finalize(success: bool, genome_name: str, fasta_input_file: Path, repeatmasker_output_file: Path,
             final_repeatmasker_output_file: Path):
    if not success:
        with open(repeatmasker_output_file, 'w') as handle:
            handle.write('#No repeats detected by repeatmasker')
    shutil.move(repeatmasker_output_file, final_repeatmasker_output_file)
    for file in fasta_input_file.parent.glob(f'{fasta_input_file.name}.*'):
        if file.is_dir():
            shutil.rmtree(file)
        else:
            file.unlink()


def words2feature(words: list[str], contig, genome_name:str) -> sqlite.Feature:
    start = int(words[5]) - 1
    end = int(words[6])
    strand = -1 if 'C' == words[8] else 1
    seq = contig['seq'][start:end]
    if strand < 0:
        seq = fasta.reverse_complement(seq)
    return sqlite.Feature(genome = genome_name,
                          contig = contig['id'],
                          start = start,
                          end = end,
                          strand = strand,
                          type = 'repeat',
                          inference = 'repeatmasker',
                          nt_seq = seq)


def _read_results(genome_name, contig_dict, db_connection, result_files) -> int:
    """(1) simple repeats, these are consecutive
       (2) unspecified repeats, these occur scattered and are identified by an id in words[9]. We only
           add those when they occur 10 or more times."""
    count = 0
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
                feature.notes = f'repeat {words[9]}'
                sqlite.add_new_feature_to_db(db_connection, feature)
                count += 1
            else:
                repeat_list = repeat_hash.setdefault(words[9], list())
                repeat_list.append(words2feature(words, contig, genome_name))
    for repeat_list in repeat_hash.values():
        if len(repeat_list) >= 10:
            for feature in repeat_list:
                feature.notes = f' (occurs {len(repeat_list)}x)'
                sqlite.add_new_feature_to_db(db_connection, feature)
                count += 1
    return count


@context.register_annotator
def run_and_read_repeatmasker():
    return ({'pipeline_position': 51,
             'annotator_key': 'repeat_masker',
             'purpose': 'repeat prediction with repeatmasker',
             'programs': ('build_lmer_table', 'RepeatScout', 'filter-stage-1.prl', 'RepeatMasker'),
             'result_files': ('repeatmasker',),
             'run': _run_programs,
             'read': _read_results})
