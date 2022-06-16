import shutil
from pathlib import Path

from metaerg.data_model import Genome, SeqRecord, SeqFeature, FeatureType
from metaerg import context
from metaerg import bioparsers


def _run_programs(genome:Genome, result_files):
    fasta_file, = bioparsers.write_genome_to_fasta_files(genome, context.spawn_file(genome, 'masked', genome.id),
                                                         mask=True)
    lmer_table_file = context.spawn_file('lmer-table', genome.id)
    repeatscout_file_raw = context.spawn_file('repeatscout-raw', genome.id)
    repeatscout_file_filtered = context.spawn_file('repeatscout-filtered', genome.id)

    context.run_external(f'build_lmer_table -sequence {fasta_file} -freq {lmer_table_file}')
    context.run_external(f'RepeatScout -sequence {fasta_file} -output {repeatscout_file_raw} -freq {lmer_table_file}')
    with open(repeatscout_file_filtered, 'w') as output, open(repeatscout_file_raw) as input:
        context.run_external('filter-stage-1.prl', stdin=input, stdout=output)
    context.run_external(f'RepeatMasker -pa {context.CPUS_PER_GENOME} -lib {repeatscout_file_filtered} -dir . '
                         f'{fasta_file}')
    repeatmasker_output_file = Path(f'{fasta_file.name}.out')  # nothing we can do about that
    shutil.move(repeatmasker_output_file, result_files[0])
    for file in Path.cwd().glob(f'{fasta_file.name}.*'):
        if file.is_dir():
            shutil.rmtree(file)
        else:
            file.unlink()


def words2feature(words: list[str], contig: SeqRecord) -> SeqFeature:
    start = int(words[5]) - 1
    end = int(words[6])
    strand = -1 if 'C' == words[8] else 1
    seq = contig.seq[start:end]
    if strand < 0:
        seq = bioparsers.reverse_complement(seq)
    return SeqFeature(start, end, strand, FeatureType.repeat, seq=seq, inference='repeatmasker')


def _read_results(genome:Genome, result_files) -> int:
    """(1) simple repeats, these are consecutive
       (2) unspecified repeats, these occur scattered and are identified by an id in words[9]. We only
           add those when they occur 10 or more times."""
    repeat_count = 0
    repeat_hash = dict()
    with open(result_files[0]) as repeatmasker_handle:
        for line in repeatmasker_handle:
            words = line.split()
            if len(words) < 11:
                continue
            contig: SeqRecord = genome.contigs[words[4]]
            if 'Simple_repeat' == words[10]:
                repeat_count += 1
                feature = words2feature(words, contig)
                contig.features.append(feature)
                feature.notes.add(f'repeat {words[9]}')
            else:
                repeat_list = repeat_hash.setdefault(words[9], list())
                repeat_list.append((contig, words2feature(words, contig)))
    for repeat_list in repeat_hash.values():
        if len(repeat_list) >= 10:
            for c, f in repeat_list:
                repeat_count += 1
                c.features.append(f)
                f.notes.add(f' (occurs {len(repeat_list)}x)')
    return repeat_count


@context.register_annotator
def run_and_read_repeatmasker():
    return ({'pipeline_position': 51,
             'purpose': 'repeat prediction with repeatmasker',
             'programs': ('build_lmer_table', 'RepeatScout', 'filter-stage-1.prl', 'RepeatMasker'),
             'result_files': ('repeatmasker',),
             'run': _run_programs,
             'read': _read_results})
