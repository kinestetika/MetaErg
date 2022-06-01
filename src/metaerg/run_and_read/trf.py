from metaerg.data_model import MetaergGenome, MetaergSeqRecord, FeatureType
from metaerg import context


def _run_programs(genome:MetaergGenome, result_files):
    fasta_file, = genome.write_fasta_files(context.spawn_file('masked', genome.id), masked=True)
    with open(result_files[0], 'w') as output:
        context.run_external(f'trf {fasta_file} 2 7 7 80 10 50 500 -d -h -ngs', stdout=output)


def _read_results(genome:MetaergGenome, result_files) -> int:
    tr_count = 0
    with open(result_files[0]) as trf_handle:
        for line in trf_handle:
            if line.startswith("@"):
                contig: MetaergSeqRecord = genome.contigs[line[1:].strip()]
                continue
            if not contig:
                continue
            words = line.split()
            f = contig.spawn_feature(int(words[0]) - 1, int(words[1]), 1, FeatureType.repeat,
                                     inference='tandem-repeat-finder')
            f.notes.add(f'period size {words[2]}; copies {words[3]}')
            tr_count += 1
    return tr_count


@context.register_annotator
def run_and_read_trf():
    return ({'pipeline_position': 41,
             'purpose': 'tandem repeat prediction with trf',
             'programs': ('trf'),
             'result_files': ('tandem-repeat-finder',),
             'run': _run_programs,
             'read': _read_results})
