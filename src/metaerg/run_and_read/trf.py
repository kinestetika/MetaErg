import pandas as pd
from metaerg import context
from metaerg.datatypes import fasta


def _run_programs(genome_name, contig_dict, feature_data: pd.DataFrame, result_files):
    fasta_file = context.spawn_file('masked', genome_name)
    fasta.write_contigs_to_fasta(contig_dict, fasta_file, feature_data, genome_name,
                                 mask_targets=fasta.ALL_MASK_TARGETS)
    with open(result_files[0], 'w') as output:
        context.run_external(f'trf {fasta_file} 2 7 7 80 10 50 500 -d -h -ngs', stdout=output)


def _read_results(genome_name, contig_dict, feature_data: pd.DataFrame, result_files) -> tuple:
    new_features = []
    with open(result_files[0]) as trf_handle:
        for line in trf_handle:
            if line.startswith("@"):
                contig = contig_dict[line[1:].strip()]
                continue
            if not contig:
                continue
            words = line.split()
            start = int(words[0]) - 1
            end = int(words[1])
            seq = contig['seq'][start:end]
            feature = {'genome': genome_name,
                       'contig': contig['id'],
                       'start': start,
                       'end': end,
                       'strand': 1,
                       'type': 'repeat',
                       'inference': 'tandem-repeat-finder',
                       'seq': seq,
                       'notes': f'period size {words[2]}; copies {words[3]}'}
            new_features.append(feature)
    feature_data = pd.concat([feature_data, pd.DataFrame(new_features)], ignore_index=True)
    return feature_data, len(new_features)


@context.register_annotator
def run_and_read_trf():
    return ({'pipeline_position': 41,
             'purpose': 'tandem repeat prediction with trf',
             'programs': ('trf',),
             'result_files': ('tandem-repeat-finder',),
             'run': _run_programs,
             'read': _read_results})
