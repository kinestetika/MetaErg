from metaerg import context
from metaerg.datatypes import fasta
from metaerg.datatypes import sqlite


def _run_programs(genome, contig_dict, db_connection, result_files):
    fasta_file = context.spawn_file('masked', genome.name)
    fasta.write_contigs_to_fasta(contig_dict, fasta_file, db_connection, genome.name,
                                 mask_targets=fasta.ALL_MASK_TARGETS)
    with open(result_files[0], 'w') as output:
        context.run_external(f'trf {fasta_file} 2 7 7 80 10 50 500 -d -h -ngs', stdout=output)

def _read_results(genome, contig_dict, db_connection, result_files) -> int:
    count = 0
    with open(result_files[0]) as trf_handle:
        for line in trf_handle:
            if line.startswith("@"):
                contig_name = line[1:].strip()
                if space_index := contig_name.find(' ') > 0:
                    contig_name = contig_name[:space_index]
                contig = contig_dict[contig_name]
                continue
            if not contig:
                continue
            words = line.split()
            start = int(words[0]) - 1
            end = int(words[1])
            seq = contig['seq'][start:end]
            feature = sqlite.Feature(genome = genome.name,
                                     contig = contig['id'],
                                     start = start,
                                     end = end,
                                     strand = 1,
                                     type = 'repeat',
                                     inference = 'tandem-repeat-finder',
                                     nt_seq = seq,
                                     notes = f'period size {words[2]}; copies {words[3]}')
            sqlite.add_new_feature_to_db(db_connection, feature)
            count += 1
    return count

@context.register_annotator
def run_and_read_trf():
    return ({'pipeline_position': 41,
             'annotator_key': 'trf',
             'purpose': 'tandem repeat prediction with trf',
             'programs': ('trf',),
             'result_files': ('tandem-repeat-finder',),
             'run': _run_programs,
             'read': _read_results})
