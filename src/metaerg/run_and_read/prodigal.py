import re

import datatypes.sqlite
from metaerg import context
from metaerg.datatypes import fasta
from metaerg.datatypes import sqlite


def _run_programs(genome_name, contig_dict, db_connection, result_files):
    fasta_file = context.spawn_file('masked', genome_name)
    # no masking here becasuse we want to arbitrate with repeatmasker results
    #fasta.write_contigs_to_fasta(contig_dict, fasta_file, feature_data, genome_name,
    #                             mask_targets=fasta.ALL_MASK_TARGETS)
    context.run_external(f'prodigal -g {context.TRANSLATION_TABLE} -m -f gff -q -i {fasta_file} -a {result_files[0]}'
                         f' -d {result_files[1]}')


def _read_results(genome_name, contig_dict, db_connection, result_files) -> int:
    ORF_ID_PATTERN = re.compile(r'_(\d+?)$')
    nucl_seq_hash = {}
    with fasta.FastaParser(result_files[1], cleanup_seq=False) as fasta_reader:
        for seq_rec in fasta_reader:
            nucl_seq_hash[seq_rec['id']] = seq_rec['seq']
    count = 0
    feature_count_before_arbritration = sum(1 for f in sqlite.read_all_features(db_connection))
    with fasta.FastaParser(result_files[0], cleanup_seq=False) as fasta_reader:
        rejected_cds_count = 0
        for seq_rec in fasta_reader:
            x_count = seq_rec['seq'].count('X')
            if x_count > 0.2 * len(seq_rec['seq'] or x_count >= 10):
                rejected_cds_count += 1
                continue
            words = seq_rec['descr'].split('#')
            try:
                m = ORF_ID_PATTERN.search(seq_rec['id'])
                contig_id = seq_rec['id'][0:m.start()]
                contig = contig_dict[contig_id]
            except KeyError:
                context.log(f'({genome_name}) Warning: Failed to find contig with "{seq_rec["id"]}"')
                continue
            start = int(words[1].strip()) - 1
            end = int(words[2].strip())

            # reconciliation with repeats
            overlapping_features = [f for f in sqlite.read_all_features(db_connection, contig=contig_id, location=(start, end))]

            if len(overlapping_features):
                overlap = 0
                for f in overlapping_features:
                    overlap += min(f.end, end) - max(start, f.start)
                non_repeat_features = [f for f in overlapping_features if f.type in datatypes.sqlite.RNA_TARGETS]
                if len(non_repeat_features):
                    rejected_cds_count += 1
                    continue
                elif overlap < 0.33 * (end-start):
                    for overlapping_feature in overlapping_features:
                        sqlite.drop_feature(db_connection, overlapping_feature)
                else:
                    rejected_cds_count += 1
                    continue

            strand = int(words[3].strip())
            if seq_rec['seq'].endswith('*'):
                seq_rec['seq'] = seq_rec['seq'][:-1]
            feature = sqlite.Feature(genome = genome_name,
                       contig = contig_id,
                       start = start,
                       end = end,
                       strand = strand,
                       type = 'CDS',
                       inference = 'prodigal',
                       aa_seq = seq_rec['seq'],
                       nt_seq = nucl_seq_hash[seq_rec['id']])
            if 'partial=01' in seq_rec['descr'] or 'partial=01' in seq_rec['descr'] or 'partial=11' in seq_rec['descr']:
                feature.notes = 'partial protein'
            sqlite.add_new_feature_to_db(db_connection, feature)

        feature_count_after_arbritration = sum(1 for f in sqlite.read_all_features(db_connection))
        context.log(f'({genome_name}) Dropped {feature_count_before_arbritration - feature_count_after_arbritration} repeats and'
                    f' rejected {rejected_cds_count} CDS during arbitration of prodigal results.')
        return count


@context.register_annotator
def run_and_read_prodigal():
    return ({'pipeline_position': 61,
             'annotator_key': 'prodigal',
             'purpose': 'coding sequence prediction with prodigal',
             'programs': ('prodigal',),
             'result_files': ('prodigal','prodigal-nucl'),
             'run': _run_programs,
             'read': _read_results})
