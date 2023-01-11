import shutil
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor

from metaerg import context
from metaerg.datatypes import fasta
from metaerg.datatypes import sqlite


def _run_programs(genome_name, contig_dict, db_connection, result_files):
    cds_aa_file = context.spawn_file('cds.faa', genome_name)
    if context.CPUS_PER_GENOME > 1:
        split_fasta_files = fasta.write_features_to_fasta(db_connection, 'aa', cds_aa_file, context.CPUS_PER_GENOME,
                                                          targets=('CDS',))
        split_signalp_files = [result_files[0].parent / f'{result_files[0].name}.{i}'
                               for i in range(len(split_fasta_files))]
        with ProcessPoolExecutor(max_workers=context.CPUS_PER_GENOME) as executor:
            for split_input, split_output in zip(split_fasta_files, split_signalp_files):
                executor.submit(context.run_external, f'signalp6 --fastafile {split_input} --output_dir '
                                                    f'{split_output} --format none --organism other')

        result_files[0].mkdir(exist_ok=True)
        with open(result_files[0] / 'prediction_results.txt', 'wb') as output:
            for split_cds_aa_file, split_signalp_dir in zip(split_fasta_files, split_signalp_files):
                signalp_result_file = Path(split_signalp_dir, 'prediction_results.txt')
                if signalp_result_file.exists():
                    with open(signalp_result_file, 'rb') as input:
                        shutil.copyfileobj(input, output)
                else:
                    context.log(f'({genome_name}) WARNING - missing part of signalp output!')
                if split_signalp_dir.exists():
                    shutil.rmtree(split_signalp_dir)
                split_cds_aa_file.unlink(missing_ok=True)
    else:
        context.run_external(f'signalp6 --fastafile {cds_aa_file} --output_dir {result_files[0]} --format none --organism other')


def _read_results(genome_name, contig_dict, db_connection, result_files) -> int:
    count = 0
    signalp_result_file = result_files[0] / 'prediction_results.txt'
    if signalp_result_file.exists():
        with open(signalp_result_file) as signalp_handle:
            for line in signalp_handle:
                if line.startswith("#"):
                    continue
                words = line.split("\t")
                if "OTHER" == words[1]:
                    continue
                feature_id = words[0].split()[0]
                feature = sqlite.read_feature_by_id(db_connection, feature_id)
                if not feature:
                    raise Exception(f'Found results for unknown feature {feature_id}, '
                                    f'may need to rerun metaerg with --force')
                feature.signal_peptide = words[1]
                sqlite.update_feature_in_db(db_connection, feature)
                count += 1
    return count


@context.register_annotator
def run_and_read_signalp():
    return ({'pipeline_position': 121,
             'annotator_key': 'signalp',
             'purpose': 'signal peptide prediction with signalp',
             'programs': ('signalp6',),
             'result_files': ('signalp',),
             'run': _run_programs,
             'read': _read_results})
