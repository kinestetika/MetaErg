import re
import shutil
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor

from metaerg import context
from metaerg.datatypes import fasta
from metaerg.datatypes import sqlite


def _run_programs(genome, contig_dict, db_connection, result_files):
    cds_aa_file = context.spawn_file('cds.faa', genome.name)
    context.run_external(f'signalp6 --fastafile {cds_aa_file} --output_dir {result_files[0]} --format none '
                         f'--organism other -wp 1')

def _read_results(genome, contig_dict, db_connection, result_files) -> int:
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
                    raise Exception(f'Found signalp result for unknown feature {feature_id}, '
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
