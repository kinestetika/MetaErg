import shutil
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor

from metaerg.run_and_read.data_model import MetaergSeqFeature, FeatureType, MetaergGenome
from metaerg.run_and_read.context import register, spawn_file, run_external, CPUS_PER_GENOME, log


def _run_programs(genome:MetaergGenome, result_files):
    cds_aa_file = spawn_file('cds.faa', genome.id)
    if CPUS_PER_GENOME > 1:
        split_fasta_files = genome.write_fasta_files(cds_aa_file, CPUS_PER_GENOME, target=FeatureType.CDS)
        split_signalp_files = [Path(result_files[0].parent, f'{result_files[0].name}.{i}')
                               for i in range(len(split_fasta_files))]
        with ProcessPoolExecutor(max_workers=CPUS_PER_GENOME) as executor:
            for split_input, split_output in zip(split_fasta_files, split_signalp_files):
                executor.submit(run_external, f'signalp6 --fastafile {split_input} --output_dir '
                                                    f'{split_output} --format none --organism other')

        result_files[0].mkdir()
        with open(Path(result_files[0], 'prediction_results.txt'), 'wb') as output:
            for split_cds_aa_file, split_signalp_dir in zip(split_fasta_files, split_signalp_files):
                signalp_result_file = Path(split_signalp_dir, 'prediction_results.txt')
                if signalp_result_file.exists():
                    with open(signalp_result_file, 'rb') as input:
                        shutil.copyfileobj(input, output)
                else:
                    log(f'({genome.id}) WARNING - missing part of signalp output!')
                shutil.rmtree(split_signalp_dir)
                split_cds_aa_file.unlink()
    else:
        run_external(f'signalp6 --fastafile {cds_aa_file} --output_dir {result_files[0]} --format none --organism other')


def _read_results(genome:MetaergGenome, result_files) -> int:
    count = 0
    with open(Path(result_files[0], 'prediction_results.txt')) as signalp_handle:
        for line in signalp_handle:
            if line.startswith("#"):
                continue
            words = line.split("\t")
            if "OTHER" == words[1]:
                continue
            feature: MetaergSeqFeature = genome.get_feature(words[0].split()[0])
            feature.signal_peptide = words[1]
            count += 1
    return count


@register
def run_and_read_signalp():
    return ({'pipeline_position': 121,
             'purpose': 'signal peptide prediction with signalp',
             'programs': ('signalp6'),
             'result_files': ('signalp',),
             'run': _run_programs,
             'read': _read_results})
