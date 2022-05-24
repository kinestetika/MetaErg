from pathlib import Path
import shutil
from concurrent.futures import ProcessPoolExecutor

from metaerg.run_and_read.data_model import MetaergSeqFeature, MetaergGenome, FeatureType, TabularBlastParser, DBentry
from metaerg.run_and_read.context import register_annotator, spawn_file, run_external, DATABASE_DIR, CPUS_PER_GENOME, log

def _run_programs(genome:MetaergGenome, result_files):
    cds_aa_file = spawn_file('cds.faa', genome.id)
    cdd_database = Path(DATABASE_DIR, "Cdd")
    if CPUS_PER_GENOME > 1:
        split_fasta_files = genome.write_fasta_files(cds_aa_file, CPUS_PER_GENOME, target=FeatureType.CDS)
        split_cdd_files = [Path(result_files[0].parent, f'{result_files[0].name}.{i}')
                           for i in range(len(split_fasta_files))]
        with ProcessPoolExecutor(max_workers=CPUS_PER_GENOME) as executor:
            for split_input, split_output in zip(split_fasta_files, split_cdd_files):
                executor.submit(run_external, f'rpsblast -db {cdd_database} -query {split_input} '
                                              f'-out {split_output} -outfmt 6 -evalue 1e-7')

        with open(result_files[0], 'wb') as output:
            for split_input_file, split_output_file in zip(split_fasta_files, split_cdd_files):
                with open(split_output_file, 'rb') as input:
                    shutil.copyfileobj(input, output)
                split_input_file.unlink()
                split_output_file.unlink()
    else:
        run_external(f'rpsblast -db {cdd_database} -query {cds_aa_file} -out {result_files[0]} -outfmt 6 -evalue 1e-7')



def _read_results(genome:MetaergGenome, result_files) -> int:
    # load cdd index
    cdd = {}
    cdd_index = Path(DATABASE_DIR, 'cddid.tbl')
    with open(cdd_index) as db_handle:
        for line in db_handle:
            words = line.split("\t")
            cdd[int(words[0])] = DBentry(words[0], words[1], words[2], '', int(words[4]), 0)
    log(f'Parsed {len(cdd)} entries from conserved domain database.')
    # parse cdd results
    cdd_result_count = 0
    with TabularBlastParser(result_files, 'BLAST', lambda i: cdd[int(i[4:])]) as handle:
        for blast_result in handle:
            feature: MetaergSeqFeature = genome.get_feature(blast_result.query)
            cdd_result_count += 1
            feature.cdd = blast_result
            genome.subsystems.match(feature, (h.hit.descr for h in blast_result.hits
                                              if h.aligned_length / h.hit.length >= 0.8))
            top_entry = blast_result.hits[0].hit
            feature.product = f'{top_entry.id}|{top_entry.gene} {top_entry.descr}'
            if len(feature.product) > 35:
                    feature.product = feature.product[:35] + '...'
    return cdd_result_count


@register_annotator
def run_and_read_cdd():
    return ({'pipeline_position': 71,
             'purpose': 'function prediction using RPSBlast and the conserved domain database',
             'programs': ('rpsblast',),
             'databses': ('cddid.tbl', 'Cdd'),
             'result_files': ('cdd',),
             'run': _run_programs,
             'read': _read_results})
