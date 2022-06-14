from pathlib import Path
import shutil
import os
from concurrent.futures import ProcessPoolExecutor

from metaerg.data_model import MetaergSeqFeature, MetaergGenome, FeatureType, DBentry
from metaerg import context
from metaerg import bioparsers

def _run_programs(genome:MetaergGenome, result_files):
    cds_aa_file = context.spawn_file('cds.faa', genome.id)
    cdd_database = Path(context.DATABASE_DIR, 'cdd', 'Cdd')
    if context.CPUS_PER_GENOME > 1:
        split_fasta_files = bioparsers.write_genome_fasta_files(genome, cds_aa_file, context.CPUS_PER_GENOME,
                                                                target=FeatureType.CDS)
        split_cdd_files = [Path(result_files[0].parent, f'{result_files[0].name}.{i}')
                           for i in range(len(split_fasta_files))]
        with ProcessPoolExecutor(max_workers=context.CPUS_PER_GENOME) as executor:
            for split_input, split_output in zip(split_fasta_files, split_cdd_files):
                executor.submit(context.run_external, f'rpsblast -db {cdd_database} -query {split_input} '
                                                      f'-out {split_output} -outfmt 6 -evalue 1e-7')

        with open(result_files[0], 'wb') as output:
            for split_input_file, split_output_file in zip(split_fasta_files, split_cdd_files):
                with open(split_output_file, 'rb') as input:
                    shutil.copyfileobj(input, output)
                split_input_file.unlink()
                split_output_file.unlink()
    else:
        context.run_external(f'rpsblast -db {cdd_database} -query {cds_aa_file} -out {result_files[0]} -outfmt 6 -evalue 1e-7')



def _read_results(genome:MetaergGenome, result_files) -> int:
    # load cdd index
    cdd = {}
    cdd_index = Path(context.DATABASE_DIR, 'cdd', 'cddid.tbl')
    with open(cdd_index) as db_handle:
        for line in db_handle:
            words = line.split("\t")
            cdd[int(words[0])] = DBentry(words[0], words[1], words[2], '', int(words[4]), 0)
    context.log(f'Parsed {len(cdd)} entries from conserved domain database.')
    # parse cdd results
    cdd_result_count = 0
    with bioparsers.TabularBlastParser(result_files, 'BLAST', lambda i: cdd[int(i[4:])]) as handle:
        for blast_result in handle:
            feature: MetaergSeqFeature = genome.get_feature(blast_result.query)
            cdd_result_count += 1
            feature.cdd = blast_result
            genome.subsystems.match(feature, (h.hit.descr for h in blast_result.hits
                                              if h.aligned_length / h.hit.length >= 0.8))
            top_entry = blast_result.hits[0].hit
            feature.descr = f'{top_entry.id}|{top_entry.gene} {top_entry.descr}'
            if len(feature.descr) > 35:
                    feature.descr = feature.descr[:35] + '...'
    return cdd_result_count


@context.register_annotator
def run_and_read_cdd():
    return ({'pipeline_position': 71,
             'purpose': 'function prediction using RPSBlast and the conserved domain database',
             'programs': ('rpsblast',),
             'databases': (Path(context.DATABASE_DIR, 'cdd', 'cddid.tbl'), Path(context.DATABASE_DIR, 'cdd', 'Cdd')),
             'result_files': ('cdd',),
             'run': _run_programs,
             'read': _read_results})


@context.register_database_installer
def install_cdd_database():
    if 'C' not in context.CREATE_DB_TASKS:
       return
    cdd_dir = Path(context.DATABASE_DIR, 'cdd')
    if context.FORCE or not cdd_dir.exists():
        context.log(f'Installing the conserved domain database to {cdd_dir}...')
        cdd_dir.mkdir(exist_ok=True, parents=True)
        context.run_external(f'wget -P {cdd_dir} https://ftp.ncbi.nih.gov/pub/mmdb/cdd/cddid.tbl.gz')
        cdd_index = Path(cdd_dir, 'cddid.tbl')
        context.run_external(f'gunzip {cdd_index}.gz')

        temp_cdd_dir = Path(context.DATABASE_DIR, 'cdd-temp')
        temp_cdd_dir.mkdir(exist_ok=True, parents=True)
        context.run_external(f'wget -P {temp_cdd_dir} https://ftp.ncbi.nih.gov/pub/mmdb/cdd/cdd.tar.gz')
        context.run_external(f'tar -xf {Path(cdd_dir, "cdd.tar.gz")} -C {temp_cdd_dir}')
        current_dir = os.getcwd()
        os.chdir(temp_cdd_dir)
        context.run_external(f'makeprofiledb -title CDD.v.3.12 -in {Path(temp_cdd_dir, "Cdd.pn")} -out '
                             f'{Path(cdd_dir, "Cdd")} -threshold 9.82 -scale 100.0 -dbtype rps '
                             f'-index true')
        os.chdir(current_dir)
    else:
       context.log(f'Keeping existing conserved domain database in {cdd_dir}, use --force to overwrite.')
