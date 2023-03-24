import argparse
import tarfile
from pathlib import Path
from concurrent import futures
from hashlib import md5

import pandas as pd

import metaerg.run_and_read.diamond_and_blastn
import metaerg.run_and_read.functional_genes
import metaerg.run_and_read.antismash
from metaerg.datatypes import sqlite
from metaerg.run_and_read import *
from metaerg.html import *
from metaerg import context
from metaerg import registry
from metaerg.datatypes import functional_genes
from metaerg.datatypes import fasta
from metaerg.datatypes import gbk
from metaerg.datatypes.excel import write_genomes_to_xls
from metaerg.html import html_all_genomes
from metaerg.run_and_read import tmhmm
from metaerg.installation import install_all_helper_programs

VERSION = "2.3.38"


def parse_arguments():
    parser = argparse.ArgumentParser(description='metaerg.py. (C) Marc Strous, Xiaoli Dong 2019, 2022')
    parser.add_argument('--contig_file', help='Fasta nucleotide file of the contigs, or dir that '
                                              'contains multiple fasta nucleotide files.')
    parser.add_argument('--database_dir', help='Dir that contains the annotation databases.')
    parser.add_argument('--output_dir', default='', help='Dir that will contain the annotation results.')
    parser.add_argument('--contig_mode', default=False,  action="store_true",
                        help='Annotate contigs individually instead of assuming they are part of a genome, MAG or bin.')
    parser.add_argument('--rename_contigs', default=False,  action="store_true",
                        help='Renaming contigs can improve visualization and presentation of results.')
    parser.add_argument('--rename_genomes', default=False,  action="store_true",
                        help='Renaming genomes can improve visualization and presentation of results.')
    parser.add_argument('--min_contig_length', default=0,  help='Shorter contigs will be filtered before annotaton.')
    parser.add_argument('--cpus', default=0, help='How many cpus/threads to use (default: all = 0).')
    parser.add_argument('--file_extension', default='.fna', help='When annotating multiple files in a folder, extension'
                                                                 'of the fasta nucleotide files (default: .fna).')
    parser.add_argument('--translation_table', default='', help='Which translation tables to use (default 11,25).')
    parser.add_argument('--delimiter', default='.', help='Separater character used in feature ids.')
    parser.add_argument('--prefix', default='g', help='Prefix used when renaming genomes (default: g).')
    parser.add_argument('--checkm_dir', default='checkm', help='Dir with the checkm results (default: checkm)')
    parser.add_argument('--gtdbtk_dir', default='gtdbtk', help='Dir with the gtdbtk results (default: gtdbtk).')
    parser.add_argument('--download_database', default=False, action="store_true",
                        help='Download ready-made metaerg database.')
    parser.add_argument('--create_database', default='', help='Create metaerg database from scratch (use "all", '
                                                                 'to create all components of the database.). Use '
                                                                 'any combination of PVEBRCSA to only create specific '
                                                                 'parts of the database (see README)')
    parser.add_argument('--install_deps', default='', help='Dir for installation of all dependencies '
                                                           '(helper programs). Dependencies will be installed here.')
    parser.add_argument('--path_to_signalp', default='', help='Path to signalp-6.0g.fast.tar.gz.')
    parser.add_argument('--path_to_tmhmm', default='', help='Path to tmhmm-2.0c.Linux.tar.gz.')
    parser.add_argument('--force', default='',  help='Use force to overwrite previous result files. Use "--force all" to redo '
                                                     'everything, or antismash, aragorn, cdd, cmscan, diamond_and_blastn, hmm, '
                                                     'ltr_harvest, minced, prodigal, signalp, repeat_masker, tmhmm, trf, '
                                                     'separated by commas (,) to redo specific steps')
    parser.add_argument('--update_annotations', default=False, action='store_true', help="Do not run any helper programs, only "
                                                                                "update annotations woth results from previous runs.")
    parser.add_argument('--skip_step', default='', help="Skip one or more annotation steps. Steps are: antismash, aragorn, "
                                                        "cdd, cmscan, diamond_and_blastn, hmm, ltr_harvest, minced, prodigal, "
                                                        "signalp, repeat_masker, tmhmm, trf, separated by commas (,)")

    return parser.parse_args()


def annotate_genome(genome_name, input_fasta_file: Path):
    context.log(f'Started annotation of {genome_name}...')
    current_progress = context.parse_metaerg_progress(genome_name)
    if context.ANNOTATOR_STATUS['visualization'] != context.FORCE_ANNOTATOR and \
            current_progress.get('visualization', 0) == context.PROGRESS_COMPLETE:
        context.log(f'({genome_name}) already completed!')
        return
    # (1) prepare feature sqlite3 database
    #feature_db_file = context.spawn_file('annotations.sqlite', genome_name, context.BASE_DIR)
    #sqlite.create_db(feature_db_file)
    feature_db_connection = sqlite.create_db(target='Features')
    # (2) load sequence data
    contig_dict = fasta.load_contigs(genome_name, input_fasta_file, delimiter=context.DELIMITER,
                                     min_contig_length=context.MIN_CONTIG_LENGTH, rename_contigs=context.RENAME_CONTIGS)
    if not len(contig_dict):
        context.log(f'({genome_name}) WARNING: {input_fasta_file} appears to contain no fasta data... aborting!')
        return
    genome = sqlite.Genome(name=genome_name, input_fasta_file=input_fasta_file)

    contigs: list[dict] = list(contig_dict.values())
    contigs.sort(key=lambda c: len(c['seq']), reverse=True)

    genome.size = sum(len(c['seq']) for c in contigs)
    genome.number_of_contigs = len(contigs)
    genome.fraction_gc = sum((c['seq'].count('G') + c['seq'].count('G') for c in contigs)) / \
                       (genome.size - sum((c['seq'].count('N') for c in contigs)))
    cum_size = 0
    for c in contigs:
        cum_size += len(c['seq'])
        if cum_size > + genome.size / 2:
            genome.n50_contig_length = len(c['seq'])
            break

    # (3) now annotate
    for annotator in context.sorted_annotators():
        annotator(genome, contig_dict, feature_db_connection)
    # feather_file = context.spawn_file("all_genes.feather", genome_name, context.BASE_DIR)
    # feature_data = pd.read_feather(feather_file)
    # feature_data = feature_data.set_index('id', drop=False)

    # (4) save results
    context.log(f'({genome_name}) Now writing annotations to .fasta, .gbk...')
    faa_file = context.spawn_file("faa", genome_name, context.BASE_DIR, extension='faa')
    rna_file = context.spawn_file("rna.fna", genome_name, context.BASE_DIR, extension='fna')
    gbk_file = context.spawn_file("gbk", genome_name, context.BASE_DIR, extension='gbk')
    fna_file = context.spawn_file("fna", genome_name, context.BASE_DIR, extension='fna')
    fasta.write_features_to_fasta(feature_db_connection, 'aa', faa_file, targets=('CDS',))
    fasta.write_features_to_fasta(feature_db_connection, 'nt', rna_file, targets=('rRNA tRNA tmRNA ncRNA retrotransposon'.split()))
    fasta.write_contigs_to_fasta(contig_dict, fna_file, feature_db_connection)
    with open(gbk_file, 'w') as gbk_writer:
        gbk.gbk_write_genome(gbk_writer, contig_dict, feature_db_connection)
    feature_count = sqlite.count_features(feature_db_connection)
    if feature_count < 1e6:
        context.log(f'({genome_name}) Writing {feature_count} annotations to .feather format...')
        feather_file = context.spawn_file("annotations.feather", genome_name, context.BASE_DIR, extension='pyarrow')
        rows = []
        for feature in sqlite.read_all_features(feature_db_connection):
            rows.append({k: str(v) for k, v in feature})
        pd.DataFrame(rows).to_feather(feather_file)
    else:
        context.log(f'({genome_name}) Skipping feather format, too many ({feature_count}) annotations...')
    # (5) visualize
    context.log(f'({genome_name}) Now writing final result as .html for visualization...')
    for html_writer in registry.HTML_WRITER_REGISTRY:
        html_writer(genome, feature_db_connection, context.HTML_DIR)
    # (6) update progress
    current_progress = context.parse_metaerg_progress(genome_name)
    current_progress['visualization'] = context.PROGRESS_COMPLETE
    context.write_metaerg_progress(genome_name, current_progress)
    context.log(f'({genome_name}) Completed html visualization.')
    # (7) save feature sqlite database
    feature_db_file = context.spawn_file('annotations.sqlite', genome_name, context.BASE_DIR)
    feature_db_file.unlink(missing_ok=True)
    sqlite.write_db(feature_db_connection, feature_db_file)
    feature_db_connection.close()
    return genome


def main():
    print(f'This is metaerg.py {VERSION}')
    context.init(**parse_arguments().__dict__)
    if context.METAERG_MODE == context.METAERG_MODE_CREATE_DATABASE:
        context.log(f'Creating/installing/downloading metaerg databases. Tasks: {context.DATABASE_TASKS}; ')
        for db_installer in registry.DATABASE_INSTALLER_REGISTRY:
            db_installer()
    elif context.METAERG_MODE == context.METAERG_MODE_DOWNLOAD_DATABASE:
        context.log('Downloading premade databases from https://object-arbutus.cloud.computecanada.ca...')
        database_tarbal_file = context.DATABASE_DIR / 'metaerg_db_207_v2.tar.gz'
        context.download('https://object-arbutus.cloud.computecanada.ca/metaerg/metaerg_2.25_gtdb_207_v2.tar.gz',
                         database_tarbal_file)
        md5sum = md5(open(database_tarbal_file,'rb').read()).hexdigest()
        print(md5sum)
        if '48ffe7150711dd6f982e1f3d759e4ff9' == md5sum:
            context.log(f'checksum {md5sum} as expected.')
        else:
            raise Exception('Downloaded database has incorrect checksum - download failed. Aborting...')
        context.log('Now extracting databases from tar archive...')
        database_archive = tarfile.open(database_tarbal_file)
        database_archive.extractall(context.DATABASE_DIR)
        # database_tarbal_file.unlink()
        context.DATABASE_TASKS = 'SA'
        metaerg.run_and_read.diamond_and_blastn.compile_databases()
        metaerg.run_and_read.functional_genes.install_functional_gene_databases()
        metaerg.run_and_read.antismash.format_antismash_databases()
    elif context.METAERG_MODE == context.METAERG_MODE_INSTALL_DEPS:
        install_all_helper_programs(context.BIN_DIR, context.PATH_TO_SIGNALP, context.PATH_TO_TMHMM)
    else:
        # sqlite.create_db(genome_db_file, target='Genomes')
        genome_db_connection = sqlite.create_db(target='Genomes')
        for preloader in registry.DB_PRELOADER_REGISTRY:
            preloader()
        if context.PARALLEL_ANNOTATIONS > 1:
            outcomes = {}
            with futures.ProcessPoolExecutor(max_workers=context.PARALLEL_ANNOTATIONS) as executor:
                for genome_name, contig_file in zip(context.GENOME_NAMES, context.CONTIG_FILES):
                    outcomes[genome_name] = executor.submit(annotate_genome, genome_name, contig_file)
            for genome_name, future in outcomes.items():
                try:
                    genome = genome=future.result()
                    if genome:
                        sqlite.add_new_genome_to_db(genome_db_connection, genome)
                except Exception:
                    context.log(f'({genome_name}) Error while processing:')
                    raise
        else:
            for genome_name, contig_file in zip(context.GENOME_NAMES, context.CONTIG_FILES):
                sqlite.add_new_genome_to_db(sql_connection=genome_db_connection,
                                            genome=annotate_genome(genome_name, contig_file))
        tmhmm.cleanup(context.TEMP_DIR)

        context.log('Now writing all-genomes overview to html...')
        html_all_genomes.write_html(genome_db_connection, context.HTML_DIR)

        context.log('Now writing all-genomes overview to sqlite...')
        genome_db_file = Path(context.BASE_DIR, 'genome_properties.sqlite')
        genome_db_file.unlink(missing_ok=True)
        sqlite.write_db(genome_db_connection, genome_db_file)
        context.log('Now writing all-genomes overview to excel...')
        write_genomes_to_xls(genome_db_connection)
        genome_db_connection.close()
        context.log(f'Done. Thank you for using metaerg.py {VERSION}')


if __name__ == "__main__":
    main()
