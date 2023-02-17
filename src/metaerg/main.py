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
from metaerg.datatypes.genome_properties import compute_genome_properties, write_genome_properties_to_xls
from metaerg.html import html_all_genomes
from metaerg.run_and_read import tmhmm
from metaerg.installation import install_all_helper_programs

VERSION = "2.3.21"


def parse_arguments():
    parser = argparse.ArgumentParser(description='metaerg.py. (C) Marc Strous, Xiaoli Dong 2019, 2022')
    parser.add_argument('--contig_file', help='Fasta nucleotide file of the contigs, or dir that '
                                              'contains multiple fasta nucleotide files.')
    parser.add_argument('--database_dir', help='Dir that contains the annotation databases.')
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
    parser.add_argument('--translation_table', default=0, help='Which translation table to use (default 11).')
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
    parser.add_argument('--force', default=False, action="store_true",
                        help='Use force to overwrite previous result files.')
    parser.add_argument('--read_only', default=False, action='store_true', help="Do not run any helper programs, only "
                                                                                "read results from previous runs.")
    parser.add_argument('--skip_step', default='', help="Skip one or more annotation steps. Steps are: antismash, aragorn, "
                                                        "cdd, cmscan, diamond_and_blastn, hmm, ltr_harvest, minced, prodigal, "
                                                        "signalp, repeat_masker, tmhmm, trf, separated by commas (,)")
    parser.add_argument('--run_step', default='', help="Only (re)run one or more annotation steps. Steps are: antismash, aragorn, "
                                                        "cdd, cmscan, diamond_and_blastn, hmm, ltr_harvest, minced, prodigal, "
                                                        "signalp, repeat_masker, tmhmm, trf, separated by commas (,)")

    return parser.parse_args()


def annotate_genome(genome_name, input_fasta_file: Path):
    context.log(f'Started annotation of {genome_name}...')
    current_progress = context.parse_metaerg_progress(genome_name)
    if not context.FORCE and 'visualization=complete' in current_progress:
        context.log(f'({genome_name}) already completed!')
        return
    # (1) prepare sqlite3 database
    db_file = context.spawn_file('annotations.sqlite', genome_name, context.BASE_DIR)
    sqlite.create_db(db_file)
    db_connection = sqlite.connect_to_db(db_file)
    # (2) load sequence data
    contig_dict = fasta.load_contigs(genome_name, input_fasta_file, delimiter=context.DELIMITER,
                                     min_contig_length=context.MIN_CONTIG_LENGTH, rename_contigs=context.RENAME_CONTIGS)
    # (3) now annotate
    for annotator in context.sorted_annotators():
        annotator(genome_name, contig_dict, db_connection)
    # feather_file = context.spawn_file("all_genes.feather", genome_name, context.BASE_DIR)
    # feature_data = pd.read_feather(feather_file)
    # feature_data = feature_data.set_index('id', drop=False)

    # (4) save results
    context.log(f'({genome_name}) Now writing annotations to .fasta, .gbk...')
    faa_file = context.spawn_file("faa", genome_name, context.BASE_DIR)
    rna_file = context.spawn_file("rna.fna", genome_name, context.BASE_DIR)
    gbk_file = context.spawn_file("gbk", genome_name, context.BASE_DIR)
    fna_file = context.spawn_file("fna", genome_name, context.BASE_DIR)
    fasta.write_features_to_fasta(db_connection, 'aa', faa_file, targets=('CDS',))
    fasta.write_features_to_fasta(db_connection, 'nt', rna_file, targets=('rRNA tRNA tmRNA ncRNA retrotransposon'.split()))
    fasta.write_contigs_to_fasta(contig_dict, fna_file, db_connection)
    with open(gbk_file, 'w') as gbk_writer:
        gbk.gbk_write_genome(gbk_writer, contig_dict, db_connection)
    feature_count = sqlite.count_features(db_connection)
    if feature_count < 1e6:
        context.log(f'({genome_name}) Writing {feature_count} annotations to .feather format...')
        feather_file = context.spawn_file("annotations.feather", genome_name, context.BASE_DIR)
        rows = []
        for feature in sqlite.read_all_features(db_connection):
            rows.append({k: str(v) for k, v in feature})
        pd.DataFrame(rows).to_feather(feather_file)
    else:
        context.log(f'({genome_name}) Skipping feather format, too many ({feature_count}) annotations...')
    # (5) visualize
    genome_properties = compute_genome_properties(genome_name, input_fasta_file, contig_dict, db_connection)
    context.log(f'({genome_name}) Now writing final result as .html for visualization...')
    for html_writer in registry.HTML_WRITER_REGISTRY:
        html_writer(genome_name, db_connection, genome_properties, context.HTML_DIR)
    # (5) update progress
    current_progress = context.parse_metaerg_progress(genome_name)
    current_progress += 'visualization=complete\n'
    context.write_metaerg_progress(genome_name, current_progress)
    context.log(f'({genome_name}) Completed html visualization.')
    return genome_properties


def main():
    print(f'This is metaerg.py {VERSION}')
    context.init(**parse_arguments().__dict__)
    if context.METAERG_MODE == context.METAERG_MODE_CREATE_DATABASE:
        context.log(f'Creating/installing/downloading metaerg databases. Tasks: {context.TASKS}; '
                    f'force: {context.FORCE}.')
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
        context.TASKS = 'SA'
        metaerg.run_and_read.diamond_and_blastn.compile_databases()
        metaerg.run_and_read.functional_genes.install_functional_gene_databases()
        metaerg.run_and_read.antismash.format_antismash_databases()
    elif context.METAERG_MODE == context.METAERG_MODE_INSTALL_DEPS:
        install_all_helper_programs(context.BIN_DIR, context.PATH_TO_SIGNALP, context.PATH_TO_TMHMM)
    else:
        functional_genes.init_functional_gene_config()
        annotation_summaries = {}
        if context.PARALLEL_ANNOTATIONS > 1:
            outcomes = {}
            with futures.ProcessPoolExecutor(max_workers=context.PARALLEL_ANNOTATIONS) as executor:
                for genome_name, contig_file in zip(context.GENOME_NAMES, context.CONTIG_FILES):
                    outcomes[genome_name] = executor.submit(annotate_genome, genome_name, contig_file)
            for genome_name, future in outcomes.items():
                try:
                    annotation_summaries[genome_name] = future.result()
                except Exception:
                    context.log(f'({genome_name}) Error while processing:')
                    raise
        else:
            for genome_name, contig_file in zip(context.GENOME_NAMES, context.CONTIG_FILES):
                annotation_summaries[genome_name] = annotate_genome(genome_name, contig_file)
        tmhmm.cleanup(context.BASE_DIR)
        context.log('Now writing all-genomes overview html...')
        html_all_genomes.write_html(context.HTML_DIR)
        context.log('Now writing excel files with functional genes per genome...')
        write_genome_properties_to_xls(annotation_summaries)
        context.log(f'Done. Thank you for using metaerg.py {VERSION}')


if __name__ == "__main__":
    main()
