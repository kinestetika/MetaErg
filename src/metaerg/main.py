import argparse
import tarfile
from pathlib import Path
from concurrent import futures
from hashlib import md5
from collections import Counter

import pandas as pd

import metaerg.run_and_read.diamond_and_blastn
import metaerg.run_and_read.functional_genes
import metaerg.run_and_read.antismash
from metaerg.datatypes import sqlite
from metaerg import context
from metaerg import registry
from metaerg import functional_gene_configuration
from metaerg.datatypes import fasta, gbk
from metaerg.html import html_all_genomes
from metaerg.run_and_read import tmhmm
from metaerg.installation import install_all_helper_programs
from metaerg.calculations.codon_usage_bias import compute_codon_bias_estimate_doubling_time
from metaerg.run_and_read import *
from metaerg.html import *

VERSION = "2.3.10"


def parse_arguments():
    parser = argparse.ArgumentParser(description='metaerg.py. (C) Marc Strous, Xiaoli Dong 2019, 2022')
    parser.add_argument('--contig_file', help='Fasta nucleotide file of the contigs, or dir that '
                                              'contains multiple fasta nucleotide files.')
    parser.add_argument('--database_dir', help='Dir that contains the annotation databases.')
    parser.add_argument('--contig_mode', default=False,  action="store_true",
                        help='Annotate contigs individually instead of assuming they are part of a genome, MAG or bin.')
    parser.add_argument('--skip', default='', help="Skip one or more annotation steps. Steps are: antismash, aragorn, "
                                                   "cdd, cmscan, diamond_and_blastn, hmm, ltr_harvest, minced, prodigal, "
                                                   "signalp, repeat_masker, tmhmm, trf, separated by commas (,)")
    parser.add_argument('--rename_contigs', default=False,  action="store_true",
                        help='Renaming contigs can improve visualization and presentation of results.')
    parser.add_argument('--rename_genomes', default=False,  action="store_true",
                        help='Renaming genomes can improve visualization and presentation of results.')
    parser.add_argument('--min_contig_length', default=0,  help='Shorter contigs will be filtered before annotaton.')
    parser.add_argument('--cpus', default=0, help='How many cpus/threads to use (default: all = 0).')
    parser.add_argument('--force', default=False, action="store_true",
                        help='Use force to overwrite previous result files.')
    parser.add_argument('--file_extension', default='.fna', help='When annotating multiple files in a folder, extension'
                                                                 'of the fasta nucleotide files (default: .fna).')
    parser.add_argument('--translation_table', default=11, help='Which translation table to use (default 11).')
    parser.add_argument('--delimiter', default='.', help='Separater character used in feature ids.')
    parser.add_argument('--prefix', default='g', help='Prefix used when renaming genomes (default: g).')
    parser.add_argument('--checkm_dir', default='checkm', help='Dir with the checkm results (default: checkm)')
    parser.add_argument('--gtdbtk_dir', default='gtdbtk', help='Dir with the gtdbtk results (default: gtdbtk).')
    parser.add_argument('--download_database', default=False, action="store_true",
                        help='Download ready-made metaerg database.')
    parser.add_argument('--create_database', default='all', help='Create metaerg database from scratch (default: all, '
                                                                 'to create all components of the database.). Use '
                                                                 'any combination of PVEBRCSA to only create specific '
                                                                 'parts of the database (see README)')
    parser.add_argument('--install_deps', default='', help='Dir for installation of all dependencies '
                                                           '(helper programs). Dependencies will be installed here.')
    parser.add_argument('--path_to_signalp', default='', help='Path to signalp-6.0g.fast.tar.gz.')
    parser.add_argument('--path_to_tmhmm', default='', help='Path to tmhmm-2.0c.Linux.tar.gz.')

    return parser.parse_args()


def load_contigs(genome_name, input_fasta_file, delimiter='.', rename_contigs=False, min_contig_length=500):
    if delimiter in genome_name:
        raise Exception(f'Genome id {genome_name} contains "{delimiter}"; change delimiter with '
                                f'--delimiter [new delimiter] or use --rename_genomes')
    names_done = set()
    contigs = list()
    with fasta.FastaParser(input_fasta_file) as fasta_reader:
        for c in fasta_reader:
            if len(c['seq']) < min_contig_length:
                continue
            contigs.append(c)
            if not rename_contigs:
                if c['id'] in names_done:
                    raise Exception(f'Contig id {c["id"]} not unique. Use --rename_contigs to avoid this problem.')
                names_done.add(c["id"])
    contigs.sort(key=lambda c:len(c['seq']), reverse=True)
    total_length = sum(len(c['seq']) for c in contigs)
    context.log(f'({genome_name}) Loaded {len(contigs)} contigs with total length {total_length:,} from file.')
    if rename_contigs:
        contig_name_mapping_file = context.spawn_file('contig.name.mappings', genome_name)
        context.log(f'({genome_name}) Renaming contigs (see {contig_name_mapping_file})...')
        i = 0
        with open(contig_name_mapping_file, 'w') as mapping_writer:
            for c in contigs:
                new_id = f'{genome_name}.c{i:0>4}'
                mapping_writer.write(f'{c["id"]}\t{new_id}\n')
                c['id'] = new_id
                i += 1
    else:
        for c_id in contigs:
            if delimiter in c_id:
                raise Exception(f'Contig id {c_id} contains "{delimiter}"; change delimiter with '
                                f'--delimiter [new delimiter] or use --rename_contigs')
    contigs = {c['id']: c for c in contigs}
    return contigs


def compute_genome_properties(contig_dict: dict[str, dict], db_connection) -> dict:
    properties = {'# total features':         0,
                  '# proteins':               0,
                  '# ribosomal RNA':          0,
                  '# transfer RNA':           0,
                  '# non-coding RNA':         0,
                  '# retrotransposons':       0,
                  '# CRISPR repeats':         0,
                  '# other repeats':          0,
                  '% coding':                 0.0,
                  '% repeats':                0.0,
                  'mean protein length (aa)': 0.0,}
    contigs:list[dict] = list(contig_dict.values())
    contigs.sort(key=lambda c:len(c['seq']), reverse=True)
    properties['size'] = sum(len(c['seq']) for c in contigs)
    properties['% GC'] = sum((c['seq'].count('G') + c['seq'].count('G') for c in contigs)) / \
                               (properties['size'] - sum((c['seq'].count('N') for c in contigs)))
    cum_size = 0
    for c in contigs:
        cum_size += len(c['seq'])
        if cum_size >+ properties['size'] / 2:
            properties['N50'] = len(c['seq'])
            break

    #feature_data = feature_data.assign(length=feature_data['end'] - feature_data['start'])
    #feature_data_cds = feature_data[feature_data['type'] == 'CDS']
    #feature_data_repeats = feature_data[feature_data['type'].isin(('retrotransposon', 'crispr_repeat', 'repeat'))]
    taxon_counts = Counter()
    for feature in sqlite.read_all_features(db_connection):
        properties['# total features'] += 1
        if feature.type == 'CDS':
            properties['# proteins'] += 1
            properties['% coding'] += feature.length_nt()
            properties['mean protein length (aa)'] += feature.length_aa()
            taxon_counts.update((feature.taxon,))
        elif feature.type == 'rRNA':
            properties['# ribosomal RNA'] += 1
            taxon_counts.update((feature.taxon,))
        elif feature.type == 'tRNA':
            properties['# transfer RNA'] += 1
            taxon_counts.update((feature.taxon,))
        elif feature.type == 'ncRNA':
            properties['# non-coding RNA'] += 1
            taxon_counts.update((feature.taxon,))
        elif feature.type == 'retrotransposon':
            properties['# retrotransposons'] += 1
            properties['% repeats'] += feature.length_nt()
        elif feature.type == 'crispr_repeat':
            properties['# CRISPR repeats'] += 1
            properties['% repeats'] += feature.length_nt()
        elif feature.type == 'repeat':
            properties['# other repeats'] += 1
            properties['% repeats'] += feature.length_nt()
    properties['% coding'] /= properties['size']
    properties['mean protein length (aa)'] /= properties['# proteins']
    properties['% repeats'] /= properties['size']
    properties['% CDS classified to taxon'] = 1 - taxon_counts[''] / taxon_counts.total()
    properties['dominant taxon'] = ''
    properties['% of CDS classified to dominant taxon'] = 0.0
    del taxon_counts['']
    for k, v in taxon_counts.items():
        properties['% of CDS classified to dominant taxon'] = v / taxon_counts.total()
        properties['dominant taxon'] = k
        break
    codon_usage_bias, doubling_time = compute_codon_bias_estimate_doubling_time(db_connection)
    properties['codon usage bias'] = codon_usage_bias
    properties['doubling_time (days)'] = doubling_time
    properties['subsystems'] = functional_gene_configuration.aggregate(db_connection)
    return properties


def annotate_genome(genome_name, input_fasta_file: Path):
    context.log(f'Started annotation of {genome_name}...')
    current_progress = context.parse_metaerg_progress(genome_name)
    if 'visualization=complete' in current_progress:
        context.log(f'({genome_name}) already completed!')
        return
    # (1) prepare sqlite3 database
    db_file = context.spawn_file('annotations.sqlite', genome_name, context.BASE_DIR)
    sqlite.create_db(db_file)
    db_connection = sqlite.connect_to_db(db_file)
    # (2) load sequence data
    contig_dict = load_contigs(genome_name, input_fasta_file, delimiter=context.DELIMITER,
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
    genome_properties = compute_genome_properties(contig_dict, db_connection)
    context.log(f'({genome_name}) Now writing final result as .html for visualization...')
    for html_writer in registry.HTML_WRITER_REGISTRY:
        html_writer(genome_name, db_connection, genome_properties, context.HTML_DIR)
    # (5) update progress
    current_progress = context.parse_metaerg_progress(genome_name)
    current_progress += 'visualization=complete\n'
    context.write_metaerg_progress(genome_name, current_progress)
    context.log(f'({genome_name}) Completed html visualization.')


def write_functional_genes_to_xls():
    excel_file = context.BASE_DIR / 'functional_genes.xls'
    db_files = [context.spawn_file('annotations.sqlite', genome_name, context.BASE_DIR)
                     for genome_name in context.GENOME_NAMES]
    fna_files = [context.spawn_file("fna", genome_name, context.BASE_DIR)
                 for genome_name in context.GENOME_NAMES]
    genome_properties = {}
    for db_file, contig_file in zip(db_files, fna_files):
        genome_name = contig_file.stem
        contig_dict = load_contigs(genome_name, contig_file, delimiter='xxxx')
        db_connection = sqlite.connect_to_db(db_file)
        genome_properties[genome_name] = compute_genome_properties(contig_dict, db_connection)

    all_genome_feature_data = None
    for genome_name, genome_property_hash in genome_properties.items():
        subsystems_df = genome_property_hash['subsystems'].rename(columns={'genes': genome_name})
        try:
            subsystems_df.drop('', level=0, axis=0, inplace=True)
            subsystems_df.drop('secondary-metabolites', level=0, axis=0, inplace=True)
        except Exception:
            pass

        if all_genome_feature_data is None:
            all_genome_feature_data = subsystems_df
        else:
            all_genome_feature_data[genome_name] = subsystems_df[genome_name]
        all_genome_feature_data = all_genome_feature_data.copy()
    del all_genome_feature_data['profiles']
    all_genome_feature_data.to_excel(excel_file)


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
        functional_gene_configuration.init_functional_gene_config()
        if context.PARALLEL_ANNOTATIONS > 1:
            outcomes = {}
            with futures.ProcessPoolExecutor(max_workers=context.PARALLEL_ANNOTATIONS) as executor:
                for genome_name, contig_file in zip(context.GENOME_NAMES, context.CONTIG_FILES):
                    outcomes[genome_name] = executor.submit(annotate_genome, genome_name, contig_file)
            for genome_name, future in outcomes.items():
                try:
                    r = future.result()
                except Exception:
                    context.log(f'({genome_name}) Error while processing:')
                    raise
        else:
            for genome_name, contig_file in zip(context.GENOME_NAMES, context.CONTIG_FILES):
                annotate_genome(genome_name, contig_file)
        context.log('Now writing all-genomes overview html...')
        html_all_genomes.write_html(context.HTML_DIR)
        tmhmm.cleanup(context.BASE_DIR)
        context.log('Now writing excel files with functional genes per genome...')
        write_functional_genes_to_xls()
        context.log(f'Done. Thank you for using metaerg.py {VERSION}')


if __name__ == "__main__":
    main()
