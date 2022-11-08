import argparse
import tarfile
import pandas as pd
from pathlib import Path
from concurrent import futures
from hashlib import md5

import metaerg.run_and_read.diamond_and_blastn
from metaerg import context
from metaerg import registry
from metaerg import functional_gene_configuration
from metaerg.datatypes import fasta, gbk
from metaerg.html import html_all_genomes
from metaerg.run_and_read import tmhmm
from metaerg.installation import install_all_helper_programs
from metaerg.run_and_read import *
from metaerg.html import *

VERSION = "2.2.32"


def parse_arguments():
    parser = argparse.ArgumentParser(description='metaerg.py. (C) Marc Strous, Xiaoli Dong 2019, 2022')
    parser.add_argument('--contig_file', help='Fasta nucleotide file of the contigs, or dir that '
                                              'contains multiple fasta nucleotide files.')
    parser.add_argument('--database_dir', required=True,  help='Dir that contains the annotation databases.')
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
    parser.add_argument('--create_database', default='', help='Create metaerg database from scratch.')
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


def compute_genome_properties(contig_dict: dict[str, dict], feature_data: pd.DataFrame) -> dict:
    properties = {}
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

    feature_data = feature_data.assign(length=feature_data['end'] - feature_data['start'])
    feature_data_cds = feature_data[feature_data['type'] == 'CDS']
    feature_data_repeats = feature_data[feature_data['type'].isin(('retrotransposon', 'crispr_repeat', 'repeat'))]

    properties['# proteins'] = len(feature_data_cds.index)
    properties['% coding'] = feature_data_cds['length'].sum() / properties['size']
    properties['mean protein length (aa)'] = feature_data_cds['length'].mean() / 3
    properties['# ribosomal RNA'] = len(feature_data[feature_data['type'] == 'rRNA'].index)
    properties['# transfer RNA'] = len(feature_data[feature_data['type'] == 'tRNA'].index)
    properties['# non-coding RNA'] = len(feature_data[feature_data['type'] == 'ncRNA'].index)
    properties['# retrotransposons'] = len(feature_data[feature_data['type'] == 'retrotransposon'].index)
    properties['# CRISPR repeats'] = len(feature_data[feature_data['type'] == 'crispr_repeat'].index)
    properties['# other repeats'] = len(feature_data[feature_data['type'] == 'repeat'].index)
    properties['% repeats'] = feature_data_repeats['length'].sum() / properties['size']
    properties['# total features'] = len(feature_data.index)
    taxon_counts = dict(feature_data_cds.taxon.value_counts(normalize=True))
    properties['% CDS classified to taxon'] = 1 - taxon_counts['']
    properties['dominant taxon'] = ''
    del taxon_counts['']
    for k, v in taxon_counts.items():
        properties['% of CDS classified to dominant taxon'] = v / properties['% CDS classified to taxon']
        properties['dominant taxon'] = k
        break
    properties['subsystems'] = functional_gene_configuration.aggregate(feature_data)
    return properties


def annotate_genome(genome_name, input_fasta_file: Path):
    context.log(f'Started annotation of {genome_name}...')
    # (1) prepare dataframe
    feature_data = pd.DataFrame(columns=context.DATAFRAME_COLUMNS)
    # (2) load sequence data
    contig_dict = load_contigs(genome_name, input_fasta_file, delimiter=context.DELIMITER,
                               min_contig_length=context.MIN_CONTIG_LENGTH, rename_contigs=context.RENAME_CONTIGS)
    # (3) now annotate
    for annotator in context.sorted_annotators():
        feature_data = annotator(genome_name, contig_dict, feature_data)
    # feather_file = context.spawn_file("all_genes.feather", genome_name, context.BASE_DIR)
    # feature_data = pd.read_feather(feather_file)
    # feature_data = feature_data.set_index('id', drop=False)

    # (4) save results
    context.log(f'({genome_name}) Now writing annotations to .fasta, .gbk and .feather...')
    faa_file = context.spawn_file("faa", genome_name, context.BASE_DIR)
    rna_file = context.spawn_file("rna.fna", genome_name, context.BASE_DIR)
    gbk_file = context.spawn_file("gbk", genome_name, context.BASE_DIR)
    fna_file = context.spawn_file("fna", genome_name, context.BASE_DIR)
    feather_file = context.spawn_file("all_genes.feather", genome_name, context.BASE_DIR)
    fasta.write_features_to_fasta(feature_data, faa_file, targets=('CDS',))
    fasta.write_features_to_fasta(feature_data, rna_file, targets=('rRNA tRNA tmRNA ncRNA retrotransposon'.split()))
    fasta.write_contigs_to_fasta(contig_dict, fna_file)
    with open(gbk_file, 'w') as gbk_writer:
        gbk.gbk_write_genome(gbk_writer, contig_dict, feature_data)
    feature_data = feature_data.reset_index(drop=True)
    feature_data.to_feather(feather_file)
    # (5) visualize
    genome_properties = compute_genome_properties(contig_dict, feature_data)
    context.log(f'({genome_name}) Now writing final result as .html for visualization...')
    for html_writer in registry.HTML_WRITER_REGISTRY:
        html_writer(genome_name, feature_data, genome_properties, context.HTML_DIR)
    context.log(f'({genome_name}) Completed html visualization.')


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
        metaerg.run_and_read.diamond_and_blastn.compile_databases()
    elif context.METAERG_MODE == context.METAERG_MODE_INSTALL_DEPS:
        install_all_helper_programs(context.BIN_DIR, context.DATABASE_DIR, context.PATH_TO_SIGNALP,
                                    context.PATH_TO_TMHMM)
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
    context.log(f'Done. Thank you for using metaerg.py {VERSION}')


if __name__ == "__main__":
    main()
