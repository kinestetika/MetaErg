import argparse
import pandas as pd
from pathlib import Path

from metaerg import context
from metaerg import registry
from metaerg.datatypes import fasta

VERSION = "2.2.0"
DATAFRAME_COLUMNS = 'genome contig id start end strand type inference subsystems descr taxon notes seq antismash ' \
                    'signal_peptide tmh tmh_topology blast cdd'.split()
DATAFRAME_DATATYPES = str, str, str, int, int, int, "category", str, str, str, str, str, str, str, \
                      str, int, str, str, str, str
DATAFRAME_DATATYPES = {k: v for k, v in zip(DATAFRAME_COLUMNS, DATAFRAME_DATATYPES)}


def parse_arguments():
    parser = argparse.ArgumentParser(description='metaerg.py. (C) Marc Strous, Xiaoli Dong 2019, 2022')
    parser.add_argument('--contig_file', required=True,  help='Fasta nucleotide file of the contigs, or dir that '
                                                              'contains multiple fasta nucleotide files.')
    parser.add_argument('--database_dir', required=True,  help='Dir that contains the annotation databases.')
    parser.add_argument('--rename_contigs', default=False,  action="store_true",
                        help='Renaming contigs can improve visualization. It will be enforced when original '
                             'contig names are long')
    parser.add_argument('--rename_genomes', default=False,  action="store_true",
                        help='Renaming genomes can improve visualization.')
    parser.add_argument('--min_contig_length', default=0,  help='Shorter contigs will be filtered before annotaton.')
    parser.add_argument('--cpus', default=0, help='How many cpus/threads to use (default: all = 0).')
    parser.add_argument('--force', default=False, action="store_true",
                        help='Use force to overwrite previous result files.')
    parser.add_argument('--file_extension', default='.fna', help='When annotating multiple files in a folder, extension'
                                                                 'of the fasta nucleotide files (default: .fna).')
    parser.add_argument('--translation_table', default=11, help='Which translation table to use (default 11).')
    parser.add_argument('--delimiter', default='.', help='Separater character used in feature ids.')
    parser.add_argument('--checkm_dir', default='checkm', help='Dir with the checkm results (default: checkm)')
    parser.add_argument('--gtdbtk_dir', default='gtdbtk', help='Dir with the gtdbtk results (default: gtdbtk).')
    parser.add_argument('--create_db', default='', help='Create/download metaerg database. '
                                                        'Use PVEBRCS to create all...')

    return parser.parse_args()


def load_contigs(genome_name, input_fasta_file):
    with fasta.FastaParser(input_fasta_file) as fasta_reader:
        contigs = [c for c in fasta_reader if len(c) > context.MIN_CONTIG_LENGTH]
    contigs.sort(key=len, reverse=True)
    contigs = {c.id: c for c in contigs}
    total_length = sum(len(c['seq']) for c in contigs.values())
    context.log(f'({genome_name}) Loaded {len(contigs)} contigs with total length {total_length} from file.')
    if context.RENAME_CONTIGS:
        contig_name_mapping_file = context.spawn_file('contig.name.mappings', genome_name)
        context.log(f'({genome_name}) Renaming contigs (see {contig_name_mapping_file})...')
        i = 0
        with open(contig_name_mapping_file, 'w') as mapping_writer:
            for c in contigs.values():
                new_id = f'{genome_name}.c{i:0>4}'
                mapping_writer.write(f'{c.id}\t{new_id}\n')
                c.id = new_id
                i += 1
    return contigs


def validate_ids(genome_name, contig_hash):
    if context.DELIMITER in genome_name:
        raise Exception(f'Genome id {genome_name} contains "{context.DELIMITER}"; change using --delimiter')
    for c_id in contig_hash.keys():
        if context.DELIMITER in c_id:
            raise Exception(f'Contig id {c_id} contains "{context.DELIMITER}"; change using --delimiter')


def annotate_genome(genome_name, input_fasta_file: Path):
    context.log(f'Started annotation of {genome_name}...')
    # (1) prepare dataframe
    feature_data = pd.DataFrame(columns=DATAFRAME_COLUMNS)
    feature_data = feature_data.astype(DATAFRAME_DATATYPES)
    # (2) load sequence data
    contig_dict = load_contigs(genome_name, input_fasta_file)
    validate_ids(genome_name, contig_dict)
    # (3) now annotate
    for annotator in context.sorted_annotators():
        feature_data = annotator(genome_name, contig_dict, feature_data)
    # (4) save results
    context.log(f'({genome_name}) Now writing to .gbk, .gff, and fasta...')
    faa_file = context.spawn_file("faa", genome.id, context.BASE_DIR)
    rna_file = context.spawn_file("rna.fna", genome.id, context.BASE_DIR)
    bioparsers.write_genome_to_fasta_files(genome, faa_file, targets=(FeatureType.CDS,))
    bioparsers.write_genome_to_fasta_files(genome, rna_file, targets=RNA_FEATURES)
    # (4) visualize
    genome.compute_properties()
    context.log(f'({genome.id}) Now writing final result as .html for visualization...')
    for html_writer in registry.HTML_WRITER_REGISTRY:
        html_writer(genome, Path(context.BASE_DIR, 'html'))
    context.log(f'({genome.id}) Completed html visualization.')


def main():
    print(f'This is metaerg.py {VERSION}')
    context.init(**parse_arguments().__dict__)
    # print(',\n'.join(f'{k}={eval("context." + k)}' for k in dir(context) if not k.startswith('_')))

    if context.CREATE_DB_TASKS:
        context.log(f'Creating/installing/downloading metaerg databases. Tasks: {context.CREATE_DB_TASKS}; '
                    f'force: {context.FORCE}.')
        for db_installer in registry.DATABASE_INSTALLER_REGISTRY:
            db_installer()
    else:
        for genome_name, contig_file in zip(context.GENOME_NAMES, context.CONTIG_FILES):
            annotate_genome(genome_name, contig_file)
        #with ProcessPoolExecutor(max_workers=context.PARALLEL_ANNOTATIONS) as executor:
        #    for genome_name, contig_file in zip(context.GENOME_NAMES, context.CONTIG_FILES):
        #        executor.submit(annotate_genome, genome_name, contig_file)
        # context.log('Now writing all-genomes overview html...')
        # write_html(context.HTML_DIR)
    context.log(f'Done. Thank you for using metaerg.py {VERSION}')


if __name__ == "__main__":
    main()
