import argparse
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor

from metaerg import context
from metaerg.data_model import MetaergGenome, FeatureType
from metaerg.html.html_all_genomes import write_html
from metaerg import registry
from metaerg.run_and_read import *

VERSION = "2.1.0"


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
    parser.add_argument('--checkm_dir', default='checkm', help='Dir with the checkm results (default: checkm)')
    parser.add_argument('--gtdbtk_dir', default='gtdbtk', help='Dir with the gtdbtk results (default: gtdbtk).')
    parser.add_argument('--create_db', default='', help='Create/download metaerg database. '
                                                        'Use PVEBRCS to create all...')

    return parser.parse_args()


def annotate_genome(genome_name, input_fasta_file:Path):
    context.log(f'Now starting to annotate {genome_name}...')
    # (1) create genome, this will load the contigs into memory, they will be filtered and perhaps renamed
    genome = MetaergGenome(genome_name, input_fasta_file, rename_contigs=context.RENAME_CONTIGS,
                           min_contig_length=context.MIN_CONTIG_LENGTH)
    if context.RENAME_CONTIGS:
        contig_name_mappings_file = context.spawn_file('contig.name.mappings', genome.id)
        genome.write_contig_name_mappings(contig_name_mappings_file)
    # (2) now annotate
    for annotator in context.sorted_annotators():
        annotator()
    # (3) save results
    context.log(f'({genome.id}) Now writing to .gbk, .gff, and fasta...')
    genome.write_fasta_files(context.spawn_file("faa", genome.id, context.BASE_DIR), target=FeatureType.CDS)
    genome.write_fasta_files(context.spawn_file("rna.fna", genome.id, context.BASE_DIR),
                             target=(FeatureType.ncRNA, FeatureType.rRNA, FeatureType.retrotransposon,
                                     FeatureType.tRNA, FeatureType.tmRNA))
    genome.write_gbk_gff(gbk_file=context.spawn_file("gbk", genome.id, context.BASE_DIR),
                         gff_file=context.spawn_file("gff", genome.id))
    # (4) visualize
    genome.compute_properties()
    context.log(f'({genome.id}) Now writing final result as .html for visualization...')
    for html_writer in registry.HTML_WRITER_REGISTRY:
        html_writer()
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
        with ProcessPoolExecutor(max_workers=context.PARAlLEL_ANNOTATIONS) as executor:
            for genome_name, contig_file in zip(context.GENOME_NAMES, context.CONTIG_FILES):
                executor.submit(annotate_genome, contig_file, genome_name)
        context.log('Now writing all-genomes overview html...')
        write_html(context.HTML_DIR)
    context.log(f'Done. Thank you for using metaerg.py {VERSION}')


if __name__ == "__main__":
    main()
