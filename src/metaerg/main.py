import argparse
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor

from metaerg.run_and_read.context import log, spawn_file, init, sorted_annotators, BASE_DIR, RENAME_CONTIGS, \
    MIN_CONTIG_LENGTH, PARAlLEL_ANNOTATIONS, GENOME_NAMES, CONTIG_FILES, HTML_DIR, HTML_WRITER_REGISTRY
from metaerg.run_and_read.data_model import MetaergGenome, FeatureType
from metaerg.html.html_all_genomes import write_html

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
    parser.add_argument('--file_extension', default='.fna', help='When annotating multiple files in a folder, the extension'
                                                                 'of the fasta nucleotide files (default: .fna).')
    parser.add_argument('--translation_table', default=11, help='Which translation table to use (default 11).')
    parser.add_argument('--checkm_dir', default='checkm', help='Dir with the checkm results (default: checkm)')
    parser.add_argument('--gtdbtk_dir', default='gtdbtk', help='Dir with the gtdbtk results (default: gtdbtk).')

    return parser.parse_args()


def annotate_genome(genome_name, input_fasta_file:Path):
    log(f'Now starting to annotate {genome_name}...')
    # (1) create genome, this will load the contigs into memory, they will be filtered and perhaps renamed
    genome = MetaergGenome(genome_name, input_fasta_file, rename_contigs=RENAME_CONTIGS,
                           min_contig_length=MIN_CONTIG_LENGTH)
    if RENAME_CONTIGS:
        contig_name_mappings_file = spawn_file('contig.name.mappings', genome.id)
        genome.write_contig_name_mappings(contig_name_mappings_file)
    # (2) now annotate
    for annotator in sorted_annotators():
        annotator()
    # (3) save results
    log(f'({genome.id}) Now writing to .gbk, .gff, and fasta...')
    genome.write_fasta_files(spawn_file("faa", genome.id, BASE_DIR), target=FeatureType.CDS)
    genome.write_fasta_files(spawn_file("rna.fna", genome.id, BASE_DIR),
                             target=(FeatureType.ncRNA, FeatureType.rRNA, FeatureType.retrotransposon,
                                     FeatureType.tRNA, FeatureType.tmRNA))
    genome.write_gbk_gff(gbk_file=spawn_file("gbk", genome.id, BASE_DIR),
                         gff_file=spawn_file("gff", genome.id))
    # (4) visualize
    genome.compute_properties()
    log(f'({genome.id}) Now writing final result as .html for visualization...')
    for html_writer in HTML_WRITER_REGISTRY:
        html_writer()
    log(f'({genome.id}) Completed html visualization.')


def main():
    print(f'This is metaerg.py {VERSION}')
    init(**parse_arguments().__dict__)

    with ProcessPoolExecutor(max_workers=PARAlLEL_ANNOTATIONS) as executor:
        for genome_name, contig_file in zip(GENOME_NAMES, CONTIG_FILES):
            executor.submit(annotate_genome, contig_file, genome_name)
    log('Now writing all-genomes overview html...')
    write_html(HTML_DIR)
    log(f'Done. Thank you for using metaerg.py {VERSION}')


if __name__ == "__main__":
    main()
