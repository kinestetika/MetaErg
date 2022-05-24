import argparse
import inspect
import os
import shutil
from pathlib import Path
from multiprocessing import cpu_count
from concurrent.futures import ProcessPoolExecutor

import run_and_read.context
from metaerg.run_and_read.data_model import MetaergGenome, FeatureType
from metaerg.run_and_read.context import Executor, annotator_registry
from metaerg.html.abc import html_registry
from metaerg.html.html_all_genomes import HTMLAllGenomesTable
from metaerg import utils


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


def annotate_genome(exec:Executor, genome_name, input_fasta_file:Path):
    run_and_read.execution.log(f'Now starting to annotate {genome_name}...')
    # (1) create genome, this will load the contigs into memory, they will be filtered and perhaps renamed
    genome = MetaergGenome(genome_name, input_fasta_file, rename_contigs=exec.rename_contigs,
                           min_contig_length=exec.min_contig_length)
    if exec.rename_contigs:
        contig_name_mappings_file = exec.spawn_file('contig.name.mappings', genome.id)
        genome.write_contig_name_mappings(contig_name_mappings_file)
    # (2) now annotate
    for annotator_class in annotator_registry:
        annotator = annotator_class(genome, exec)
        annotator.run_and_read()
    # (3) save results
    run_and_read.execution.log(f'({genome.id}) Now writing to .gbk, .gff, and fasta...')
    genome.write_fasta_files(exec.spawn_file("faa", genome.id, exec.base_dir), target=FeatureType.CDS)
    genome.write_fasta_files(exec.spawn_file("rna.fna", genome.id, exec.base_dir),
                             target=(FeatureType.ncRNA, FeatureType.rRNA, FeatureType.retrotransposon,
                                     FeatureType.tRNA, FeatureType.tmRNA))
    genome.write_gbk_gff(gbk_file=exec.spawn_file("gbk", genome.id, exec.base_dir),
                         gff_file=exec.spawn_file("gff", genome.id))
    # (4) visualize
    genome.compute_properties()
    run_and_read.execution.log(f'({genome.id}) Now writing final result as .html for visualization...')
    for html_writer_class in html_registry:
        html_writer = html_writer_class(genome)
        html_writer.write_html()
    run_and_read.execution.log(f'({genome.id}) Completed html visualization.')


def main():
    run_and_read.execution.log(f'This is metaerg.py {VERSION}')
    args = parse_arguments()
    exec = Executor(**args.__dict__)
    exec.prep_environment()

    with ProcessPoolExecutor(max_workers=exec.parallel_annotations) as executor:
        for genome_name, contig_file in zip(exec.genome_names, exec.contig_files):
            executor.submit(annotate_genome, contig_file, genome_name)
        os.chdir(exec.contig_file)
    run_and_read.execution.log('Now writing all-genome overview html...')
    HTMLAllGenomesTable(None, exec).write_html()
    run_and_read.execution.log(f'Done. Thank you for using metaerg.py {VERSION}')


if __name__ == "__main__":
    main()
