import argparse
import inspect
import os
import shutil
from pathlib import Path
from multiprocessing import cpu_count
from concurrent.futures import ProcessPoolExecutor

from Bio import SeqIO
from BCBio import GFF

from metaerg import run_and_read
from metaerg.run_and_read.data_model import MetaergGenome
from metaerg.run_and_read.abc import Annotator, ExecutionEnvironment
from metaerg.run_and_read.antismash import Antismash
from metaerg.run_and_read.aragorn import Aragorn
from metaerg.run_and_read.canthyd import CantHyd
from metaerg.run_and_read.cdd import CDD
from metaerg.run_and_read.cmscan import CMScan
from metaerg.run_and_read.diamond_and_blastn import DiamondAndBlastN
from metaerg.run_and_read.ltr_harvest import LTRHarvest
from metaerg.run_and_read.minced import Minced
from metaerg.run_and_read.prodigal import Prodigal
from metaerg.run_and_read.repeat_masker import RepeatMasker
from metaerg.run_and_read.signalp import SignalP
from metaerg.run_and_read.tmhmm import TMHMM
from metaerg import databases
from metaerg import predict
from metaerg import utils
from metaerg import visualization
from metaerg import subsystems


VERSION = "2.1.0"

def parse_arguments():
    parser = argparse.ArgumentParser(description='metaerg.py. (C) Marc Strous, Xiaoli Dong 2019, 2022')
    parser.add_argument('--contig_file', required=True,  help='Fasta nucleotide file of the contigs, or dir that '
                                                              'contains multiple fasta nucleotide files.')
    parser.add_argument('--database_dir', required=True,  help='Dir that contains the annotation databases.')
    parser.add_argument('--rename_contigs', default=False,  action=argparse.BooleanOptionalAction,
                        help='Renaming contigs can improve visualization. It will be enforced when original '
                             'contig names are long')
    parser.add_argument('--rename_mags', default=False,  action=argparse.BooleanOptionalAction,
                        help='Renaming mags can improve visualization. It will be enforced when original '
                             'file names are long (>7 chars)')
    parser.add_argument('--min_contig_length', default=0,  help='Shorter contigs will be filtered before annotaton.')
    parser.add_argument('--cpus', default=0, help='How many cpus/threads to use (default: all = 0).')
    parser.add_argument('--file_extension', default='.fna', help='When annotating multiple files in a folder, the extension'
                                                                 'of the fasta nucleotide files (default: .fna).')
    parser.add_argument('--translation_table', default=11, help='Which translation table to use (default 11).')
    parser.add_argument('--checkm_dir', default='checkm', help='Dir with the checkm results (default: checkm)')
    parser.add_argument('--gtdbtk_dir', default='gtdbtk', help='Dir with the gtdbtk results (default: gtdbtk).')

    return parser.parse_args()




def create_temp_dir(parent_dir:Path):
    temp_dir = Path(parent_dir, "temp")
    if temp_dir.exists():
        utils.log('Warning: may overwrite existing temp files...')
        if temp_dir.is_file():
            utils.log(f'Expected folder at {temp_dir}, found regular file, crash! Delete this file first')
            exit(1)
    else:
        os.mkdir(temp_dir)

    return temp_dir


def filter_and_rename_contigs(mag_name, input_fasta_file, rename_contigs, min_length):
    filtered_fasta_file = predict.spawn_file('filtered.fna', mag_name)
    contig_dict = SeqIO.to_dict(SeqIO.parse(input_fasta_file, "fasta"))
    if not rename_contigs: # check if contig names are short enough
        for contig_name in contig_dict.keys():
            if len(contig_name) > 5:
                utils.log(f'Contig name "{contig_name}" is too long for visualization, will rename contigs...')
                rename_contigs = True
                break
    if rename_contigs:
        rename_txt = ' renaming contigs,'
    else:
        rename_txt = ''
    contig_name_mappings = {}
    utils.log(f'Filtering contigs for length, removing gaps,{rename_txt} replacing non-IUPAC bases with N, '
              'capitalizing...')
    i = 0
    filtered_contig_dict = {}
    with open(filtered_fasta_file, 'w') as fasta_writer:
        for contig in contig_dict.values():
            if len(contig) < min_length:
                continue
            filtered_contig = utils.filter_seq(contig)
            if rename_contigs:
                new_id = f'{mag_name}.c{i:0>4}'
                contig_name_mappings[new_id] = contig.id
                filtered_contig.id = new_id
                filtered_contig.description = filtered_contig.id
            i += 1
            filtered_contig.annotations['molecule_type'] = 'DNA'
            SeqIO.write(filtered_contig, fasta_writer, "fasta")
            filtered_contig_dict[filtered_contig.id] = filtered_contig
    utils.log(f'Wrote {i} contigs of length >{min_length} nt to {filtered_fasta_file}')

    if rename_contigs:
        contig_name_mappings_file = predict.spawn_file('contig.name.mappings', mag_name)
        with open(contig_name_mappings_file, 'w') as mapping_writer:
            for key, value in contig_name_mappings.items():
                mapping_writer.write(f'{key}\t{value}\n')

    return filtered_contig_dict


def annotate_genome_2(exec:ExecutionEnvironment, genome_name, input_fasta_file:Path, rename_contigs=True, min_length=0):
    # (1) set and validate fasta .fna file, load data
    utils.log(f'Now starting to annotate {genome_name}...')
    if not input_fasta_file.exists() or input_fasta_file.is_dir():
        utils.log(f'Input file "{input_fasta_file}" is missing or not a valid file. Expecting a nt fasta file.')
        return
    genome = MetaergGenome(genome_name, input_fasta_file, rename_contigs=rename_contigs, min_contig_length=min_length)
    if rename_contigs:
        contig_name_mappings_file = exec.spawn_file('contig.name.mappings', genome.id)
        genome.write_contig_name_mappings(contig_name_mappings_file)

    # (2) prep temp output files, note that these file paths are relative to the current working dir
    working_directory = input_fasta_file.parent # eventually: os.getcwd()
    temp_dir = create_temp_dir(working_directory)  # (doesn't overwrite)
    os.chdir(temp_dir)
    # (3) Filter, rename and load contigs

    inspect.getmembers(run_and_read, inspect.isclass)
    for annotator in (Minced(genome, exec),):
        annotator.run_and_read()



def annotate_genome(input_fasta_file:Path, mag_name, rename_contigs=True, min_length=0):
    # (1) set and validate fasta .fna file, mag (genome) name,
    working_directory = input_fasta_file.parent # eventually: os.getcwd()
    utils.log(f'Now starting to annotate {mag_name}...')
    if not input_fasta_file.exists() or input_fasta_file.is_dir():
        utils.log(f'Input file "{input_fasta_file}" is missing or not a valid file. Expecting a nt fasta file.')
        return

    # (2) prep temp output files, note that these file paths are relative to the current working dir
    temp_dir = create_temp_dir(working_directory)  # (doesn't overwrite)
    os.chdir(temp_dir)
    gbk_file = predict.spawn_file("gbk", mag_name)
    gff_file = predict.spawn_file("gff", mag_name)

    # (3) Filter, rename and load contigs
    contig_dict = filter_and_rename_contigs(mag_name, input_fasta_file, rename_contigs, min_length)

    # (4) Create empty subsystems hash
    subsystem_hash = subsystems.get_empty_subsystem_hash()

    # (5) Feature prediction

    ################
    # for debugging individual predicton tools, first run the pipeline,
    # then, uncomment the following lines:
    ################
    # contig_dict = {}
    # with open(gbk_file) as handle:
    #     for gb_record in SeqIO.parse(handle, "genbank"):
    #         contig_dict[gb_record.id] = gb_record
    #         #print(gb_record.id, len(gb_record.features))
    # features.predict_retrotransposons_with_ltrharvest(fasta_file, contig_dict)
    # SeqIO.write(contig_dict.values(), gbk_file, "genbank")
    # with open(gff_file, "w") as gff_handle:
    #     GFF.write(contig_dict.values(), gff_handle)
    # exit(0)

    for prediction in (predict.predict_crisprs_with_minced,
                       predict.predict_trnas_with_aragorn,
                       predict.predict_non_coding_rna_features_with_infernal_cmscan,  # multithreaded by splitting input
                       predict.predict_retrotransposons_with_ltrharvest,
                       predict.predict_tandem_repeats_with_trf,
                       predict.predict_remaining_repeats_with_repeatmasker,  # may use >1 threads
                       predict.predict_coding_sequences_with_prodigal,
                       predict.create_ids,
                       predict.write_gene_files,
                       predict.predict_functions_and_taxa_with_diamond, # may use >1 threads
                       predict.predict_functions_and_taxa_with_blastn,
                       predict.predict_functions_with_cdd,  # ??? multithreaded by splitting input
                       predict.predict_functions_with_antismash,
                       predict.predict_hydrocarbon_genes_with_canthyd,
                       predict.predict_transmembrane_helixes,
                       predict.predict_signal_peptides,  # multithreaded by splitting input
                       predict.predict_subsystems,
                       predict.write_databases
                       ):
        prediction(mag_name, contig_dict, subsystem_hash)
        SeqIO.write(contig_dict.values(), gbk_file, "genbank")
        with open(gff_file, "w") as gff_handle:
             GFF.write(contig_dict.values(), gff_handle)
    mag_antismash_result_dir = predict.spawn_file('antismash', mag_name).absolute()

    # (7) write flat text result files
    utils.log(f'({mag_name}) Now writing to .gbk, .gff, and fasta...')
    os.chdir(working_directory)
    faa_file = predict.spawn_file("faa", mag_name)
    fna_file = predict.spawn_file("rna.fna", mag_name)
    gbk_file = predict.spawn_file("gbk", mag_name)
    gff_file = predict.spawn_file("gff", mag_name)
    with open(faa_file, 'w') as faa_handle, open(fna_file, 'w') as fna_handle:
        for contig in contig_dict.values():
            for f in contig.features:
                if f.type == 'CDS':
                    feature_seq = utils.pad_seq(f.extract(contig)).translate(table=utils.TRANSLATION_TABLE)[:-1]
                    feature_seq.id = utils.get_feature_qualifier(f, 'id')
                    feature_seq.description = f'{feature_seq.id} {visualization.make_feature_short_description(f)}'
                    SeqIO.write(feature_seq, faa_handle, "fasta")
                elif f.type in ['tRNA', 'rRNA', 'ncRNA']:
                    feature_seq = f.extract(contig)
                    feature_seq.id = utils.get_feature_qualifier(f, 'id')
                    feature_seq.description = f'{feature_seq.id} {visualization.make_feature_short_description(f)}'
                    SeqIO.write(feature_seq, fna_handle, "fasta")
    SeqIO.write(contig_dict.values(), gbk_file, "genbank")
    with open(gff_file, "w") as gff_handle:
        GFF.write(contig_dict.values(), gff_handle)
    # (8) write html result files
    utils.log(f'({mag_name}) Now writing final result as .html for visualization...')
    mag_html_dir = predict.spawn_file("html", mag_name)
    shutil.rmtree(mag_html_dir, ignore_errors=True)
    mag_html_dir.mkdir()
    os.chdir(mag_html_dir)
    if mag_antismash_result_dir.exists():
        shutil.copytree(mag_antismash_result_dir, 'antismash')
        #shutil.rmtree('antismash', ignore_errors=True)
        #shutil.move(mag_antismash_result_dir.name, 'antismash')
    genome_stats = predict.compile_genome_stats(mag_name, contig_dict)
    visualization.html_write_all(mag_name, genome_stats, contig_dict, predict.BLAST_RESULTS, subsystem_hash)


def main():
    utils.log(f'This is metaerg.py {VERSION}')
    args = parse_arguments()
    get_available_prereqs()
    utils.TRANSLATION_TABLE = args.translation_table
    # (1) set and validate database dir
    dbdir = Path(args.database_dir)
    databases.DBDIR = dbdir
    if databases.does_db_appear_valid():
        utils.log(f'Metaerg database at "{dbdir}" appears valid.')
    else:
        utils.log(f'Metaerg database at "{dbdir}" appears missing or invalid.')
        exit(1)
    databases.load_descriptions_taxonomy_cdd()
    # (2) prep subsystems
    subsystems.prep_subsystems()
    # (3) determine # of threads available and how many we're using
    cpus_available = cpu_count()
    cpus_used = int(args.cpus)
    if cpus_used > 0:
        cpus_used = min(cpus_used, cpus_available)
    else:
        cpus_used = cpus_available
    utils.log(f'Detected {cpus_available} available threads/cpus, will use {cpus_used}.')
    # (4) determine how many genomes we're annotating and annotate...
    contig_file = Path(args.contig_file).absolute()
    if contig_file.is_dir():
        contig_files = [x.absolute() for x in sorted(contig_file.glob(f'*{args.file_extension}'))]
        if not len(contig_files):
            utils.log(f'Did not find any contig files with extension "{args.file_extension}" in dir "{contig_file}"')
            exit(1)
        predict.MULTI_MODE = True
        predict.THREADS_PER_GENOME = max(1, int(cpus_used / len(contig_files)))
        utils.log(f'Ready to annotate {len(contig_files)} genomes in dir "{contig_file}" with '
                  f'{predict.THREADS_PER_GENOME} threads per genome.')
        with ProcessPoolExecutor(max_workers=cpus_used) as executor:
            count = 0
            with open(Path(contig_file, 'mag.name.mapping.txt'), 'w') as mapping_file:
                for f in contig_files:
                    new_mag_name = f'g{count:0>4}'
                    utils.log(f'Now submitting {f} for annotation as {new_mag_name}...')
                    executor.submit(annotate_genome, f, new_mag_name, True, args.min_contig_length)
                    mapping_file.write(f'{new_mag_name}\t{f.name}\n')
                    count += 1
        os.chdir(contig_file)
        visualization.html_write_mag_table(open(Path(contig_file, 'html', 'index.html'), 'w'),
                                           contig_file, gtdbtk_dir=args.gtdbtk_dir, checkm_dir=args.checkm_dir)
    else:
        # annotate a single genome
        utils.log(f'Ready to annotate genome in nucleotide fasta file "{contig_file}".')
        predict.THREADS_PER_GENOME = cpus_used
        mag_name = contig_file.stem
        if args.rename_mags or len(mag_name) > 7:
            new_mag_name = 'g00001'
            utils.log(f'Genome name "{mag_name}" is too long for visualization, genome renamed to {new_mag_name}...')
            mag_name = new_mag_name
        annotate_genome(contig_file, mag_name, rename_contigs=args.rename_contigs, min_length=args.min_contig_length)
    utils.log(f'Done. Thank you for using metaerg.py {VERSION}')


if __name__ == "__main__":
    main()
