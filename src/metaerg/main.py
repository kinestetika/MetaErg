import argparse
import os
from pathlib import Path

from Bio import SeqIO
from BCBio import GFF

from metaerg import databases
from metaerg import predict
from metaerg import utils
from metaerg import visualization
from metaerg import subsystems

VERSION = "2.0.15"


def parse_arguments():
    parser = argparse.ArgumentParser(description='metaerg.py. (C) Marc Strous, Xiaoli Dong 2019, 2021')
    parser.add_argument('--contig_file', required=True,  help='fasta nucleotide file of the contigs')
    parser.add_argument('--database_dir', required=True,  help='dir that contains the annotation databases')
    parser.add_argument('--rename_contigs', default=False,  action=argparse.BooleanOptionalAction,
                        help='renaming contigs improves visualization. It will be enforced when original '
                             'contig names are long')
    parser.add_argument('--rename_mags', default=False,  action=argparse.BooleanOptionalAction,
                        help='renaming mags improves visualization. It will be enforced when original '
                             'file names are long (>7 chars)')
    parser.add_argument('--min_contig_length', default=0,  help='shorter contigs will be filtered before annotaton.')

    args = parser.parse_args()
    return args


def create_output_dir(parent_dir:Path):
    output_dir = Path(parent_dir, "temp")
    if output_dir.exists():
        utils.log('Warning: may overwrite existing temp files...')
        if output_dir.is_file():
            utils.log(f'Expected folder at {output_dir}, found regular file, crash! Delete this file first')
            exit(1)
    else:
        os.mkdir(output_dir)

    html_dir = Path(parent_dir, "html")
    if html_dir.exists():
        utils.log('Warning: may overwrite existing html files...')
        if html_dir.is_file():
            utils.log(f'Expected folder at {html_dir}, found regular file, crash! Delete this file first')
            exit(1)
    else:
        os.mkdir(html_dir)
    return output_dir, html_dir


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


def annotate_genome(contig_file, genome_id=0, rename_contigs=True, rename_mags=True, min_length=0):
    # (2) set and validate fasta .fna file, mag (genome) name,
    input_fasta_file = Path(contig_file).absolute()
    working_directory = input_fasta_file.parent # eventually: os.getcwd()
    if not input_fasta_file.exists() or input_fasta_file.is_dir():
        utils.log(f'Input file "{input_fasta_file}" is missing or not a valid file. Expecting a nt fasta file.')
        exit(1)
    if rename_mags:
        mag_name = f'g{genome_id:0>4}'
    else:
        mag_name = input_fasta_file.stem
        if len(mag_name) > 5:
            new_mag_name = f'g{genome_id:0>4}'
            utils.log(f'Genome name "{mag_name}" is too long for visualization, genome renamed to {new_mag_name}...')
            mag_name = new_mag_name

    # (3) set output files, note that these file paths are relative to the current working dir
    output_dir, html_dir = create_output_dir(working_directory)
    os.chdir(output_dir)
    gbk_file = predict.spawn_file("gbk", mag_name)
    gff_file = predict.spawn_file("gff", mag_name)
    faa_file = predict.spawn_file("faa", mag_name)
    fna_file = predict.spawn_file("rna.fna", mag_name)

    # (4) Filter, rename and load contigs
    contig_dict = filter_and_rename_contigs(mag_name, input_fasta_file, rename_contigs, min_length)

    # (5) Create empty subsystems hash
    subsystem_hash = subsystems.get_empty_subsystem_hash()

    # (6) Feature prediction

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
                       predict.predict_non_coding_rna_features_with_infernal,
                       predict.predict_retrotransposons_with_ltrharvest,
                       predict.predict_tandem_repeats_with_trf,
                       predict.predict_remaining_repeats_with_repeatmasker,
                       predict.predict_coding_sequences_with_prodigal,
                       predict.create_ids,
                       databases.load_descriptions_taxonomy_cdd,
                       predict.write_gene_files,
                       predict.predict_functions_and_taxa_with_diamond,
                       predict.predict_functions_and_taxa_with_blastn,
                       predict.predict_functions_with_cdd,
                       predict.predict_functions_with_antismash,
                       predict.predict_transmembrane_helixes,
                       predict.predict_signal_peptides,
                       predict.predict_subsystems,
                       predict.write_databases
                       ):
        prediction(mag_name, contig_dict, subsystem_hash)
        SeqIO.write(contig_dict.values(), gbk_file, "genbank")
        with open(gff_file, "w") as gff_handle:
             GFF.write(contig_dict.values(), gff_handle)

    # (7) write main result files
    utils.log("Now writing final result files...")
    os.chdir(working_directory)
    with open(faa_file, 'w') as faa_handle, open(fna_file, 'w') as fna_handle:
        for contig in contig_dict.values():
            for f in contig.features:
                if f.type == 'CDS':
                    feature_seq = utils.pad_seq(f.extract(contig)).translate(table=predict.TRANSLATION_TABLE)[:-1]
                    feature_seq.id = utils.get_feature_qualifier(f, 'id')
                    feature_seq.description = f'{feature_seq.id} {visualization.make_feature_short_description(f)}'
                    if '*' in feature_seq:
                        utils.log(f'Warning, internal stop codon(s) in CDS {utils.get_feature_qualifier(f, "id")} {feature_seq.seq}')
                    SeqIO.write(feature_seq, faa_handle, "fasta")
                elif f.type in ['tRNA', 'rRNA', 'ncRNA']:
                    feature_seq = f.extract(contig)
                    feature_seq.id = utils.get_feature_qualifier(f, 'id')
                    feature_seq.description = f'{feature_seq.id} {visualization.make_feature_short_description(f)}'
                    SeqIO.write(feature_seq, fna_handle, "fasta")
    SeqIO.write(contig_dict.values(), gbk_file, "genbank")
    with open(gff_file, "w") as gff_handle:
        GFF.write(contig_dict.values(), gff_handle)
    os.chdir(html_dir)
    genome_stats = predict.compile_genome_stats(mag_name, contig_dict, subsystem_hash)
    visualization.html_save_all(mag_name, genome_stats, contig_dict, predict.BLAST_RESULTS, subsystem_hash)
    utils.log(f'Done. Thank you for using metaerg.py {VERSION}')


def main():
    utils.log(f'This is metaerg.py {VERSION}')
    args = parse_arguments()
    # (1) set and validate database dir
    dbdir = Path(args.database_dir)
    databases.DBDIR = dbdir
    if not databases.does_db_appear_valid():
        utils.log(f'Metaerg database at "{dbdir}" appears missing or invalid.')
        exit(1)
    # (2) prep subsystems
    subsystems.prep_subsystems()
    # (3) annotate..
    annotate_genome(args.contig_file, rename_contigs=args.rename_contigs, rename_mags=args.rename_mags,
                    min_length=args.min_contig_length)

if __name__ == "__main__":
    main()
