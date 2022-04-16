import argparse
import os
from pathlib import Path

from Bio import SeqIO
from BCBio import GFF

from metaerg import databases
from metaerg import predict
from metaerg import utils

VERSION = "2.0.15"


def parse_arguments():
    parser = argparse.ArgumentParser(description='metaerg.py. (C) Marc Strous, Xiaoli Dong 2019, 2021')
    parser.add_argument('--contig_file', required=True,  help='fasta nucleotide file of the contigs')
    parser.add_argument('--database_dir', required=True,  help='dir that contains the annotation databases')
    args = parser.parse_args()
    return args


def create_output_dir(parent_dir:Path):
    output_dir = Path(parent_dir, "metaerg")
    if output_dir.exists():
        utils.log('Warning: may overwrite existing output files...')
        if output_dir.is_file():
            utils.log(f'Expected folder at {output_dir}, found regular file, crash! Delete this fiel first')
            exit(1)
    else:
        os.mkdir(output_dir)
    return output_dir


def get_tmp_file(tmp_dir):
    return os.path.join()


def main():
    utils.log(f'This is metaerg.py {VERSION}')
    args = parse_arguments()

    # (1) set and validate database dir
    dbdir = Path(args.database_dir)
    databases.DBDIR = dbdir
    if not databases.does_db_appear_valid():
        utils.log(f'Metaerg database at "{dbdir}" appears missing or invalid.')
        exit(1)

    # (2) set and validate fasta .fna file, mag (genome) name,
    input_fasta_file = Path(args.contig_file).absolute()
    if not input_fasta_file.exists() or input_fasta_file.is_dir():
        utils.log(f'Input file "{input_fasta_file}" is missing or not a valid file. Expecting a nt fasta file.')
        exit(1)
    mag_name = input_fasta_file.stem

    # (3) set output files
    output_dir = create_output_dir(input_fasta_file.parent)
    os.chdir(output_dir)
    gbk_file = predict.spawn_file("gbk", mag_name)
    gff_file = predict.spawn_file("gff", mag_name)

    # (4) Filter and load contigs
    filtered_fasta_file = predict.spawn_file('filtered.fna', mag_name)
    utils.create_filtered_contig_fasta_file(fasta_file_in=input_fasta_file, fasta_file_out=filtered_fasta_file)
    contig_dict = SeqIO.to_dict(SeqIO.parse(filtered_fasta_file, "fasta"))
    for contig in contig_dict.values():
        contig.annotations['molecule_type'] = 'DNA'


    # Feature prediction

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
                       predict.predict_transmembrane_helixes,
                       predict.predict_signal_peptides,
                       predict.predict_functions_with_antismash,
                       predict.create_files
                       ):
        prediction(mag_name, contig_dict)
        SeqIO.write(contig_dict.values(), gbk_file, "genbank")
        with open(gff_file, "w") as gff_handle:
             GFF.write(contig_dict.values(), gff_handle)

    utils.log(f'Done. Thank you for using metaerg.py {VERSION}')


if __name__ == "__main__":
    main()
