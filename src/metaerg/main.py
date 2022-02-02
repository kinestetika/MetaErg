import argparse
import os
import utils
import features
from pathlib import Path
from Bio import SeqIO
from BCBio import GFF

from src.metaerg import databases

VERSION = 0.0


def parse_arguments():
    parser = argparse.ArgumentParser(description='metaerg.py. (C) Marc Strous, Xiaoli Dong 2019, 2021')
    parser.add_argument('--contig_file', required=True,  help='fasta nucleotide file of the contigs')
    parser.add_argument('--database_dir', required=True,  help='dir that contains the annotation databases')
    args = parser.parse_args()
    return args


def prep_output_dir(args):
    output_dir = Path(Path(args.contig_file).parent, "metaerg")
    if os.path.exists(output_dir):
        utils.log('Warning: may overwrite existing output files...')
        if os.path.isfile(output_dir):
            utils.log(f'Expected folder at {output_dir}, found regular file, crash! Delete this fiel first')
            exit(1)
    else:
        os.mkdir(output_dir)
    os.chdir(output_dir)
    return output_dir


def get_tmp_file(tmp_dir):
    return os.path.join()


def main():
    utils.log(f'This is metaerg.py {VERSION}')
    args = parse_arguments()
    prep_output_dir(args)
    # Filter and load contigs
    input_fasta_file = Path(args.contig_file)
    fasta_file = Path(input_fasta_file.name)
    gbk_file = Path(input_fasta_file.stem + ".gbk")
    gff_file = Path(input_fasta_file.stem + ".gff")

    databases.DBDIR = Path(args.database_dir)
    # Feature prediction
    utils.create_filtered_contig_fasta_file(fasta_file_in=input_fasta_file, fasta_file_out=fasta_file)
    contig_dict = SeqIO.to_dict(SeqIO.parse(fasta_file, "fasta"))
    for contig in contig_dict.values():
        contig.annotations['molecule_type'] = 'DNA'
    contig_dict = {}
    with open(gbk_file) as handle:
        for gb_record in SeqIO.parse(handle, "genbank"):
            contig_dict[gb_record.id] = gb_record
            print(gb_record.id, len(gb_record.features))
    features.discover_transmembrane_helixes(fasta_file, contig_dict)
    # for prediction in (features.predict_crisprs_with_minced,
    #                    features.predict_trnas_with_aragorn,
    #                    features.predict_non_coding_rna_features_with_infernal,
    #                    features.predict_retrotransposons_with_ltrharvest,
    #                    features.predict_tandem_repeats_with_trf,
    #                    features.predict_remaining_repeats_with_repeatmasker,
    #                    features.predict_coding_sequences_with_prodigal,
    #                    features.create_ids,
    #                    databases.load_descriptions_and_taxonomy,
    #                    features.annotate_features_by_homology_diamond,
    #                    features.annotate_features_by_homology_blastn,
    #                    features.annotate_features_by_homology_cdd,
    #                    features.annotate_features_by_homology_antismash
    #                    features.discover_transmembrane_helixes
    #                    ):
    #     prediction(fasta_file, contig_dict)
    #     SeqIO.write(contig_dict.values(), gbk_file, "genbank")
        #except:
         #  utils.log(f'Error running {prediction}; moving to next step.')
    # contig_dict = {}
    # with open(gbk_file) as handle:
    #     for gb_record in SeqIO.parse(handle, "genbank"):
    #         contig_dict[gb_record.id] = gb_record
    #         print(gb_record.id, len(gb_record.features))
    #
    # #SeqIO.write(contig_dict.values(), gbk_file, "genbank")
    # with open(gff_file, "w") as gff_handle:
    #     GFF.write(contig_dict.values(), gff_handle)

    utils.log(f'Done. Thank you for using metaerg.py {VERSION}')


if __name__ == "__main__":
    main()
