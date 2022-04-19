import argparse
import os
from pathlib import Path

from Bio import SeqIO
from BCBio import GFF

from metaerg import databases
from metaerg import predict
from metaerg import utils
from metaerg import visualization

VERSION = "2.0.15"


def parse_arguments():
    parser = argparse.ArgumentParser(description='metaerg.py. (C) Marc Strous, Xiaoli Dong 2019, 2021')
    parser.add_argument('--contig_file', required=True,  help='fasta nucleotide file of the contigs')
    parser.add_argument('--database_dir', required=True,  help='dir that contains the annotation databases')
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
    working_directory = input_fasta_file.parent # eventually: os.getcwd()
    if not input_fasta_file.exists() or input_fasta_file.is_dir():
        utils.log(f'Input file "{input_fasta_file}" is missing or not a valid file. Expecting a nt fasta file.')
        exit(1)
    mag_name = input_fasta_file.stem

    # (3) set output files, note that these file paths are relative to the current working dir
    output_dir, html_dir = create_output_dir(working_directory)
    os.chdir(output_dir)
    gbk_file = predict.spawn_file("gbk", mag_name)
    gff_file = predict.spawn_file("gff", mag_name)
    faa_file = predict.spawn_file("faa", mag_name)
    fna_file = predict.spawn_file("rna.fna", mag_name)

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
                       #predict.predict_remaining_repeats_with_repeatmasker,
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
                       predict.write_databases
                       ):
        prediction(mag_name, contig_dict)
        SeqIO.write(contig_dict.values(), gbk_file, "genbank")
        with open(gff_file, "w") as gff_handle:
             GFF.write(contig_dict.values(), gff_handle)

    #load blast results
    blast_results = {'diamond': predict.spawn_file('diamond', mag_name),
                     'blastn': predict.spawn_file('blastn', mag_name),
                     'cdd': predict.spawn_file('cdd', mag_name)}
    for key in blast_results.keys():
        with utils.TabularBlastParser(blast_results[key]) as handle:
            blast_result_hash = {}
            for br in handle:
                blast_result_hash[br[0]] = br[1]
        blast_results[key] = blast_result_hash
    # write main result files
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
    visualization.html_save_all(mag_name, contig_dict, blast_results)
    utils.log(f'Done. Thank you for using metaerg.py {VERSION}')


if __name__ == "__main__":
    main()
