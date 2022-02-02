import argparse
import os
import gzip
import re
import warnings
from pathlib import Path
from collections import deque, Counter

import utils
import subprocess
import shutil
from tqdm import tqdm

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio import BiopythonParserWarning

import ncbi.datasets

VERSION = 0.0
RELEVANT_RNA_GENES = 'rRNA tRNA RNase_P_RNA SRP_RNA riboswitch snoRNA ncRNA tmRNA antisense_RNA binding_site ' \
                     'hammerhead_ribozyme scRNA mobile_element mobile_genetic_element misc_RNA'.split()
IGNORED_FEATURES = 'gene pseudogene exon direct_repeat region sequence_feature pseudogenic_tRNA pseudogenic_rRNA ' \
                   'repeat_region ribosome_entry_site minus_10_signal minus_35_signal protein_binding_site regulatory' \
                   ' transcript mRNA gap assembly_gap source misc_feature variation misc_structure transcript' \
                   ' 3\'UTR 5\'UTR intron'.split()
EUK_ROOT_TAXA = 'Amoebozoa Ancyromonadida Apusozoa Breviatea CRuMs Cryptophyceae Discoba Glaucocystophyceae ' \
                'Haptista Hemimastigophora Malawimonadida Metamonada Rhodelphea Rhodophyta Sar Aphelida ' \
                'Choanoflagellata Filasterea Fungi Ichthyosporea Rotosphaerida'.split()
DBDIR = Path('.')
DESCRIPTIONS = {'p': [], 'e': [], 'v': []}
TAXONOMY  = {'p': [], 'e': [], 'v': []}


def load_descriptions_and_taxonomy(x=0, y=0):
    with open(Path(DBDIR, 'db_descriptions.txt')) as descr_handle:
        for line in descr_handle:
            words = line.split('\t')
            DESCRIPTIONS[words[0]].append(words[2].strip())
    utils.log(f'Parsed ({len(DESCRIPTIONS["p"])}, {len(DESCRIPTIONS["e"])}, {len(DESCRIPTIONS["v"])}) gene descriptions '
              f'from db for (prokaryotes, eukaryotes and viruses) respectively. ')
    # load taxonomy
    with open(Path(DBDIR, 'db_taxonomy.txt')) as taxon_handle:
        for line in taxon_handle:
            words = line.split('\t')
            TAXONOMY[words[0]].append(words[2].strip())
    utils.log(f'Parsed ({len(TAXONOMY["p"])}, {len(TAXONOMY["e"])}, {len(TAXONOMY["v"])}) taxa from db for (prokaryotes,'
          f'eukaryotes and viruses) respectively.')


warnings.simplefilter('ignore', BiopythonParserWarning)


def parse_arguments():
    parser = argparse.ArgumentParser(description='metaerg.py. (C) Marc Strous, Xiaoli Dong 2019, 2021')
    parser.add_argument('--target_dir', required=True,  help='where to create the database')
    parser.add_argument('--gtdbtk_dir', required=True,  help='where to create the database')
    parser.add_argument('--tasks', required=False,  default='FPVERB', help='F = folders, P = prokaryotes, V = viruses, '
                                                                           'E = eukaryotes, R = rfam, B = build_db')

    args = parser.parse_args()
    return args


def prep_database_folder(settings):
    if os.path.exists(settings["db_dir"]):
        utils.log('Warning: may overwrite existing database files...')
        if os.path.isfile(settings["db_dir"]):
            utils.log(f'Expected folder at {settings["db_dir"]}, found regular file, terminating now.')
            exit(1)
    else:
        os.mkdir(settings["db_dir"])
        os.mkdir(settings["ncbi_cache"])
        os.mkdir(settings["ncbi_cache_pro"])
        os.mkdir(settings["ncbi_cache_euk"])
        os.mkdir(settings["ncbi_cache_vir"])


def unlink_files(files):
    for f in files:
        try:
            os.unlink(f)
        except FileNotFoundError:
            pass


def rank_assembly_level(level):
    if 'Scaffold' == level:
        return 2
    if 'Chromosome' == level:
        return 3
    if 'Complete Genome' == level:
        return 4
    return 0


def parse_euk_feature_description(feature):
    if 'product' in feature.qualifiers:
        return feature.qualifiers['product'][0]
    elif 'gene' in feature.qualifiers:
        return feature.qualifiers['gene'][0]
    elif 'note' in feature.qualifiers:
        return feature.qualifiers['note'][0]
    elif 'mobile_element_type' in feature.qualifiers:
        return feature.qualifiers['mobile_element_type'][0]
    elif 'ncRNA_class' in feature.qualifiers:
        return feature.qualifiers['ncRNA_class'][0]
    else:
        return f'Unknown {feature.type}'


def update_descriptions_and_get_id(description, dictionary, file_handle, kingdom):
    try:
        descr_id = dictionary[description]
    except KeyError:
        descr_id = len(dictionary)
        dictionary[description] = descr_id
        file_handle.write(f'{kingdom}\t{descr_id}\t{description}\n')
    return descr_id


def extract_proteins_and_rna_prok(settings, t):
    contig_dict = dict()
    descr_dict = settings["description_dict"]
    with gzip.open(t["gtdb_seq_file"], "rt") as handle:
        for seq in SeqIO.parse(handle, "fasta"):
            contig_dict[seq.id] = seq
    with gzip.open(t["cached_gff_file"], "rt") as gff_handle, \
            open(os.path.join(settings['ncbi_cache_pro'], settings["faa_db_name"]), 'a') as prot_fasta_out_handle, \
            open(os.path.join(settings['ncbi_cache_pro'], settings["fna_db_name"]), 'a') as rna_fasta_out_handle, \
            open(os.path.join(settings['ncbi_cache_pro'], settings["descr_db_name"]), 'a') as description_handle:
        gene_counter = 0
        for line in gff_handle:
            if line.startswith("#"):
                continue
            if 'similar to' in line or 'identical to' in line:
                continue
            words = line.split('\t')
            if len(words) < 9:
                continue
            if 'CDS' == words[2]:
                gene_counter += 1
                feature = utils.gff_words_to_seqfeature(words)
                translation_table = int(feature.qualifiers['transl_table'])
                contig = contig_dict[words[0]]
                feature_seq = utils.pad_seq(feature.extract(contig))
                prot_seq = feature_seq.translate(table=translation_table)
                prot_seq = prot_seq[:-1]
                if '*' in prot_seq:
                    continue

                prot_seq.description = feature.qualifiers["product"]
                descr_id = update_descriptions_and_get_id(prot_seq.description, descr_dict, description_handle, 'p')
                prot_seq.id = f'{t["accession"]}~{feature.qualifiers["id"]}~p~{gene_counter}~{descr_id}~{t["id"]}~{len(prot_seq.seq)}'
                SeqIO.write(prot_seq, prot_fasta_out_handle, "fasta")
            elif words[2] in RELEVANT_RNA_GENES:
                gene_counter += 1
                feature = utils.gff_words_to_seqfeature(words)
                contig = contig_dict[words[0]]
                feature_seq = feature.extract(contig)
                try:
                    feature_seq.description = feature.qualifiers["product"]
                except KeyError:
                    if "note" in feature.qualifiers.keys():
                        feature_seq.description = feature.qualifiers["note"]
                    else:
                        continue
                descr_id = update_descriptions_and_get_id(feature_seq.description, descr_dict, description_handle, 'p')
                feature_seq.id = f'{t["accession"]}~{feature.qualifiers["id"]}~p~{gene_counter}~{descr_id}~{t["id"]}~{len(feature_seq.seq)}'
                SeqIO.write(feature_seq, rna_fasta_out_handle, "fasta")
            elif words[2] in IGNORED_FEATURES:
                pass
            else:
                utils.log(f'  Warning: unknown feature type in line "{line}"')


def download_genomes_with_datasets(settings, taxon_list):
    accessions = [tt['accession'] for tt in taxon_list]
    if not len(accessions):
        return
    for attempt in range(3):
        try:
            genome_api = ncbi.datasets.GenomeApi(ncbi.datasets.ApiClient())
            api_response = genome_api.download_assembly_package(
                accessions=accessions,
                exclude_sequence=True,
                include_annotation_type=['GENOME_GFF'],
                _preload_content=False
            )
            zipfilename = os.path.join(settings["db_dir"], 'ncbi_genomes.zip')
            with open(zipfilename, 'wb') as f:
                f.write(api_response.data)
            if not os.path.exists(zipfilename):
                return
            try:
                utils.run_external(f'unzip -qq -o -d {settings["db_dir"]} {zipfilename}', log_cmd=False)
            except:
                return
            unlink_files((zipfilename, ))
            for t in taxon_list:
                gff_file = os.path.join(settings["db_dir"], 'ncbi_dataset', 'data', t['accession'], 'genomic.gff')
                # surprisingly, refseq entries may have disappeared from NCBI, so not all gff files will
                # actually be downloaded
                if os.path.exists(gff_file):
                    with open(gff_file, 'rb') as f_in, gzip.open(t['cached_gff_file'], 'wb') as f_out:
                        f_out.writelines(f_in)
                    t['retrieved'] = True
                    extract_proteins_and_rna_prok(settings, t)

            return
        except ncbi.datasets.openapi.exceptions.ServiceException:
            utils.log(f'NCBI error for accessions {accessions}, aborting.')
        #except:
        #    utils.log("Error, perhaps while decompressing zipfile, retrying...")
    utils.log("Three failed retrieval attempts - giving up for this batch.")


def prep_prokaryote_database(settings):
    # (1) Delete some previous files
    fasta_protein_db = os.path.join(settings['ncbi_cache_pro'], settings["faa_db_name"])
    fasta_nt_db = os.path.join(settings['ncbi_cache_pro'], settings["fna_db_name"])
    taxon_db = os.path.join(settings['ncbi_cache_pro'], settings["taxon_db_name"])
    descr_db = os.path.join(settings['ncbi_cache_pro'], settings["descr_db_name"])
    unlink_files((fasta_protein_db, fasta_nt_db, descr_db, taxon_db))

    # (2) Load and prep list of target taxa from GTDB
    taxon_list = []
    with open(settings["gtdbtk_taxonomy_file"]) as taxonomy_handle:
        for line in taxonomy_handle:
            words = line.split()
            accession = words[0][3:]
            taxonomy = re.sub("\w__", "", words[1]).split(';')
            is_gbk = line.startswith('GB_')
            retrieved = False
            seq_file = os.path.join(settings["gtdbtk_genome_sequences_root"], accession[0:3], accession[4:7],
                                    accession[7:10], accession[10:13], f'{accession}_genomic.fna.gz')
            future_gff_file = os.path.join(settings["ncbi_cache_pro"], f'{accession}.gff.gz')
            taxon_list.append({'id': len(taxon_list), 'accession': accession, 'taxonomy': taxonomy,
                               'retrieved': retrieved, 'is_gbk': is_gbk, 'cached_gff_file': future_gff_file,
                               'gtdb_seq_file': seq_file})

    # (3) Download corresponding gff files from NCBI and extracting seqs
    utils.log(f'Acquiring {len(taxon_list)} genomes from GTDB template.')
    current_targets_for_download = []
    for t in tqdm(taxon_list):
        if t["is_gbk"]:
            continue
        if os.path.exists(t["cached_gff_file"]):
            t['retrieved'] = True
            extract_proteins_and_rna_prok(settings, t)
            continue

        current_targets_for_download.append(t)
        if len(current_targets_for_download) >= 100:
            download_genomes_with_datasets(settings, current_targets_for_download)
            current_targets_for_download.clear()
    download_genomes_with_datasets(settings, current_targets_for_download)

    # (4) save taxonomy (includes status)
    utils.log(f'Writing taxonomy file.')
    with open(taxon_db, "w") as taxon_handle:
        for t in taxon_list:
            taxon_handle.write(f'p\t{t["id"]}\t{"~".join(t["taxonomy"])}\t{t["retrieved"]}\t{t["accession"]}\n')
    # (5) clean up
    unlink_files((os.path.join(settings["db_dir"], 'ncbi_genomes.zip'), os.path.join(settings["db_dir"],
                                                                                     '../../README.md')))


def prep_viral_database(settings):
    fasta_protein_db = os.path.join(settings['ncbi_cache_vir'], settings["faa_db_name"])
    fasta_nt_db = os.path.join(settings['ncbi_cache_vir'], settings["fna_db_name"])
    descr_db = os.path.join(settings['ncbi_cache_vir'], settings["descr_db_name"])
    taxon_db = os.path.join(settings['ncbi_cache_vir'], settings["taxon_db_name"])
    unlink_files((fasta_protein_db, fasta_nt_db, descr_db, taxon_db))

    descr_dict = dict()
    taxon_dict = dict()

    pattern = re.compile(r'\s*\[(.+?)\]$')
    with open(fasta_protein_db, 'w') as prot_fasta_out_handle, \
            open(descr_db, 'w') as desrc_handle,\
            open(taxon_db, 'w') as taxon_handle:
        for i in range(1,5):
            f = os.path.join(settings["ncbi_cache_vir"], f'viral.{i}.protein.gpff.gz')
            if not os.path.exists(f):
                url = 'https://ftp.ncbi.nlm.nih.gov/refseq/release/viral/viral'
                utils.run_external(f'wget -P {settings["ncbi_cache_vir"]} {url}.{i}.protein.gpff.gz')
            gene_count = 0
            with gzip.open(f, 'rt') as file_handle:
                for gb_record in SeqIO.parse(file_handle, "genbank"):
                    if gb_record.annotations['molecule_type'] != 'protein':
                        utils.log("warning: skipping non-coding gene in viral refseq")
                        continue
                    match = pattern.search(gb_record.description)
                    if match:
                        gb_record.annotations['taxonomy'].append(match.group(1))
                        gb_record.description = gb_record.description[0:match.start()]
                    else:
                        utils.log(f'Warning, Failed to parse viral species from {gb_record.description}')

                    taxon_str = "~".join(gb_record.annotations['taxonomy'])

                    try:
                        taxon_id = taxon_dict[taxon_str]
                    except KeyError:
                        taxon_id = len(taxon_dict)
                        taxon_dict[taxon_str] = taxon_id
                        taxon_handle.write(f'v\t{taxon_id}\t{taxon_str}\n')
                    descr_id = update_descriptions_and_get_id(gb_record.description, descr_dict, desrc_handle, 'v')
                    seq_record = SeqRecord(Seq(gb_record.seq))
                    seq_record.id = f'{taxon_id}~{gb_record.id}~v~{gene_count}~{descr_id}~{taxon_id}~{len(seq_record.seq)}'
                    seq_record.description = gb_record.description
                    SeqIO.write(seq_record, prot_fasta_out_handle, "fasta")
                    gene_count += 1


def prep_eukaryote_database(settings):
    fasta_protein_db = os.path.join(settings['ncbi_cache_euk'], settings["faa_db_name"])
    fasta_nt_db = os.path.join(settings['ncbi_cache_euk'], settings["fna_db_name"])
    descr_db = os.path.join(settings['ncbi_cache_euk'], settings["descr_db_name"])
    taxon_db = os.path.join(settings['ncbi_cache_euk'], settings["taxon_db_name"])
    unlink_files((fasta_protein_db, fasta_nt_db, descr_db, taxon_db))

    # (1) probe the NCBI for available eukaryotic assemblies
    api_instance = ncbi.datasets.GenomeApi(ncbi.datasets.ApiClient())
    targets = dict()
    taxa_done = dict()
    for tax_name in EUK_ROOT_TAXA:
        genome_summary = api_instance.assembly_descriptors_by_taxon(taxon=tax_name, page_size=1000)
        count = 0
        if genome_summary.assemblies:
            for assembly in map(lambda d: d.assembly, genome_summary.assemblies):
                if assembly.assembly_accession[:3] == 'GCF':
                    if assembly.org.tax_id in taxa_done.keys():
                        new_level = rank_assembly_level(assembly.assembly_level)
                        if new_level > taxa_done[assembly.org.tax_id]:
                            #print(f'{new_level} > {taxa_done[assembly.org.tax_id]}')
                            taxa_done[assembly.org.tax_id] = new_level
                            targets[assembly.org.tax_id] = assembly.assembly_accession.strip()
                        else:
                            continue
                    else:
                        count += 1
                        taxa_done[assembly.org.tax_id] = rank_assembly_level(assembly.assembly_level)
                        targets[assembly.org.tax_id] = assembly.assembly_accession.strip()
        try:
            total_count = int(genome_summary.total_count)
        except TypeError:
            total_count = 0
        utils.log(f"({count}/{total_count}) taxon-unique assemblies for {tax_name}")
    total_target_count = len(targets)
    utils.log(f'Found  {total_target_count} targets for download in total.')
    # (2) Download new assemblies not yet in local cache
    targets = [accession for accession in targets.values() if
               not os.path.exists(os.path.join(settings['ncbi_cache_euk'], f'{accession}.gbk.gz'))]
    utils.log(f'Of those, {total_target_count - len(targets)} already in local cache. ')

    api_response = api_instance.download_assembly_package(targets, exclude_sequence = False,
                                                          include_annotation_type = ['GENOME_GBFF'],
                                                          _preload_content = False
    )
    zipfilename = os.path.join(settings["db_dir"], 'ncbi_genomes.zip')
    unlink_files((zipfilename,))
    with open(zipfilename, 'wb') as f:
        f.write(api_response.data)
    if not os.path.exists(zipfilename):
        return
    utils.run_external(f'unzip -qq -o -d {settings["db_dir"]} {zipfilename}', log_cmd=False)
    # (3) move gbk file into local cache gzipped
    count=0
    for accession in targets:
        count += 1
        src_dir = os.path.join(settings["db_dir"], 'ncbi_dataset', 'data', accession)
        dest_file = os.path.join(settings['ncbi_cache_euk'], f'{accession}.gbk.gz')
        if os.path.exists(dest_file):
            continue
        utils.log(f'({count}/{len(targets)}) Now extracting {accession} into ncbi-cache as {dest_file}...')
        gbff_count = 0
        for file in os.listdir(src_dir):
            if file.endswith('.gbff'):
                with open(os.path.join(src_dir, file), 'rb') as f_in, gzip.open(dest_file, 'wb') as f_out:
                    f_out.writelines(f_in)
                gbff_count += 1
        if gbff_count > 1:
            utils.log(f'WARNING: {src_dir} has {gbff_count} gbff files')
        elif gbff_count == 0:
            utils.log(f'WARNING: {src_dir} has no gbff files')
    unlink_files((zipfilename, os.path.join(settings["db_dir"], '../../README.md')))
    # (4) build blast databases
    descr_dict = dict()
    taxon_dict = dict()
    with open(fasta_protein_db, 'w') as prot_fasta_out_handle,open(fasta_nt_db, 'w')  as rna_fasta_out_handle, \
            open(descr_db, 'w') as desrc_handle, open(taxon_db, 'w') as taxon_handle:
        targets = os.listdir(settings['ncbi_cache_euk'])
        genome_count = 0
        for file in targets:
            genome_count += 1
            if not file.endswith('.gbk.gz'):
                continue
            utils.log(f'({genome_count}/{len(targets)}) Now extracting "{file}"')
            file_path = os.path.join(settings['ncbi_cache_euk'], file)
            feature_success_counter = 0
            feature_total_counter = 0
            without_translation_counter = 0
            with gzip.open(file_path, "rt") as handle:
                first = True
                for gb_record in SeqIO.parse(handle, "genbank"):
                    if first:
                        taxon = gb_record.annotations['taxonomy']
                        taxon.append(gb_record.annotations['organism'])
                        taxon_str = "~".join(taxon)
                        try:
                            taxon_id = taxon_dict[taxon_str]
                        except KeyError:
                            taxon_id = len(taxon_dict)
                            taxon_dict[taxon_str] = taxon_id
                            taxon_handle.write(f'e\t{taxon_id}\t{taxon_str}\n')
                        first = False
                    feature_total_counter += len(gb_record.features)
                    for feature in gb_record.features:
                        if feature.type == 'CDS':
                            if not 'translation' in feature.qualifiers:
                                without_translation_counter += 1
                                continue
                            seq_record = SeqRecord(Seq(feature.qualifiers['translation'][0]))
                            seq_record.description = parse_euk_feature_description(feature)
                            if not seq_record.description:
                                utils.log(f"WARNING: No description for feature {feature}")
                                continue
                            descr_id = update_descriptions_and_get_id(seq_record.description, descr_dict, desrc_handle,
                                                                      'e')
                            seq_record.id = f'{gb_record.id}~0~e~{feature_success_counter}~{descr_id}~{taxon_id}~{len(seq_record.seq)}'
                            SeqIO.write(seq_record, prot_fasta_out_handle, "fasta")
                            feature_success_counter += 1
                        elif feature.type in RELEVANT_RNA_GENES:
                            seq_record = SeqRecord(feature.extract(gb_record.seq))
                            seq_record.description = parse_euk_feature_description(feature)
                            if not seq_record.description:
                                utils.log(f"WARNING: No description for feature {feature}")
                                continue
                            descr_id = update_descriptions_and_get_id(seq_record.description, descr_dict, desrc_handle,
                                                                      'e')
                            seq_record.id = f'{gb_record.id}~0~e~{feature_success_counter}~{descr_id}~{taxon_id}~{len(seq_record.seq)}'
                            SeqIO.write(seq_record, rna_fasta_out_handle, "fasta")
                            feature_success_counter += 1
                        elif feature.type in IGNORED_FEATURES:
                            pass
                        else:
                            utils.log(f'  Warning: unknown feature "{feature}"')
            utils.log(f'  ... extracted ({feature_success_counter}/{feature_total_counter}) proteins and RNAs. '
                      f'{without_translation_counter} proteins skipped for lack of translation.')
    unlink_files((os.path.join(settings["db_dir"], 'ncbi_genomes.zip'), os.path.join(settings["db_dir"],
                                                                                     '../../README.md')))


def build_blast_db(settings):
    for filename in [settings['faa_db_name'], settings['fna_db_name'], settings['taxon_db_name'],
                    settings['descr_db_name']]:
        with open(os.path.join(settings['db_dir'], filename), mode="wb") as destination:
            for dir in [settings['ncbi_cache_pro'], settings['ncbi_cache_euk'], settings['ncbi_cache_vir']]:
                src_file = os.path.join(dir, filename)
                if os.path.exists(src_file):
                    with open(os.path.join(dir, filename), mode="rb") as source:
                        shutil.copyfileobj(source, destination)
    fasta_protein_db = os.path.join(settings['db_dir'], settings['faa_db_name'])
    utils.run_external(f'diamond makedb --in {fasta_protein_db} --db {fasta_protein_db}')
    fasta_nt_db = os.path.join(settings['db_dir'], settings['fna_db_name'])
    utils.run_external(f'makeblastdb -in {fasta_nt_db} -dbtype nucl')


def decipher_database_id(id):
    words = id.split('~') # org_acc gene_acc [pev] gene# decr# taxon#
    taxon = TAXONOMY[words[2]][int(words[5])]
    descr = DESCRIPTIONS[words[2]][int(words[4])]
    length = int(words[6])
    pos = int(words[3])
    return {'taxon': taxon,
            'descr': descr,
            'length': length,
            'gene_number': pos,
            'descr_id': (words[2], int(words[4]))}


def prep_rfam(settings):
    rfam_file = settings["rfam"]
    if not os.path.exists(rfam_file):
        utils.log(f'Installing the RFAM database to {rfam_file}...')
        utils.run_external(f'wget -P {settings["db_dir"]} http://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/Rfam.cm.gz')
        utils.run_external(f'gunzip {rfam_file}.gz')
    else:
        utils.log('Keeping previously installed RFAM database...')
    utils.log(f'Running cmpress...')
    if not os.path.exists(os.path.join(settings["db_dir"], "Rfam.cm.i1f")):
        utils.run_external(f'cmpress -F {rfam_file}')
    else:
        utils.log('Skipping cmpress for previously cmpressed RFAM database...')


def prep_cdd(settings):
    cdd_dir = settings["cdd"]
    if not os.path.exists(cdd_dir):
        utils.log(f'Installing the conserved domain database to {cdd_dir}...')
        os.mkdir(cdd_dir)
        os.chdir(cdd_dir)
        utils.run_external(f'wget https://ftp.ncbi.nih.gov/pub/mmdb/cdd/cddid.tbl.gz')
        cdd_index = os.path.join(cdd_dir, 'cddid.tbl')
        utils.run_external(f'gunzip {cdd_index}.gz')
        utils.run_external(f'wget -P {cdd_dir} https://ftp.ncbi.nih.gov/pub/mmdb/cdd/cdd.tar.gz')
        utils.run_external(f'tar -xf cdd.tar.gz')
        utils.run_external(f'makeprofiledb -title CDD.v.3.12 -in Cdd.pn -out Cdd'
                           f' -threshold 9.82 -scale 100.0 -dbtype rps -index true')
    else:
        utils.log('Keeping previously installed conserved domain database...')


def main():
    utils.log(f'This is metaerg.py\'s make database script {VERSION}')
    args = parse_arguments()
    gtdbtk_root = args.gtdbtk_dir
    for f in os.listdir(args.gtdbtk_dir):
        if 'release' in f.lower():
            gtdbtk_root = os.path.join(args.gtdbtk_dir, f)
            break
    settings = {'db_dir': args.target_dir,
                'rfam': os.path.join(args.target_dir, "Rfam.cm"),
                'cdd': os.path.join(args.target_dir, "cdd"),
                'gtdbtk_root': gtdbtk_root,
                'ncbi_cache': os.path.join(args.target_dir, "ncbi-cache"),
                'ncbi_cache_pro': os.path.join(args.target_dir, "ncbi-cache", "pro"),
                'ncbi_cache_euk': os.path.join(args.target_dir, "ncbi-cache", "euk"),
                'ncbi_cache_vir': os.path.join(args.target_dir, "ncbi-cache", "vir"),
                'gtdbtk_taxonomy_file' : os.path.join(gtdbtk_root, "taxonomy", "gtdb_taxonomy.tsv"),
                'gtdbtk_genome_sequences_root': os.path.join(gtdbtk_root, "fastani", "database"),
                'faa_db_name': "db_protein.faa",
                'fna_db_name': "db_rna.fna",
                'descr_db_name': "db_descriptions.txt",
                'taxon_db_name': "db_taxonomy.txt",
                'description_dict': dict()}

    if 'F' in args.tasks:
        utils.log('Creating folders...')
        prep_database_folder(settings)
    if 'P' in args.tasks:
        utils.log('Now adding prokaryote proteins and rna to search databases from GTDBtk data and NCBI annotations...')
        prep_prokaryote_database(settings)
    if 'V' in args.tasks:
        utils.log('Now adding viral proteins to search databases from viral refseq...')
        prep_viral_database(settings)
    if 'E' in args.tasks:
        utils.log('Now adding eukaryote proteins and rna to search databases from NCBI...')
        prep_eukaryote_database(settings)
    if 'B' in args.tasks:
        utils.log('Now building diamond and blast databases...')
        build_blast_db(settings)
    if 'R' in args.tasks:
        utils.log('Now downloading and installing rfam...')
        prep_rfam(settings)
    if 'C' in args.tasks:
        utils.log('Now downloading and installing cdd...')
        prep_cdd(settings)
    utils.log("Done.")

if __name__ == '__main__':
    main()
