import argparse
import os
import gzip
import re
import warnings
import time
from pathlib import Path
from ftplib import FTP

from metaerg import utils
import shutil
from tqdm import tqdm

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio import BiopythonParserWarning

import ncbi.datasets

VERSION = "2.0.17"
RELEVANT_RNA_GENES = 'rRNA tRNA RNase_P_RNA SRP_RNA riboswitch snoRNA ncRNA tmRNA antisense_RNA binding_site ' \
                     'hammerhead_ribozyme scRNA mobile_element mobile_genetic_element misc_RNA autocatalytically_spliced_intron'.split()
IGNORED_FEATURES = 'gene pseudogene exon direct_repeat region sequence_feature pseudogenic_tRNA pseudogenic_rRNA ' \
                   'repeat_region ribosome_entry_site minus_10_signal minus_35_signal protein_binding_site regulatory' \
                   ' transcript mRNA gap assembly_gap source misc_feature variation misc_structure transcript' \
                   ' 3\'UTR 5\'UTR intron signal_peptide_region_of_CDS sequence_alteration'.split()
EUK_ROOT_TAXA = 'Amoebozoa Ancyromonadida Apusozoa Breviatea CRuMs Cryptophyceae Discoba Glaucocystophyceae ' \
                'Haptista Hemimastigophora Malawimonadida Metamonada Rhodelphea Rhodophyta Sar Aphelida ' \
                'Choanoflagellata Filasterea Fungi Ichthyosporea Rotosphaerida'.split()
DB_DESCR_FILENAME = 'db_descriptions.txt'
DB_TAXON_FILENAME = 'db_taxonomy.txt'
CDD_INDEX_FILENAME = "cddid.tbl"
DBDIR = Path('metaerg')
DESCRIPTIONS = {'p': {}, 'e': {}, 'v': {}}
TAXONOMY  = {'p': {}, 'e': {}, 'v': {}}
CDD = {}
DESCRIPTIONS_CACHE = {'p': set(), 'e': set(), 'v': set()}
TAXONOMY_CACHE = {'p': set(), 'e': set(), 'v': set()}
CDD_CACHE = set()
CANTHYD_TRUSTED_CUTOFFS = {}
CANTHYD_DESCR = {'AlkB': 'alkane hydrolase',
                 'AlmA_GroupI':	'flavin-binding alkane monooxygenase',
                 'AlmA_GroupIII': 'flavin-binding alkane monooxygenase',
                 'CYP153': 'alkane oxidizing cytochrome P450',
                 'LadA_alpha': 'long-chain alkane hydrolase',
                 'LadA_beta': 'long-chain alkane hydrolase',
                 'LadB': 'long-chain alkane hydrolase',
                 'pBmoA': 'membrane-bound alkane monooxygenase subunit A',
                 'pBmoB': 'membrane-bound alkane monooxygenase subunit B',
                 'pBmoC': 'membrane-bound alkane monooxygenase subunit C',
                 'PrmA': 'propane 2-monooxygenase large subunit',
                 'PrmC': 'propane 2-monooxygenase small subunit',
                 'sBmoX': 'soluble alkane monooxygenase subunit A',
                 'sBmoY': 'soluble alkane monooxygenase subunit B',
                 'DmpO': 'phenol/toluene 2-monooxygenase (NADH dependent)',
                 'DszC': 'dibenzothiophene monooxygenase',
                 'MAH_alpha': 'benzene/toluene/naphtalene dioxygenase subunit alpha',
                 'MAH_beta': 'benzene/toluene/naphtalene dioxygenase subunit beta',
                 'NdoB': 'benzene/toluene/naphtalene dioxygenase subunit alpha',
                 'non_NdoB_type': 'similar to benzene/toluene/naphtalene dioxygenase subunit alpha',
                 'NdoC': 'benzene/toluene/naphtalene dioxygenase subunit beta',
                 'TmoA_BmoA': 'toluene monooxygenase subunit A',
                 'TmoB_BmoB': 'toluene monooxygenase subunit B',
                 'TmoE': 'toluene monooxygenase system protein E',
                 'TomA1': 'phenol/toluene monooxygenase/hydroxylase (NADH dependent)',
                 'TomA3': 'phenol/toluene monooxygenase/hydroxylase (NADH dependent)',
                 'TomA4': 'phenol/toluene monooxygenase/hydroxylase (NADH dependent)',
                 'ahyA': 'molybdopterin-family alkane C2 methylene hydroxylase',
                 'AssA': 'alkylsuccinate synthase',
                 'AbcA_1': 'benzene carboxylase',
                 'BssA': 'benzylsuccinate synthase',
                 'CmdA': 'molybdopterin-family ethylbenzene dehydrogenase subunit alpha',
                 'EbdA': 'molybdopterin-family ethylbenzene dehydrogenase subunit alpha',
                 'K27540': 'naphtalene carboxylase',
                 'NmsA': 'naphtylmethyl succinate synthase'}

warnings.simplefilter('ignore', BiopythonParserWarning)


def does_db_appear_valid():
    return Path(DBDIR, DB_DESCR_FILENAME).exists() and Path(DBDIR, DB_TAXON_FILENAME).exists() \
        and Path(DBDIR, CDD_INDEX_FILENAME).exists()

def load_descriptions_taxonomy_cdd():
    # load descriptions
    with open(Path(DBDIR, DB_DESCR_FILENAME)) as descr_handle:
        for line in descr_handle:
            words = line.split('\t')
            DESCRIPTIONS[words[0]][int(words[1])] = words[2].strip()
    utils.log(f'Parsed ({len(DESCRIPTIONS["p"])}, {len(DESCRIPTIONS["e"])}, {len(DESCRIPTIONS["v"])}) gene descriptions '
              f'from db for (prokaryotes, eukaryotes and viruses) respectively. ')
    # load taxonomy
    with open(Path(DBDIR, DB_TAXON_FILENAME)) as taxon_handle:
        for line in taxon_handle:
            words = line.split('\t')
            TAXONOMY[words[0]][int(words[1])] = words[2].strip()
    utils.log(f'Parsed ({len(TAXONOMY["p"])}, {len(TAXONOMY["e"])}, {len(TAXONOMY["v"])}) taxa from db for (prokaryotes,'
          f'eukaryotes and viruses) respectively.')
    # load cdd
    cdd_descriptions_file = Path(DBDIR, CDD_INDEX_FILENAME)
    with open(cdd_descriptions_file) as cdd_descr_handle:
        for line in cdd_descr_handle:
            words = line.split("\t")
            # id, id_name, gene_name, descr, length
            CDD[int(words[0])] = (words[1], words[2], words[3], int(words[4]))
    # load canthyd
    canthyd_file =  Path(DBDIR, 'canthyd', 'CANT-HYD.hmm')
    if canthyd_file.exists() and canthyd_file.stat().st_size:
        current_name = None
        with open(canthyd_file) as handle:
            for line in handle:
                if line.startswith('NAME'):
                    current_name = line.split()[1]
                elif line.startswith('TC'):
                    CANTHYD_TRUSTED_CUTOFFS[current_name] = int(line.split()[1])


def decipher_database_id(id, add_to_cache=False):
    words = id.split('~') # org_acc gene_acc [pev] gene# decr# taxon#
    taxon = TAXONOMY[words[2]][int(words[5])]
    descr = DESCRIPTIONS[words[2]][int(words[4])]
    length = int(words[6])
    pos = int(words[3])

    if add_to_cache:
        DESCRIPTIONS_CACHE[words[2]].add(int(words[4]))
        TAXONOMY_CACHE[words[2]].add(int(words[5]))

    return {'taxon': taxon.replace('~', '~ '),
            'descr': descr,
            'length': length,
            'gene_number': pos,
            'descr_id': (words[2], int(words[4]))}


def save_cache(dir):
    with open(Path(dir, DB_DESCR_FILENAME), 'w') as file_handle:
        for kingdom in DESCRIPTIONS_CACHE:
            for db_id in DESCRIPTIONS_CACHE[kingdom]:
                file_handle.write(f'{kingdom}\t{db_id}\t{DESCRIPTIONS[kingdom][db_id]}\n')
    with open(Path(dir, DB_TAXON_FILENAME), 'w') as file_handle:
        for kingdom in TAXONOMY_CACHE:
            for db_id in TAXONOMY_CACHE[kingdom]:
                file_handle.write(f'{kingdom}\t{db_id}\t{TAXONOMY[kingdom][db_id]}\n')
    with open(Path(dir, CDD_INDEX_FILENAME), 'w') as file_handle:
        for cdd_id in CDD_CACHE:
            cdd_item = CDD[cdd_id]
            file_handle.write(f'{cdd_id}\t{cdd_item[0]}\t{cdd_item[1]}\t{cdd_item[2]}\t{cdd_item[3]}\n')


def parse_arguments():
    parser = argparse.ArgumentParser(description='metaerg.py. (C) Marc Strous, Xiaoli Dong 2019, 2021')
    parser.add_argument('--target_dir', required=True,  help='where to create the database')
    parser.add_argument('--gtdbtk_dir', required=True,  help='where to create the database')
    parser.add_argument('--tasks', required=False,  default='FPVERBCS', help='F = create folders, '
                        'P = download prokaryotes, V = download viruses, E = download eukaryotes, '
                        'B = build P,V,E blast databases, R = install rfam, C = install cdd, S = install specialized')

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
    for dir in (settings["ncbi_cache"], settings["ncbi_cache_pro"], settings["ncbi_cache_euk"],
                settings["ncbi_cache_vir"]):
        P = Path(dir)
        P.mkdir(exist_ok=True)


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
        file_handle.write_html(f'{kingdom}\t{descr_id}\t{description}\n')
    return descr_id


def extract_proteins_and_rna_prok(settings, t):
    contig_dict = dict()
    descr_dict = settings["description_dict"]
    try:
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
                    try:
                        contig = contig_dict[words[0]]
                    except KeyError:
                        continue
                    feature_seq = utils.pad_seq(feature.extract(contig))
                    prot_seq = feature_seq.translate(table=translation_table)
                    prot_seq = prot_seq[:-1]
                    if '*' in prot_seq:
                        continue
                    product = utils.get_feature_qualifier(feature, "product")
                    if not product:
                        continue
                    prot_seq.description = product
                    descr_id = update_descriptions_and_get_id(prot_seq.description, descr_dict, description_handle, 'p')
                    prot_seq.id = f'{t["accession"]}~{feature.qualifiers["id"]}~p~{gene_counter}~{descr_id}~{t["id"]}~{len(prot_seq.seq)}'
                    SeqIO.write(prot_seq, prot_fasta_out_handle, "fasta")
                elif words[2] in RELEVANT_RNA_GENES:
                    gene_counter += 1
                    feature = utils.gff_words_to_seqfeature(words)
                    try:
                        contig = contig_dict[words[0]]
                    except KeyError:
                        continue
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
    except EOFError:
        utils.log(f'Incomplete file detected for {t["cached_gff_file"]}. Deleting...')
        t["cached_gff_file"].unlink()


def download_genomes_with_datasets(settings, taxon_list):
    accessions = [tt['accession'] for tt in taxon_list]
    if not len(accessions):
        return
    for attempt in range(3):
        time.sleep(10)
        try:
            genome_api = ncbi.datasets.GenomeApi(ncbi.datasets.ApiClient())
            api_response = genome_api.download_assembly_package(
                accessions=accessions,
                exclude_sequence=True,
                include_annotation_type=['GENOME_GFF'],
                _preload_content=False,
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
            break
        except ncbi.datasets.openapi.exceptions.ApiException:
            utils.log(f'NCBI has flagged us for too many requests, taking a 20 s breath.')
            time.sleep(10)

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
    taxa_count = 0
    # determine line count for tqdm
    with open(settings["gtdbtk_taxonomy_file"]) as taxonomy_handle:
        for line in taxonomy_handle:
            taxa_count += 1
    taxon_list = []
    with open(settings["gtdbtk_taxonomy_file"]) as taxonomy_handle:
        # determine line count for tqdm
        ftp = FTP('ftp.ncbi.nlm.nih.gov')
        ftp.login()
        success_count = 0
        genomes_in_cash_count = 0
        count = 0
        for line in taxonomy_handle:
            count += 1
            words = line.split()
            accession = words[0][3:]
            taxonomy = re.sub("\w__", "", words[1]).split(';')
            seq_file = Path(settings["gtdbtk_genome_sequences_root"], accession[0:3], accession[4:7],
                            accession[7:10], accession[10:13], f'{accession}_genomic.fna.gz')
            future_gff_file = Path(settings["ncbi_cache_pro"], f'{accession}.gff.gz')
            taxon = {'id': len(taxon_list), 'accession': accession, 'taxonomy': taxonomy,
                    'cached_gff_file': future_gff_file, 'gtdb_seq_file': seq_file, 'in_local_cache': False}
            taxon_list.append(taxon)
            download_status = " "
            if os.path.exists(future_gff_file):
                taxon['in_local_cache'] = True
                download_status = "+"
            else:
                continue
                # print(f'now retrieving {taxon}')
                # try:
                #     download_status = "*"
                #     ftp.cwd('/genomes/all/')
                #     for acc_part in (accession[0:3], accession[4:7], accession[7:10], accession[10:13]):
                #         ftp.cwd(acc_part)
                #     ftp_dir_list = []
                #     ftp.dir('.', ftp_dir_list.append)
                #     ftp.cwd(ftp_dir_list[0].split()[-1])
                #     ftp_dir_list = []
                #     ftp.dir('.', ftp_dir_list.append)
                #     target = None
                #     for l in ftp_dir_list:
                #         filename = l.split()[-1]
                #         if filename.endswith('_genomic.gff.gz'):
                #             target = filename
                #             break
                #     if target:
                #         success_count += 1
                #         with open(future_gff_file, "wb") as local_handle:
                #             ftp.retrbinary("RETR " + target, local_handle.write)
                #         taxon['in_local_cache'] = True
                #     else:
                #         continue
                # except EOFError:
                #     utils.log('FTP Error at NCBI - resetting connection')
                #     ftp = FTP('ftp.ncbi.nlm.nih.gov')
                #     ftp.login()
                #     continue
                # except BrokenPipeError:
                #     utils.log('FTP Error at NCBI - resetting connection')
                #     ftp = FTP('ftp.ncbi.nlm.nih.gov')
                #     ftp.login()
                #     continue
            genomes_in_cash_count += 1
            utils.log(f'({count}/{taxa_count}) {download_status} {future_gff_file.name} {taxonomy}')
            extract_proteins_and_rna_prok(settings, taxon)
    utils.log(f'downloaded {success_count} new genomes. Total genomes in cache: {genomes_in_cash_count}.')
    # (4) save taxonomy (includes status)
    utils.log(f'Writing taxonomy file.')
    with open(taxon_db, "w") as taxon_handle:
        for t in taxon_list:
            taxon_handle.write(f'p\t{t["id"]}\t{"~".join(t["taxonomy"])}\t{t["in_local_cache"]}\t{t["accession"]}\n')
    # (5) clean up
    #unlink_files((os.path.join(settings["db_dir"], 'ncbi_genomes.zip'), os.path.join(settings["db_dir"],
    #                                                                                 '../../README.md')))


def prep_prokaryote_database_o(settings):
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
                    seq_record = SeqRecord(gb_record.seq)
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
            if file.endswith(".gbff"):
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
        cdd_index = Path(cdd_dir, CDD_INDEX_FILENAME)
        utils.run_external(f'gunzip {cdd_index}.gz')
        utils.run_external(f'cp {cdd_index} {settings["db_dir"]}')
        utils.run_external(f'wget -P {cdd_dir} https://ftp.ncbi.nih.gov/pub/mmdb/cdd/cdd.tar.gz')
        utils.run_external(f'tar -xf cdd.tar.gz')
        utils.run_external(f'makeprofiledb -title CDD.v.3.12 -in Cdd.pn -out Cdd'
                           f' -threshold 9.82 -scale 100.0 -dbtype rps -index true')
    else:
        utils.log('Keeping previously installed conserved domain database...')


def prep_other(settings):
    canthyd_dir = settings["canthyd"]
    if not os.path.exists(canthyd_dir):
        utils.log(f'Installing the conserved domain database to {canthyd_dir}...')
        os.mkdir(canthyd_dir)
        os.chdir(canthyd_dir)
        utils.run_external(f'wget https://github.com/dgittins/CANT-HYD-HydrocarbonBiodegradation/raw/main/HMMs/concatenated%20HMMs/CANT-HYD.hmm')
        utils.run_external(f'hmmpress -f CANT-HYD.hmm')


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
                'canthyd': os.path.join(args.target_dir, "canthyd"),
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
    if 'S' in args.tasks:
        utils.log('Now downloading and installing specialized databases...')
        prep_other(settings)
    utils.log("Done.")

if __name__ == '__main__':
    main()
