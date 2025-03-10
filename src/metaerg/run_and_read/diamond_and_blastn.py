import re
import gzip
import shutil
from collections import Counter
from os import chdir
from ftplib import FTP
from pathlib import Path
import zlib

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import ncbi.datasets

from metaerg import context
from metaerg.datatypes import fasta
from metaerg.datatypes import sqlite
from metaerg.datatypes.blast import BlastResult, DBentry, TabularBlastParser
from metaerg.datatypes import ncbi_ftp

DB_DESCRIPTIONS_FILENAME = 'db_descriptions.txt'
DB_TAXONOMY_FILENAME = 'db_taxonomy.txt'
DB_PROTEINS_FILENAME = 'db_protein.faa'
DB_RNA_FILENAME = 'db_rna.fna'

RELEVANT_FEATURE_TYPES = 'CDS rRNA tRNA RNase_P_RNA SRP_RNA riboswitch snoRNA ncRNA tmRNA antisense_RNA binding_site ' \
                         'hammerhead_ribozyme scRNA mobile_element mobile_genetic_element misc_RNA ' \
                         'autocatalytically_spliced_intron'.split()
IGNORED_FEATURE_TYPES = 'gene pseudogene exon direct_repeat region sequence_feature pseudogenic_tRNA pseudogenic_rRNA ' \
                        'repeat_region ribosome_entry_site minus_10_signal minus_35_signal protein_binding_site ' \
                        'regulatory transcript mRNA gap assembly_gap source misc_feature variation misc_structure ' \
                        'transcript 3\'UTR 5\'UTR intron signal_peptide_region_of_CDS sequence_alteration'.split()

TAXONOMY = {'p': {}, 'e': {}, 'v': {}}
DESCRIPTIONS = {'p': {}, 'e': {}, 'v': {}}
ANNOTATOR_KEY = 'diamond_and_blastn'

def _run_programs(genome, contig_dict, db_connection, result_files):
    rna_nt_file = context.spawn_file('rna.fna', genome.name)
    blastn_result_file = context.spawn_file('blastn', genome.name)
    if rna_nt_file.exists() and rna_nt_file.stat().st_size:
        blastn_db = Path(context.DATABASE_DIR, DB_RNA_FILENAME)
        context.run_external(f'blastn -db {blastn_db} -query {rna_nt_file} -out {blastn_result_file} -max_target_seqs 10 '
                             f'-outfmt 6')
    else:
        context.log(f'({genome.name}) Skipping blastn, fasta file with RNA genes missing or empty.')

    cds_aa_file = context.spawn_file('cds.faa', genome.name)
    diamond_db = Path(context.DATABASE_DIR, DB_PROTEINS_FILENAME)
    context.run_external(f'diamond blastp -d {diamond_db} -q {cds_aa_file} -o {result_files[0]} -f 6 '
                         f'--threads {context.CPUS_PER_GENOME} --fast --max-target-seqs 10')


def _read_results(genome, contig_dict, db_connection, result_files) -> int:
    taxon_counts = Counter()
    def process_blast_result(blast_result: BlastResult):
        feature = sqlite.read_feature_by_id(db_connection, blast_result.query())
        if not feature:
            context.log(f'({genome.name}) FATAL ERROR: Found {ANNOTATOR_KEY} result for unknown feature {blast_result.query()}, '
                        f'may need to rerun metaerg with --force')
            raise Exception(f'Found diamond_and_blastn result for unknown feature {blast_result.query()}, '
                            f'may need to rerun metaerg with --force')
        feature.blast = blast_result
        feature.descr = blast_result.hits[0].hit.descr
        feature.taxon =  blast_result.hits[0].hit.taxon
        sqlite.update_feature_in_db(db_connection, feature)
        taxon_counts.update((feature.taxon,))

    def dbentry_from_string(db_id: str) -> DBentry:
        w = db_id.split('~')
        return DBentry(domain=w[0], taxon=TAXONOMY[w[0]][int(w[1])], descr=DESCRIPTIONS[w[0]][int(w[2])],
                       accession=w[3], gene=w[4], length=int(w[5]), pos=int(w[6]))

    blast_result_count = 0
    with TabularBlastParser(result_files[0], 'BLAST', dbentry_from_string) as handle:
        for blast_result in handle:
            blast_result_count += 1
            process_blast_result(blast_result)

    blastn_result_file = context.spawn_file('blastn', genome.name)
    with TabularBlastParser(blastn_result_file, 'BLAST', dbentry_from_string) as handle:
        for blast_result in handle:
            blast_result_count += 1
            process_blast_result(blast_result)
    # (2) update genome
    genome.fraction_classified = taxon_counts.total() / (genome.number_of_features - genome.number_of_crispr_repeats -
                                                         genome.number_of_other_repeats - genome.number_of_retrotransposons)
    for k, v in taxon_counts.most_common(1):
        genome.fraction_classified_to_top_taxon = v / taxon_counts.total()
        genome.top_taxon = k
        break

    return blast_result_count


def preload_db():
    context.log('Parsing blast database indexes...')
    # (1) load databse descriptions
    db_descr = Path(context.DATABASE_DIR, DB_DESCRIPTIONS_FILENAME)
    entries_without_descr = 0
    with open(db_descr) as descr_handle:
        for line in descr_handle:
            words = line.split('\t')
            DESCRIPTIONS[words[0]][int(words[1])] = words[2].strip()
            if not DESCRIPTIONS[words[0]][int(words[1])]:
                entries_without_descr += 1
    context.log(f'Parsed ({len(DESCRIPTIONS["p"])}, {len(DESCRIPTIONS["e"])}, {len(DESCRIPTIONS["v"])})'
                f' gene descriptions from db for (prokaryotes, eukaryotes and viruses) respectively. Empty'
                f' descriptions: {entries_without_descr}.')
    # (2) load database taxonomy
    db_taxonomy = Path(context.DATABASE_DIR, DB_TAXONOMY_FILENAME)
    with open(db_taxonomy) as taxon_handle:
        for line in taxon_handle:
            words = line.split('\t')
            TAXONOMY[words[0]][int(words[1])] = words[2].strip().replace('~', '; ')
    context.log(f'Parsed ({len(TAXONOMY["p"])}, {len(TAXONOMY["e"])}, {len(TAXONOMY["v"])}) '
                f'taxa from db for (prokaryotes, eukaryotes and viruses) respectively.')


@context.register_annotator
def run_and_read_diamond_blastn():
    return ({'pipeline_position': 81,
             'annotator_key': ANNOTATOR_KEY,
             'purpose': 'function prediction and taxonomic classification of genes with diamond and blastn',
             'programs': ('diamond', 'blastn'),
             'databases': (Path(DB_PROTEINS_FILENAME + '.dmnd'),
                           Path(DB_RNA_FILENAME + '.nhr'),
                           Path(DB_DESCRIPTIONS_FILENAME),
                           Path(DB_TAXONOMY_FILENAME)),
             'result_files': ('diamond',),
             'run': _run_programs,
             'read': _read_results,
             'preload_db': preload_db})


def init_pristine_db_dir(dir) -> tuple[Path, Path, Path, Path]:
    fasta_protein_db = Path(dir, DB_PROTEINS_FILENAME)
    fasta_nt_db = Path(dir, DB_RNA_FILENAME)
    descr_db = Path(dir, DB_DESCRIPTIONS_FILENAME)
    taxon_db = Path(dir, DB_TAXONOMY_FILENAME)
    fasta_protein_db.unlink(missing_ok=True)
    fasta_nt_db.unlink(missing_ok=True)
    descr_db.unlink(missing_ok=True)
    taxon_db.unlink(missing_ok=True)
    return fasta_protein_db, fasta_nt_db, descr_db, taxon_db


def update_db_descriptions_get_db_id(description, dictionary, file_handle, kingdom):
    if description.startswith('MAG'):
        colon_index = description.find(':')
        if colon_index > 0:
            description = description[colon_index+2:]
        else:
            description = re.compile('MAG\d+ family').sub('hypothetical', description)
    if 'hypothetical protein' in description:
        description = 'hypothetical protein'
    try:
        descr_id = dictionary[description]
    except KeyError:
        descr_id = len(dictionary)
        dictionary[description] = descr_id
        file_handle.write(f'{kingdom}\t{descr_id}\t{description}\n')
    return descr_id


@context.register_database_installer
def install_viral_database():
    if 'V' not in context.DATABASE_TASKS:
        return
    context.log('Downloading viral refseq from the NCBI...')
    VIR_DB_DIR = Path(context.DATABASE_DIR, 'ncbi-cache', 'vir')
    VIR_DB_DIR.mkdir(exist_ok=True, parents=True)
    if Path(VIR_DB_DIR, DB_DESCRIPTIONS_FILENAME).exists() and not context.FORCE_INSTALLATION_OF_DB:
        context.log(f'Keeping existing viral database in {VIR_DB_DIR}, use --force to overwrite.')
        return
    fasta_protein_db, fasta_nt_db, descr_db, taxon_db = init_pristine_db_dir(VIR_DB_DIR)

    descr_dict = dict()
    taxon_dict = dict()

    pattern = re.compile(r'\s*\[(.+?)]$')  # \]
    with open(fasta_protein_db, 'w') as prot_fasta_out_handle, \
            open(descr_db, 'w') as descr_handle,\
            open(taxon_db, 'w') as taxon_handle:
        for i in range(1,5):
            f = Path(VIR_DB_DIR, f'viral.{i}.protein.gpff.gz')
            if context.FORCE_INSTALLATION_OF_DB or not f.exists() or not f.stat().st_size:
                try:
                    url = 'ftp.ncbi.nlm.nih.gov/refseq/release/viral/viral'
                    context.run_external(f'wget -P {VIR_DB_DIR} {url}.{i}.protein.gpff.gz')
                except:
                    break
            gene_count = 0
            with gzip.open(f, 'rt') as file_handle:
                for gb_record in SeqIO.parse(file_handle, "genbank"):
                    if gb_record.annotations['molecule_type'] != 'protein':
                        context.log("Warning: skipping non-coding gene in viral refseq")
                        continue
                    match = pattern.search(gb_record.description)
                    if match:
                        gb_record.annotations['taxonomy'].append(match.group(1))
                        gb_record.description = gb_record.description[0:match.start()]
                    else:
                        context.log(f'Warning, Failed to parse viral species from {gb_record.description}')
                    taxon_str = "~".join(gb_record.annotations['taxonomy'])
                    try:
                        taxon_id = taxon_dict[taxon_str]
                    except KeyError:
                        taxon_id = len(taxon_dict)
                        taxon_dict[taxon_str] = taxon_id
                        taxon_handle.write(f'v\t{taxon_id}\t{taxon_str}\n')
                    descr_id = update_db_descriptions_get_db_id(gb_record.description, descr_dict, descr_handle, 'v')
                    seq_record = SeqRecord(gb_record.seq)
                    seq_record.id = '{}~{}~{}~{}~{}~{}~{}'.format('v', taxon_id, descr_id,
                                                                  gb_record.id.replace('~', '-'), '',
                                                                  len(seq_record.seq), gene_count)
                    seq_record.description = gb_record.description
                    SeqIO.write(seq_record, prot_fasta_out_handle, "fasta")
                    gene_count += 1
        context.log(f"Successfully downloaded {gene_count} viral proteins from the Viral Refseq database.")


@context.register_database_installer
def install_eukaryote_database():
    if 'E' not in context.DATABASE_TASKS:
        return
    context.log('Downloading taxon-unique eukaryotic genomes using NCBI Datasets...')
    EUK_DB_DIR = Path(context.DATABASE_DIR, 'ncbi-cache', 'euk')
    EUK_DB_DIR.mkdir(exist_ok=True, parents=True)
    if Path(EUK_DB_DIR, DB_DESCRIPTIONS_FILENAME).exists() and not context.FORCE_INSTALLATION_OF_DB:
        context.log(f'Keeping existing eukaryote database in {EUK_DB_DIR}, use --force to overwrite.')
        return
    fasta_protein_db, fasta_nt_db, descr_db, taxon_db = init_pristine_db_dir(EUK_DB_DIR)
    chdir(context.DATABASE_DIR)
    # (1) probe the NCBI for available eukaryotic assemblies
    api_instance = ncbi.datasets.GenomeApi(ncbi.datasets.ApiClient())
    ASSEMBLY_LEVELS = 'Scaffold', 'Chromosome', 'Complete Genome'
    targets = dict()
    taxa_done = dict()
    for tax_name in 'Amoebozoa Ancyromonadida Apusozoa Breviatea CRuMs Cryptophyceae Discoba Glaucocystophyceae ' \
                'Haptista Hemimastigophora Malawimonadida Metamonada Rhodelphea Rhodophyta Sar Aphelida ' \
                'Choanoflagellata Filasterea Fungi Ichthyosporea Rotosphaerida'.split():
        genome_summary = api_instance.assembly_descriptors_by_taxon(taxon=tax_name, page_size=1000)
        count = 0
        if genome_summary.assemblies:
            for assembly in map(lambda d: d.assembly, genome_summary.assemblies):
                if assembly.assembly_accession[:3] == 'GCF':
                    level_of_new_assembly = assembly.assembly_level
                    if assembly.org.tax_id in taxa_done.keys():
                        if level_of_new_assembly in ASSEMBLY_LEVELS and ASSEMBLY_LEVELS.index(level_of_new_assembly) \
                                > taxa_done[assembly.org.tax_id]:
                            taxa_done[assembly.org.tax_id] = ASSEMBLY_LEVELS.index(level_of_new_assembly)
                            targets[assembly.org.tax_id] = assembly.assembly_accession.strip()
                        else:
                            continue
                    else:
                        count += 1
                        taxa_done[assembly.org.tax_id] = ASSEMBLY_LEVELS.index(level_of_new_assembly) \
                            if level_of_new_assembly in ASSEMBLY_LEVELS else 0
                        targets[assembly.org.tax_id] = assembly.assembly_accession.strip()
        try:
            total_count = int(genome_summary.total_count)
        except TypeError:
            total_count = 0
        context.log(f"({count}/{total_count}) taxon-unique assemblies for {tax_name}")
    total_target_count = len(targets)
    context.log(f'Found  {total_target_count} targets for download in total.')
    # (2) Download new assemblies not yet in local cache
    targets = [accession for accession in targets.values() if not Path(EUK_DB_DIR, f'{accession}.gbk.gz').exists()]
    context.log(f'Of those, {total_target_count - len(targets)} already in local cache. ')

    api_response = api_instance.download_assembly_package(targets, exclude_sequence = False,
                                                          include_annotation_type = ['GENOME_GBFF'],
                                                          _preload_content = False)
    zipfile = Path(context.DATABASE_DIR, 'ncbi_genomes.zip')
    zipfile.unlink(missing_ok=True)
    with open(zipfile, 'wb') as f:
        f.write(api_response.data)
    if not zipfile.exists():
        context.log('No success in downloading eukaryote genomes, aborting...')
        return
    context.run_external(f'unzip -qq -o -d {context.DATABASE_DIR} {zipfile}', log_cmd=False)
    # (3) move gbk file into local cache gzipped
    count=0
    for accession in targets:
        count += 1
        src_dir = Path(context.DATABASE_DIR, 'ncbi_dataset', 'data', accession)
        dest_file = Path(EUK_DB_DIR, f'{accession}.gbk.gz')
        if dest_file.exists():
            continue
        context.log(f'({count}/{len(targets)}) Now extracting {accession} into ncbi-cache as {dest_file}...')
        gbff_count = 0
        for file in src_dir.glob('*.gbff'):
            with open(Path(src_dir, file), 'rb') as f_in, gzip.open(dest_file, 'wb') as f_out:
                f_out.writelines(f_in)
            gbff_count += 1
        if gbff_count > 1:
            context.log(f'WARNING: {src_dir} has {gbff_count} gbff files')
        elif gbff_count == 0:
            context.log(f'WARNING: {src_dir} has no gbff files')
    for file in zipfile, Path(context.DATABASE_DIR, '../../README.md'):
        file.unlink(missing_ok=True)
    # (4) build blast databases
    descr_dict = dict()
    taxon_dict = dict()
    with open(fasta_protein_db, 'w') as prot_fasta_out_handle,open(fasta_nt_db, 'w') as rna_fasta_out_handle, \
            open(descr_db, 'w') as descr_handle, open(taxon_db, 'w') as taxon_handle:
        genome_count = 0
        for file in EUK_DB_DIR.glob('*.gbk.gz'):
            genome_count += 1
            context.log(f'({genome_count}/{len(targets)}) Now extracting "{file}"')
            file_path = Path(EUK_DB_DIR, file)
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
                        if feature.type in RELEVANT_FEATURE_TYPES:
                            if feature.type == 'CDS':
                                if not 'translation' in feature.qualifiers:
                                    without_translation_counter += 1
                                    continue
                                seq_record = SeqRecord(Seq(feature.qualifiers['translation'][0]))
                                file_handle = prot_fasta_out_handle
                            else:
                                seq_record = SeqRecord(feature.extract(gb_record.seq))
                                file_handle = rna_fasta_out_handle
                            if 'product' in feature.qualifiers:
                                seq_record.description = feature.qualifiers['product'][0]
                            elif 'gene' in feature.qualifiers:
                                seq_record.description = feature.qualifiers['gene'][0]
                            elif 'note' in feature.qualifiers:
                                seq_record.description = feature.qualifiers['note'][0]
                            elif 'mobile_element_type' in feature.qualifiers:
                                seq_record.description = feature.qualifiers['mobile_element_type'][0]
                            elif 'ncRNA_class' in feature.qualifiers:
                                seq_record.description = feature.qualifiers['ncRNA_class'][0]
                            else:
                                context.log(f"WARNING: No description for feature {feature}")
                                seq_record.description = f'Unknown {feature.type}'
                            descr_id = update_db_descriptions_get_db_id(seq_record.description, descr_dict,
                                                                        descr_handle, 'e')
                            seq_record.id = '{}~{}~{}~{}~{}~{}~{}'.format('e', taxon_id, descr_id,
                                                                          gb_record.id.replace('~', '-'), '',
                                                                          len(seq_record.seq), feature_success_counter)
                            SeqIO.write(seq_record, file_handle, "fasta")
                            feature_success_counter += 1
                        elif feature.type in IGNORED_FEATURE_TYPES:
                            pass
                        else:
                            context.log(f'  Warning: unknown feature "{feature}"')
            context.log(f'  ... extracted ({feature_success_counter}/{feature_total_counter}) proteins and RNAs. '
                        f'{without_translation_counter} proteins skipped for lack of translation.')


@context.register_database_installer
def install_prokaryote_database():
    if 'P' not in context.DATABASE_TASKS:
        return
    context.log('Downloading a gff file from the NCBI FTP server for each prokaryotic genome in gtdbtk...')
    PROK_DB_DIR = Path(context.DATABASE_DIR, 'ncbi-cache', 'pro')
    PROK_DB_DIR.mkdir(exist_ok=True, parents=True)
    PROK_DB_DIR_FAA = Path(PROK_DB_DIR, 'faa')
    PROK_DB_DIR_FAA.mkdir(exist_ok=True)
    PROK_DB_DIR_FNA = Path(PROK_DB_DIR, 'fna')
    PROK_DB_DIR_FNA.mkdir(exist_ok=True)
    if Path(PROK_DB_DIR, DB_DESCRIPTIONS_FILENAME).exists() and not context.FORCE_INSTALLATION_OF_DB:
        context.log(f'Keeping existing prokaryote database in {PROK_DB_DIR}, use --force to overwrite.')
        return
    RNA_DESCR_RE = re.compile(r'\[product=(.+?)]')
    AA_DESCR_RE = re.compile(r'\[(.+?)]$')
    fasta_protein_db, fasta_nt_db, descr_db, taxon_db = init_pristine_db_dir(PROK_DB_DIR)
    descr_dict = {}
    kingdom = 'p'
    # (2) Load and prep list of target taxa from GTDB
    taxa_count = 0
    gtdbtk_taxonomy_file = Path(context.GTDBTK_DIR, 'taxonomy', 'gtdb_taxonomy.tsv')
    # determine line count
    with open(gtdbtk_taxonomy_file) as taxonomy_handle:
        for line in taxonomy_handle:
            taxa_count += 1
    taxon_list = []
    context.log(f'GTDBTK comprises {taxa_count} taxa...')
    with open(gtdbtk_taxonomy_file) as taxonomy_handle, open(descr_db, 'w') as descr_handle:
        success_count = 0
        genomes_in_cache_count = 0
        count = 0
        for line in taxonomy_handle:
            count += 1
            words = line.split()
            accession = words[0][3:]
            taxonomy = re.sub("\w__", "", words[1]).split(';')
            future_faa_file = Path(PROK_DB_DIR_FAA, f'{accession}.faa.gz')
            future_rna_file = Path(PROK_DB_DIR_FNA, f'{accession}.rna.fna.gz')
            taxon = {'id': len(taxon_list), 'accession': accession, 'taxonomy': taxonomy, 'in_local_cache': False}
            taxon_list.append(taxon)
            if future_faa_file.exists() and future_rna_file.exists():
                taxon['in_local_cache'] = True
                download_status = '++'
                genomes_in_cache_count += 1
            # else:
            #     try:
            #         download_status = ''
            #         targets = []
            #         if not future_faa_file.exists():
            #             targets.append(('_protein.faa.gz', future_faa_file))
            #         if not future_rna_file.exists():
            #             targets.append(('_rna_from_genomic.fna.gz', future_rna_file))
            #         success = ncbi_ftp.fetch(accession, targets, context.TEMP_DIR)
            #         if success is not None and future_faa_file in success.keys() and success[future_faa_file]:
            #             download_status += '*'
            #             success_count += 1
            #             taxon['in_local_cache'] = True
            #             genomes_in_cache_count += 1
            #         if success is not None and future_rna_file in success.keys() and success[future_rna_file]:
            #             download_status += '*'
            #         if not download_status:
            #             continue
            #     except EOFError:
            #         context.log('FTP Error at NCBI - resetting connection')
            #         ftp = FTP('ftp.ncbi.nlm.nih.gov')
            #         ftp.login()
            #         continue
            #     except BrokenPipeError:
            #         context.log('FTP Error at NCBI - resetting connection')
            #         ftp = FTP('ftp.ncbi.nlm.nih.gov')
            #         ftp.login()
            #         continue
            context.log(f'({count}/{taxa_count}) {download_status} {future_faa_file.name} {taxonomy}')
            # extract_proteins_and_rna_prok(taxon, descr_dict, PROK_DB_DIR)
            for input, output in ((future_faa_file, fasta_protein_db), (future_rna_file, fasta_nt_db)):
                if not input.exists():
                    continue
                try:
                    with fasta.FastaParser(input, cleanup_seq=True) as reader, open(output, 'a') as writer:
                        gene_counter = 0
                        for f in reader:
                            if m:= RNA_DESCR_RE.search(f['descr']):
                                descr = m.group(1)
                            elif m:= AA_DESCR_RE.search(f['descr']):
                                descr = f['descr'][0:m.start()].strip()
                                descr = descr.replace('MULTISPECIES: ', '')
                            descr_id = update_db_descriptions_get_db_id(descr, descr_dict, descr_handle, kingdom)
                            f['id'] = '{}~{}~{}~{}~{}~{}~{}'.format(kingdom, taxon['id'], descr_id, f['id'].replace('~', '-'),
                                                                 '', len(f['seq']), gene_counter)
                            f['descr'] = descr
                            gene_counter += 1
                            fasta.write_fasta(writer, f)
                except zlib.error:
                    context.log(f'Error parsing {input}, aborted midway, only part of the genes could be used.')
                    pass


    context.log(f'downloaded {int(success_count)} new genomes. Total prok genomes in cache: {genomes_in_cache_count}.')
    # (4) save taxonomy (includes status)
    context.log(f'Now writing taxonomy file for prok.')
    with open(taxon_db, "w") as taxon_handle:
        for t in taxon_list:
            taxon_handle.write(f'p\t{t["id"]}\t{"~".join(t["taxonomy"])}\t{t["in_local_cache"]}\t{t["accession"]}\n')


def compile_databases():
    fasta_protein_db = Path(context.DATABASE_DIR, DB_PROTEINS_FILENAME)
    context.run_external(f'diamond makedb --in {fasta_protein_db} --db {fasta_protein_db}')
    fasta_nt_db = Path(context.DATABASE_DIR, DB_RNA_FILENAME)
    context.run_external(f'makeblastdb -in {fasta_nt_db} -dbtype nucl')

@context.register_database_installer
def format_blast_databases():
    if 'B' not in context.DATABASE_TASKS:
        return
    context.log('Concatenating and formatting blast databases for (prok, euk, vir)...')
    PROK_DB_DIR = Path(context.DATABASE_DIR, 'ncbi-cache', 'pro')
    EUK_DB_DIR = Path(context.DATABASE_DIR, 'ncbi-cache', 'euk')
    VIR_DB_DIR = Path(context.DATABASE_DIR, 'ncbi-cache', 'vir')

    for filename in DB_PROTEINS_FILENAME, DB_RNA_FILENAME, DB_DESCRIPTIONS_FILENAME, DB_TAXONOMY_FILENAME:
        with open(Path(context.DATABASE_DIR, filename), mode="wb") as destination:
            for dir in EUK_DB_DIR, VIR_DB_DIR, PROK_DB_DIR:
                src_file = Path(dir, filename)
                if src_file.exists():
                    with open(src_file, mode="rb") as source:
                        shutil.copyfileobj(source, destination)
    compile_databases()