import re
import gzip
from os import chdir
from ftplib import FTP
from pathlib import Path

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import ncbi.datasets

from metaerg.run_and_read.data_model import MetaergGenome, MetaergSeqFeature, BlastResult, DBentry, TabularBlastParser,\
    FeatureType
from metaerg.run_and_read.context import register_annotator, spawn_file, run_external, DATABASE_DIR, CPUS_PER_GENOME, \
    log, register_database_installer, GTDBTK_DIR

DB_DESCRIPTIONS_FILENAME = 'db_descriptions.txt'
DB_TAXONOMY_FILENAME = 'db_taxonomy.txt'
DB_PROTEINS_FILENAME = 'db_protein.faa'
DB_RNA_FILENAME = 'db_rna.fna'


def _run_programs(genome:MetaergGenome, result_files):
    cds_aa_file = spawn_file('cds.faa', genome.id)
    rna_nt_file = spawn_file('rna.nt', genome.id)
    blastn_db = Path(DATABASE_DIR, DB_RNA_FILENAME)
    diamond_db = Path(DATABASE_DIR, DB_PROTEINS_FILENAME)
    run_external(f'diamond blastp -d {diamond_db} -q {cds_aa_file} -o {result_files[0]} -f 6 '
                 f'--threads {CPUS_PER_GENOME} --max-target-seqs 10')
    run_external(f'blastn -db {blastn_db} -query {rna_nt_file} -out {result_files[1]} -max_target_seqs 10 -outfmt 6')


def _read_results(genome:MetaergGenome, result_files) -> int:
    # (1) load databse descriptions
    db_descr = Path(DATABASE_DIR, DB_DESCRIPTIONS_FILENAME)
    descriptions = {'p': {}, 'e': {}, 'v': {}}
    with open(db_descr) as descr_handle:
        for line in descr_handle:
            words = line.split('\t')
            descriptions[words[0]][int(words[1])] = words[2].strip()
    log(f'Parsed ({len(descriptions["p"])}, {len(descriptions["e"])}, {len(descriptions["v"])}) '
        f'gene descriptions from db for (prokaryotes, eukaryotes and viruses) respectively. ')
    # (2) load database taxonomy
    db_taxonomy = Path(DATABASE_DIR, DB_TAXONOMY_FILENAME)
    taxonomy = {'p': {}, 'e': {}, 'v': {}}
    with open(db_taxonomy) as taxon_handle:
        for line in taxon_handle:
            words = line.split('\t')
            taxonomy[words[0]][int(words[1])] = words[2].strip().replace('~', '; ')
    log(f'Parsed ({len(taxonomy["p"])}, {len(taxonomy["e"])}, {len(taxonomy["v"])}) '
        f'taxa from db for (prokaryotes, eukaryotes and viruses) respectively.')
    # (3) parse diamond blast results
    def get_db_entry(db_id):
        words = db_id.split('~')  # org_acc gene_acc [pev] gene# decr# taxon#
        return DBentry(db_id, '', descriptions[words[2]][int(words[4])], taxonomy[words[2]][int(words[5])],
                       int(words[6]), int(words[3]))

    def process_blast_result(blast_result: BlastResult):
        feature: MetaergSeqFeature = genome.get_feature(blast_result.query())
        feature.blast = blast_result
        feature.product = blast_result.summary()
        genome.subsystems.match(feature, (h.hit.descr for h in blast_result.hits if h.aligned_length / h.hit.length >= 0.8))

    blast_result_count = 0
    with TabularBlastParser(result_files[0], 'BLAST', get_db_entry) as handle:
        for blast_result in handle:
            blast_result_count += 1
            process_blast_result(blast_result)
    with TabularBlastParser(result_files[1], 'BLAST', get_db_entry) as handle:
        for blast_result in handle:
            blast_result_count += 1
            process_blast_result(blast_result)
    return blast_result_count


@register_annotator
def run_and_read_diamond_blastn():
    return ({'pipeline_position': 81,
             'purpose': 'function prediction and taxonomic classification of genes with diamond and blastn',
             'programs': ('diamond', 'blastn'),
             'databases': (DB_PROTEINS_FILENAME, DB_RNA_FILENAME, DB_DESCRIPTIONS_FILENAME, DB_TAXONOMY_FILENAME),
             'result_files': ('diamond', 'blastn'),
             'run': _run_programs,
             'read': _read_results})


def init_pristine_db_dir(dir) -> tuple[Path, Path, Path, Path]:
    fasta_protein_db = Path(dir, DB_PROTEINS_FILENAME)
    fasta_nt_db = Path(dir, DB_RNA_FILENAME)
    descr_db = Path(dir, DB_DESCRIPTIONS_FILENAME)
    taxon_db = Path(dir, DB_TAXONOMY_FILENAME)
    fasta_protein_db.unlink()
    fasta_nt_db.unlink()
    descr_db.unlink()
    taxon_db.unlink()
    return fasta_protein_db, fasta_nt_db, descr_db, taxon_db


def update_descriptions_and_get_id(description, dictionary, file_handle, kingdom):
    try:
        descr_id = dictionary[description]
    except KeyError:
        descr_id = len(dictionary)
        dictionary[description] = descr_id
        file_handle.write(f'{kingdom}\t{descr_id}\t{description}\n')
    return descr_id


@register_database_installer
def install_viral_database():
    vir_db_dir = Path(DATABASE_DIR, 'ncbi-cache', 'vir')
    vir_db_dir.mkdir(exist_ok=True, parents=True)
    fasta_protein_db, fasta_nt_db, descr_db, taxon_db = init_pristine_db_dir(vir_db_dir)

    descr_dict = dict()
    taxon_dict = dict()

    pattern = re.compile(r'\s*\[(.+?)]$')  # \]
    with open(fasta_protein_db, 'w') as prot_fasta_out_handle, \
            open(descr_db, 'w') as descr_handle,\
            open(taxon_db, 'w') as taxon_handle:
        for i in range(1,5):
            f = Path(vir_db_dir, f'viral.{i}.protein.gpff.gz')
            if not f.exists() or not f.stat().st_size:
                url = 'https://ftp.ncbi.nlm.nih.gov/refseq/release/viral/viral'
                run_external(f'wget -P {vir_db_dir} {url}.{i}.protein.gpff.gz')
            gene_count = 0
            with gzip.open(f, 'rt') as file_handle:
                for gb_record in SeqIO.parse(file_handle, "genbank"):
                    if gb_record.annotations['molecule_type'] != 'protein':
                        log("warning: skipping non-coding gene in viral refseq")
                        continue
                    match = pattern.search(gb_record.description)
                    if match:
                        gb_record.annotations['taxonomy'].append(match.group(1))
                        gb_record.description = gb_record.description[0:match.start()]
                    else:
                        log(f'Warning, Failed to parse viral species from {gb_record.description}')
                    taxon_str = "~".join(gb_record.annotations['taxonomy'])
                    try:
                        taxon_id = taxon_dict[taxon_str]
                    except KeyError:
                        taxon_id = len(taxon_dict)
                        taxon_dict[taxon_str] = taxon_id
                        taxon_handle.write(f'v\t{taxon_id}\t{taxon_str}\n')
                    descr_id = update_descriptions_and_get_id(gb_record.description, descr_dict, descr_handle, 'v')
                    seq_record = SeqRecord(gb_record.seq)
                    seq_record.id = f'{taxon_id}~{gb_record.id}~v~{gene_count}~{descr_id}~{taxon_id}~{len(seq_record.seq)}'
                    seq_record.description = gb_record.description
                    SeqIO.write(seq_record, prot_fasta_out_handle, "fasta")
                    gene_count += 1


RELEVANT_RNA_GENES = 'rRNA tRNA RNase_P_RNA SRP_RNA riboswitch snoRNA ncRNA tmRNA antisense_RNA binding_site ' \
                     'hammerhead_ribozyme scRNA mobile_element mobile_genetic_element misc_RNA autocatalytically_spliced_intron'.split()
IGNORED_FEATURES = 'gene pseudogene exon direct_repeat region sequence_feature pseudogenic_tRNA pseudogenic_rRNA ' \
                   'repeat_region ribosome_entry_site minus_10_signal minus_35_signal protein_binding_site regulatory' \
                   ' transcript mRNA gap assembly_gap source misc_feature variation misc_structure transcript' \
                   ' 3\'UTR 5\'UTR intron signal_peptide_region_of_CDS sequence_alteration'.split()


@register_database_installer
def install_eukaryote_database():
    euk_db_dir = Path(DATABASE_DIR, 'ncbi-cache', 'euk')
    euk_db_dir.mkdir(exist_ok=True, parents=True)
    fasta_protein_db, fasta_nt_db, descr_db, taxon_db = init_pristine_db_dir(euk_db_dir)
    chdir(DATABASE_DIR)
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
        log(f"({count}/{total_count}) taxon-unique assemblies for {tax_name}")
    total_target_count = len(targets)
    log(f'Found  {total_target_count} targets for download in total.')
    # (2) Download new assemblies not yet in local cache
    targets = [accession for accession in targets.values() if not Path(euk_db_dir, f'{accession}.gbk.gz').exists()]
    log(f'Of those, {total_target_count - len(targets)} already in local cache. ')

    api_response = api_instance.download_assembly_package(targets, exclude_sequence = False,
                                                          include_annotation_type = ['GENOME_GBFF'],
                                                          _preload_content = False)
    zipfile = Path(DATABASE_DIR, 'ncbi_genomes.zip')
    zipfile.unlink(missing_ok=True)
    with open(zipfile, 'wb') as f:
        f.write(api_response.data)
    if not zipfile.exists():
        log('No success in downloading eukaryote genomes, aborting...')
        return
    run_external(f'unzip -qq -o -d {DATABASE_DIR} {zipfile}', log_cmd=False)
    # (3) move gbk file into local cache gzipped
    count=0
    for accession in targets:
        count += 1
        src_dir = Path(DATABASE_DIR, 'ncbi_dataset', 'data', accession)
        dest_file = Path(euk_db_dir, f'{accession}.gbk.gz')
        if dest_file.exists():
            continue
        log(f'({count}/{len(targets)}) Now extracting {accession} into ncbi-cache as {dest_file}...')
        gbff_count = 0
        for file in src_dir.glob('*.gbff'):
            with open(Path(src_dir, file), 'rb') as f_in, gzip.open(dest_file, 'wb') as f_out:
                f_out.writelines(f_in)
            gbff_count += 1
        if gbff_count > 1:
            log(f'WARNING: {src_dir} has {gbff_count} gbff files')
        elif gbff_count == 0:
            log(f'WARNING: {src_dir} has no gbff files')
    for file in zipfile, Path(DATABASE_DIR, '../../README.md'):
        file.unlink()
    # (4) build blast databases
    descr_dict = dict()
    taxon_dict = dict()
    with open(fasta_protein_db, 'w') as prot_fasta_out_handle,open(fasta_nt_db, 'w')  as rna_fasta_out_handle, \
            open(descr_db, 'w') as desrc_handle, open(taxon_db, 'w') as taxon_handle:
        genome_count = 0
        for file in euk_db_dir.glob('*.gbk.gz'):
            genome_count += 1
            log(f'({genome_count}/{len(targets)}) Now extracting "{file}"')
            file_path = Path(euk_db_dir, file)
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
                        if feature.type == 'CDS' or feature.type in RELEVANT_RNA_GENES:
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
                                log(f"WARNING: No description for feature {feature}")
                                seq_record.description = f'Unknown {feature.type}'
                            descr_id = update_descriptions_and_get_id(seq_record.description, descr_dict, desrc_handle,
                                                                      'e')
                            seq_record.id = f'{gb_record.id}~0~e~{feature_success_counter}~{descr_id}~{taxon_id}~{len(seq_record.seq)}'
                            SeqIO.write(seq_record, file_handle, "fasta")
                            feature_success_counter += 1
                        elif feature.type in IGNORED_FEATURES:
                            pass
                        else:
                            log(f'  Warning: unknown feature "{feature}"')
            log(f'  ... extracted ({feature_success_counter}/{feature_total_counter}) proteins and RNAs. '
                f'{without_translation_counter} proteins skipped for lack of translation.')


@register_database_installer
def install_prokaryote_database():
    prok_db_dir = Path(DATABASE_DIR, 'ncbi-cache', 'euk')
    prok_db_dir.mkdir(exist_ok=True, parents=True)
    fasta_protein_db, fasta_nt_db, descr_db, taxon_db = init_pristine_db_dir(prok_db_dir)
    descr_dict = {}
    # (2) Load and prep list of target taxa from GTDB
    taxa_count = 0
    gtdbtk_taxonomy_file = Path(GTDBTK_DIR, 'taxonomy', 'gtdb_taxonomy.tsv')
    # determine line count for tqdm
    with open(gtdbtk_taxonomy_file) as taxonomy_handle:
        for line in taxonomy_handle:
            taxa_count += 1
    taxon_list = []
    with open(gtdbtk_taxonomy_file) as taxonomy_handle:
        # determine line count for tqdm
        ftp = FTP('ftp.ncbi.nlm.nih.gov')
        ftp.login()
        success_count = 0
        genomes_in_cache_count = 0
        count = 0
        for line in taxonomy_handle:
            count += 1
            words = line.split()
            accession = words[0][3:]
            taxonomy = re.sub("\w__", "", words[1]).split(';')
            seq_file = Path(GTDBTK_DIR, 'fastani', 'database', accession[0:3], accession[4:7],
                            accession[7:10], accession[10:13], f'{accession}_genomic.fna.gz')
            future_gff_file = Path(prok_db_dir, f'{accession}.gff.gz')
            taxon = {'id': len(taxon_list), 'accession': accession, 'taxonomy': taxonomy,
                    'cached_gff_file': future_gff_file, 'gtdb_seq_file': seq_file, 'in_local_cache': False}
            taxon_list.append(taxon)
            download_status = " "
            if future_gff_file.exists():
                taxon['in_local_cache'] = True
                download_status = "+"
            else:
                continue
            genomes_in_cache_count += 1
            log(f'({count}/{taxa_count}) {download_status} {future_gff_file.name} {taxonomy}')
            extract_proteins_and_rna_prok(taxon, descr_dict, prok_db_dir)
            log(f'downloaded {success_count} new genomes. Total genomes in cache: {genomes_in_cache_count}.')
    # (4) save taxonomy (includes status)
    log(f'Writing taxonomy file.')
    with open(taxon_db, "w") as taxon_handle:
        for t in taxon_list:
            taxon_handle.write(f'p\t{t["id"]}\t{"~".join(t["taxonomy"])}\t{t["in_local_cache"]}\t{t["accession"]}\n')


def extract_proteins_and_rna_prok(taxon, descr_dict, prok_db_dir):
    try:
        genome: MetaergGenome = MetaergGenome(contig_file=taxon["gtdb_seq_file"], rename_contigs=False,
                                              min_contig_length=0, id =taxon["gtdb_seq_file"].name)
        with gzip.open(taxon["cached_gff_file"], "rt") as gff_handle, \
                open(Path(prok_db_dir, DB_PROTEINS_FILENAME), 'a') as prot_fasta_out_handle, \
                open(Path(prok_db_dir, DB_RNA_FILENAME), 'a') as rna_fasta_out_handle, \
                open(Path(prok_db_dir, DB_DESCRIPTIONS_FILENAME), 'a') as description_handle:
            gene_counter = 0
            line: str
            for line in gff_handle:
                if line.startswith("#"):
                    continue
                if 'similar to' in line or 'identical to' in line:
                    continue
                words = line.split('\t')
                if len(words) < 9:
                    continue
                contig = genome.contigs[words[0]]
                if words[2] =='CDS'  or words[2] in RELEVANT_RNA_GENES:
                    feature: MetaergSeqFeature = contig.spawn_feature(int(words[3]) - 1,
                                                                      int(words[4]),
                                                                      -1 if words[6] == '-' else 1,
                                                                      FeatureType(words[2]))
                    gene_counter += 1
                    feature = utils.gff_words_to_seqfeature(words)
                    translation_table = int(feature.qualifiers['transl_table'])
                    try:
                        contig = contig_dict[words[0]]
                    except KeyError:
                        continue
                    if 'CDS' == words[2] or words[2] in RELEVANT_RNA_GENES:
                        seq = run_and_read.data_model.pad_seq(feature.extract(contig))
                        seq = seq.translate(table=translation_table)
                        seq = seq[:-1]
                        if '*' in seq:
                            continue
                        handle = prot_fasta_out_handle
                    else:
                        seq = feature.extract(contig)
                        handle = rna_fasta_out_handle
                    try:
                        seq.description = feature.qualifiers["product"]
                    except KeyError:
                        if "note" in feature.qualifiers.keys():
                            seq.description = feature.qualifiers["note"]
                        else:
                            continue
                    descr_id = update_descriptions_and_get_id(seq.description, descr_dict, description_handle, 'p')
                    seq.id = f'{taxon["accession"]}~{feature.qualifiers["id"]}~p~{gene_counter}~{descr_id}~{taxon["id"]}~{len(seq.seq)}'
                    SeqIO.write(seq, handle, "fasta")
                elif words[2] in IGNORED_FEATURES:
                    pass
                else:
                    log(f'  Warning: unknown feature type in line "{line}"')
    except EOFError:
        log(f'Incomplete file detected for {taxon["cached_gff_file"]}. Deleting...')
        taxon["cached_gff_file"].unlink()