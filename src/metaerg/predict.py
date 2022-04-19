import re
import shutil
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation

from metaerg import utils
from metaerg import databases

FORCE = False
MULTI_MODE = False
FILENAME_BASE = "pha"
FILTERED_CONTIGS = 0
MASKED_SEQ = 0
TOTAL_SEQ = 0
FEATURE_INFERENCES = ['minced', 'aragorn', 'cmscan', 'ltrharvest', 'tandem-repeat-finder', 'repeatscout', 'prodigal']
FEATURE_ID_TAGS = ['crispr', 'trna', 'rna', 'ltr', 'tr', 'repeat', 'cds']
TRANSLATION_TABLE = 11
SOURCE = 'meta'


def spawn_file(program_name, mag_name):
    if MULTI_MODE:
        dir =  Path(program_name)
        if not dir.exists():
            dir.mkdir()
        elif dir.is_file():
            if FORCE:
                dir.unlink()
                dir.mkdir()
            else:
                raise Exception("Use force to overwrite existing results")
        return Path(dir, mag_name)
    else:
        file = Path(f'{mag_name}.{program_name}')
        if file.exists() and file.is_dir():
            if FORCE:
                shutil.rmtree(file)
        return file


def create_masked_contig_fasta_file(mag_name, contig_dict:dict, exceptions=None, min_mask_length=50):
    if exceptions is None:
        exceptions = set()
    global MASKED_SEQ, TOTAL_SEQ
    (MASKED_SEQ, TOTAL_SEQ) = (0,0)
    masked_fasta_file = spawn_file('masked', mag_name)
    seq_iterator = (mask_seq(record, exceptions=exceptions, min_mask_length=min_mask_length)
                    for record in contig_dict.values())
    SeqIO.write(seq_iterator, masked_fasta_file, "fasta")
    utils.log(f'Masked {MASKED_SEQ/TOTAL_SEQ*100:.1f}% of sequence data.')
    return masked_fasta_file


def mask_seq(record:SeqRecord, exceptions=None, min_mask_length=50):
    global MASKED_SEQ, TOTAL_SEQ
    seq = str(record.seq)
    TOTAL_SEQ += len(seq)
    for f in record.features:
        inference = utils.get_feature_qualifier(f, 'inference').lower()
        if inference in exceptions:
            continue
        if len(f.location) < min_mask_length:
            continue
        fl:FeatureLocation = f.location
        MASKED_SEQ += fl.end - fl.start
        seq = seq[:fl.start] + 'N' * (fl.end - fl.start) + seq[fl.end:]
    return SeqRecord(Seq(seq), id=record.id, description=record.description)


def create_ids(fasta_file:Path, contig_dict):
    for contig in contig_dict.values():
        contig.features = sorted(contig.features, key=lambda f: f.location.start)
        new_id_number = 0
        old_to_new_id_map = dict()
        for feature in contig.features:
            id_tag = 'error'
            utils.set_feature_qualifier(feature, "inference", utils.get_feature_qualifier(feature, "inference").lower())

            for i in range(len(FEATURE_INFERENCES)):
                if FEATURE_INFERENCES[i] in utils.get_feature_qualifier(feature, 'inference'):
                    id_tag = FEATURE_ID_TAGS[i]
                    break
            new_id = f'{contig.id}_{new_id_number:05d}_{id_tag}'
            new_id_number += 1
            old_to_new_id_map[utils.get_feature_qualifier(feature, "id")] = new_id
            utils.set_feature_qualifier(feature, "id", new_id)
            utils.set_feature_qualifier(feature, "protein_id", new_id)
            utils.set_feature_qualifier(feature, "locus_tag", new_id)
        #for feature in contig.features:
        #    utils.set_feature_qualifier(feature, 'parent', old_to_new_id_map[
        #        utils.get_feature_qualifier(feature, "parent")])


def add_homology_search_results_to_feature(blast_result,  contig_dict, alphabet):
    deciph_feat_id = utils.decipher_metaerg_id(blast_result[0])
    target_feature = contig_dict[deciph_feat_id['contig_id']].features[deciph_feat_id['gene_number']]
    gene_function = None
    identical_function_count = 1
    hit_count = 0
    for hit in blast_result[1]:
        hit_count += 1
        deciph_db_id = databases.decipher_database_id(hit['hit_id'], add_to_cache=True)
        if not gene_function:
            gene_function = deciph_db_id["descr_id"]
        elif gene_function == deciph_db_id["descr_id"]:
            identical_function_count += 1
    top_hit = blast_result[1][0]
    deciph_db_id = databases.decipher_database_id(top_hit['hit_id'], add_to_cache=True)
    utils.set_feature_qualifier(target_feature, "taxonomy", deciph_db_id["taxon"])
    utils.set_feature_qualifier(target_feature, "product", f'[{top_hit["aligned_length"]}/{deciph_db_id["length"]}'
                                           f'{alphabet}@{top_hit["percent_id"]}%] [{identical_function_count}/{hit_count}]'
                                           f' {deciph_db_id["descr"]}')



def predict_crisprs_with_minced(mag_name, contig_dict):
    # Sequence tool CRISPR crispr.start crispr.end numRepeat_CRISPR
    # NODE_883_length_34292_cov_11.1505	minced:0.3.2	repeat_region	33687	33944	4	.	.	ID=CRISPR3;bin=91;rpt_family=CRISPR;rpt_type=direct;rpt_unit_seq=GTCGCACTGGGCTTCTAAAGCCCATGAGGATTGAAAC
    # NODE_883_length_34292_cov_11.1505	minced:0.3.2	repeat_unit	33687	33723	1	.	.	ID=DR.CRISPR3.1;Parent=CRISPR3;bin=91
    fasta_file = create_masked_contig_fasta_file(mag_name, contig_dict)
    minced_file = spawn_file('minced', mag_name)

    utils.log('Predicting CRISPR arrays with minced...')
    if not minced_file.exists() or FORCE:
        utils.run_external(f'minced -gffFull {fasta_file} {minced_file}')
    else:
        utils.log("Reusing existing result file.")

    crispr_region_count = 0
    with open(minced_file) as crispr_handle:
        for line in crispr_handle:
            if line.startswith('#'):
                continue
            words = line.split('\t')
            if len(words) < 9:
                continue
            contig: SeqRecord = contig_dict[words[0]]
            if 'repeat_region' == words[2]:
                crispr_region_count += 1
            feature = utils.gff_words_to_seqfeature(words, words[1])
            feature.type = 'crispr_repeat'
            contig.features.append(feature)
    utils.log(f'CRISPR prediction complete. Found {crispr_region_count} repeat regions.')


def predict_trnas_with_aragorn(mag_name, contig_dict):
    # >C1997739
    # 0 genes found
    # >C1997905
    # 2 genes found
    # 1   tRNA-Tyr               [170962,171044]      35      (gta)
    # 2   tRNA-???               [172294,172367]      35      (aa)
    # 1   tRNA-Leu                       [-3,82]      35      (caa)
    # 1   tmRNA*                       c[-31,312]      208,252 ANDNAQTGAVALAA*
    # ...
    # >end    957850 sequences 11875 tRNA genes
    fasta_file = create_masked_contig_fasta_file(mag_name, contig_dict)
    aragorn_file = spawn_file('aragorn', mag_name)

    utils.log('Predicting tRNAs with aragorn...')
    if not aragorn_file.exists() or FORCE:
        utils.run_external(f'aragorn -l -t -gc{TRANSLATION_TABLE} {fasta_file} -w -o {aragorn_file}')
    else:
        utils.log("Reusing existing result file.")

    trna_count = 0
    current_contig = None
    with open(aragorn_file) as aragorn_handle:
        for line in aragorn_handle:
            if line.startswith('>end'):
                break
            if line.startswith('>'):
                current_contig = contig_dict[line[1:].strip()]
            if not current_contig:
                continue
            words = line.split()
            if len(words) < 5:
                continue
            trna_count += 1
            strand = 1
            if words[2].startswith('c'):
                strand = -1
                words[2] = words[2][1:]
            pos_str = words[2][1:-1].split(',')
            start = max(0, int(pos_str[0]) - 1)
            end = min(len(current_contig.seq), int(pos_str[1]))
            f = SeqFeature(location=FeatureLocation(start, end, strand=strand), type='tRNA',
                           qualifiers={'name': [f'{words[1]}-{words[4]}'],
                                       'inference': ['aragorn']})
            current_contig.features.append(f)
    utils.log(f'tRNA prediction complete. Found {trna_count} tRNAs.')


def predict_non_coding_rna_features_with_infernal(mag_name, contig_dict):
    fasta_file = create_masked_contig_fasta_file(mag_name, contig_dict)
    cmscan_file = spawn_file('cmscan', mag_name)

    utils.log('Predicting non-coding RNAs (rRNA, tRNA, CRISPRs, etc.) with infernal (cmscan)...')
    if not cmscan_file.exists() or FORCE:
        utils.run_external(f'cmscan --tblout {cmscan_file} {Path(databases.DBDIR, "Rfam.cm")} {fasta_file}')
    else:
        utils.log("Reusing existing result file.")

    hits = []
    with open(cmscan_file) as hmm_handle:
        for line in hmm_handle:
            if line.startswith('#'):
                continue
            words = line.split()
            if len(words) < 18 or '?' == words[16]:
                continue
            words[17] = ' '.join(words[17:])
            hit = {'query_id': words[2],
                   'hit_id': words[0],
                   'hit_start': int(words[5]),
                   'hit_end': int(words[6]),
                   'query_start': int(words[7]),
                   'query_end': int(words[8]),
                   'query_strand': 1,
                   'partial': words[10],
                   'score': float(words[14]),
                   'evalue': float(words[15]),
                   'descr': words[17]}
            if '-' == words[9]:
                (hit['query_strand'], hit['query_end'], hit['query_start']) = (-1, hit['query_start'], hit['query_end'])
            overlap = None
            for prev_hit in hits:
                if hit['query_id'] == prev_hit['query_id'] and hit['query_start'] < prev_hit['query_end'] \
                        and hit['query_end'] > prev_hit['query_start']:
                    # overlap detected
                    overlap = hit
                    break
            if overlap:
                if hit['evalue'] > overlap['evalue']:
                    hits.remove(overlap)
                    hits.append(hit)
            else:
                hits.append(hit)
    for hit in hits:
        utils.non_coding_rna_hmm_hit_to_seq_record(hit, contig_dict)
    utils.log(f'Non-coding RNA prediction complete. Found {len(hits)} ncRNAs.')


def predict_retrotransposons_with_ltrharvest(mag_name, contig_dict):
    fasta_file = create_masked_contig_fasta_file(mag_name, contig_dict)
    ltr_index_file = spawn_file('ltr_index', mag_name)
    ltr_harvest_file = spawn_file('ltr_harvest', mag_name)

    utils.log('Predicting retrotransposons with genometools/ltrharvest...')
    if not ltr_harvest_file.exists() or FORCE:
        utils.run_external(f'gt suffixerator -db {fasta_file} -indexname {ltr_index_file} -tis -suf -lcp -des -ssp -sds -dna')
        utils.run_external(f'gt ltrharvest -index {ltr_index_file} -gff3 {ltr_harvest_file} -seqids')
        # remove index files
        for file in ltr_harvest_file.parent.glob(f'{mag_name}.ltr_index*'):
            file.unlink()
    else:
        utils.log("Reusing existing result file.")

    retrotransposon_count = 0
    with open(ltr_harvest_file) as ltr_handle:
        for line in ltr_handle:
            if line.startswith('#'):
                continue
            words = line.split('\t')
            if len(words) < 9:
                continue
            #print(words[2])
            if 'repeat_region' == words[2]:
                retrotransposon_count += 1
                contig: SeqRecord = contig_dict[words[0]]
                gbk_feature = utils.gff_words_to_seqfeature(words, 'LTRharvest')
                gbk_feature.type = 'retrotransposon'
                contig.features.append(gbk_feature)
    utils.log(f'Retrotransposon prediction complete. Found {retrotransposon_count} repeat regions.')


def predict_tandem_repeats_with_trf(mag_name, contig_dict):
    fasta_file = create_masked_contig_fasta_file(mag_name, contig_dict)
    trf_file = spawn_file('tandem-repeat-finder', mag_name)

    utils.log('Predicting tandem repeats with trf...')
    if not trf_file.exists() or FORCE:
        with open(trf_file, 'w') as output:
            utils.run_external(f'trf {fasta_file} 2 7 7 80 10 50 500 -d -h -ngs', stdout=output)
    else:
        utils.log("Reusing existing result file.")

    tdr_count = 0
    with open(trf_file) as trf_handle:
        for line in trf_handle:
            if line.startswith("@"):
                contig: SeqRecord = contig_dict[line[1:].strip()]
                continue
            if not contig:
                continue
            words = line.split()
            loc = FeatureLocation(int(words[0]) - 1, int(words[1]))
            qal = {'inference': ['tandem-repeat-finder'],
                   'note': [f'period size {words[2]}; copies {words[3]}']}
            contig.features.append(SeqFeature(location=loc, type='repeat_region', qualifiers=qal))
            tdr_count += 1
    utils.log(f'Tandem repeat prediction complete. Found {tdr_count} repeat regions.')


def predict_remaining_repeats_with_repeatmasker(mag_name, contig_dict):
    fasta_file = create_masked_contig_fasta_file(mag_name, contig_dict)
    lmer_table_file = spawn_file('lmer-table', mag_name)
    repeatscout_file_raw = spawn_file('repeatscout-raw', mag_name)
    repeatscout_file_filtered = spawn_file('repeatscout-filtered', mag_name)
    repeatmasker_file = spawn_file('repeatmasker', mag_name)

    utils.log('Predicting remaining repeats with repeatmasker...')
    if not repeatmasker_file.exists() or FORCE:
        utils.run_external(f'build_lmer_table -sequence {fasta_file} -freq {lmer_table_file}')
        utils.run_external(f'RepeatScout -sequence {fasta_file} -output {repeatscout_file_raw} -freq {lmer_table_file}')
        with open(repeatscout_file_filtered, 'w') as output, open(repeatscout_file_raw) as input:
            utils.run_external('filter-stage-1.prl', stdin=input, stdout=output)
        utils.run_external(f'RepeatMasker -lib {repeatscout_file_filtered} -dir . {fasta_file}')
        # repeatmasker result file will be called '{mag_name}.masked.out', nothing we can do about that
        utils.run_external(f'mv {mag_name}.masked.out {repeatmasker_file}')
        for file in repeatmasker_file.parent.glob(f'{mag_name}.masked*'):
            if file.is_dir():
                shutil.rmtree(file)
            else:
                file.unlink()
    else:
        utils.log("Reusing existing result file.")

    repeat_count = 0
    repeat_hash = dict()
    with open(repeatmasker_file) as repeatmasker_handle:
        for line in repeatmasker_handle:
            words = line.split()
            if len(words) < 11:
                continue
            try:
                contig: SeqRecord = contig_dict[words[4]]
            except KeyError:
                continue
            strand = 1
            if 'C' == words[8]:
                strand = -1
            loc = FeatureLocation(int(words[5]) - 1, int(words[6]), strand=strand)
            qal = {'note': [f'repeat {words[9]}'],
                   'inference': ['repeatscout, repeatmasker']}
            if 'Simple_repeat' == words[10]:
                contig.features.append(SeqFeature(location=loc, type='repeat_region', qualifiers=qal))
                repeat_count += 1
            else:
                try:
                    repeat_list = repeat_hash[words[9]]
                except KeyError:
                    repeat_list = []
                    repeat_hash[words[9]] = repeat_list
                qal['contig'] = [contig]
                repeat_list.append(SeqFeature(location=loc, type='repeat_region', qualifiers=qal))
    for repeat_list in repeat_hash.values():
        if len(repeat_list) >= 10:
            for f in repeat_list:
                utils.get_feature_qualifier(f, "contig").features.append(f)
                utils.set_feature_qualifier(f, "note", utils.get_feature_qualifier(f, "note") + f' (occurs {len(repeat_list)}x)')
                del f.qualifiers['contig']
                repeat_count += 1
    utils.log(f'Remaining repeat prediction complete. Found {repeat_count} repeat regions.')


def predict_coding_sequences_with_prodigal(mag_name, contig_dict):
    fasta_file = create_masked_contig_fasta_file(mag_name, contig_dict)
    #faa_file = Path(fasta_file.stem + '.cds.faa')
    #fna_file = Path(fasta_file.stem + '.cds.fna')
    prodigal_file = spawn_file('prodigal', mag_name)
    c = ''
    if 'meta' != SOURCE:
        c = f'-g {TRANSLATION_TABLE} '
        utils.log(f'... (and using -g {TRANSLATION_TABLE}) ...')

    utils.log(f'Predicting coding sequences with prodigal (using -p {SOURCE}) ...')
    if not prodigal_file.exists() or FORCE:
        utils.run_external(f'prodigal -p {SOURCE} {c}-m -f gff -q -i {fasta_file} -o {prodigal_file}')
    else:
        utils.log("Reusing existing result file.")

    cds_found = 0
    cds_id_re = re.compile(r'^ID=[0-9]+_[0-9]+;')

    with open(prodigal_file) as prodigal_handle:
        for line in prodigal_handle:
            if line.startswith("#"):
                continue
            words = line.split('\t')
            if len(words) < 9:
                continue
            contig: SeqRecord = contig_dict[words[0]]
            words[8] = cds_id_re.sub('', words[8], 1)
            feature = utils.gff_words_to_seqfeature(words, words[1].lower())
            partial = ''
            if int(utils.get_feature_qualifier(feature, 'partial')) > 0:
                partial = ' partial'
            rbsmotif = ''
            if utils.get_feature_qualifier(feature, "rbs_motif") != 'None':
                rbsmotif = f' rbs motif {utils.get_feature_qualifier(feature, "rbs_motif")}'
            rbsspacer = ''
            if utils.get_feature_qualifier(feature, "rbs_spacer") != 'None':
                rbsspacer = f' rbs spacer {utils.get_feature_qualifier(feature, "rbs_spacer")}'
            note = f'start codon {utils.get_feature_qualifier(feature, "start_type")}{partial}{rbsmotif}{rbsspacer}'
            del feature.qualifiers['partial'], feature.qualifiers["rbs_motif"] , feature.qualifiers["rbs_spacer"], \
                feature.qualifiers['start_type']
            contig.features.append(feature)
            cds_found += 1
    utils.log(f'Prediction of coding sequences complete, found {cds_found} CDS.')


def write_gene_files(mag_name, contig_dict):
    cds_aa_file = spawn_file('cds.faa', mag_name)
    rna_nt_file = spawn_file('rna.nt', mag_name)

    utils.log(f'Translating CDS and writing CDS and RNA genes to file...')
    if not cds_aa_file.exists() or not rna_nt_file.exists() or FORCE:
        with open(cds_aa_file, 'w') as faa_handle, open(rna_nt_file, 'w') as fna_handle:
            for contig in contig_dict.values():
                for f in contig.features:
                    if f.type == 'CDS':
                        feature_seq = utils.pad_seq(f.extract(contig)).translate(table=TRANSLATION_TABLE)[:-1]
                        utils.set_feature_qualifier(f, 'translation', feature_seq.seq)
                        feature_seq.id = utils.get_feature_qualifier(f, 'id')
                        feature_seq.description = utils.get_feature_qualifier(f, 'product')
                        if '*' in feature_seq:
                            utils.log(f'Warning, internal stop codon(s) in CDS {utils.get_feature_qualifier(f, "id")} {feature_seq.seq}')
                        SeqIO.write(feature_seq, faa_handle, "fasta")
                    elif f.type in ['tRNA', 'rRNA', 'ncRNA']:
                        feature_seq = f.extract(contig)
                        feature_seq.id = utils.get_feature_qualifier(f, 'id')
                        feature_seq.description = ''
                        SeqIO.write(feature_seq, fna_handle, "fasta")
    else:
        utils.log("Reusing existing gene files.")


def predict_functions_and_taxa_with_diamond(mag_name, contig_dict):
    diamond_file = spawn_file('diamond', mag_name)
    cds_aa_file = spawn_file('cds.faa', mag_name)
    utils.log(f'Performing homology searches with diamond...')
    if not diamond_file.exists() or FORCE:
        utils.run_external(f'diamond blastp -d {Path(databases.DBDIR, "db_protein.faa")} -q {cds_aa_file} -o {diamond_file} -f 6')
    else:
        utils.log("Reusing existing result file.")

    with utils.TabularBlastParser(diamond_file) as handle:
        for blast_result in handle:
            add_homology_search_results_to_feature(blast_result, contig_dict, 'aa')
    utils.log('Diamond search complete.')


def predict_functions_and_taxa_with_blastn(mag_name, contig_dict):
    blastn_file = spawn_file('blastn', mag_name)
    rna_nt_file = spawn_file('rna.nt', mag_name)

    utils.log(f'Performing homology searches with blastn ...')
    if not blastn_file.exists() or FORCE:
        utils.run_external(f'blastn -db {Path(databases.DBDIR, "db_rna.fna")} -query {rna_nt_file} -out {blastn_file} -max_target_seqs 25 -outfmt 6')
    else:
        utils.log("Reusing existing result file.")

    with utils.TabularBlastParser(blastn_file) as handle:
        for blast_result in handle:
            add_homology_search_results_to_feature(blast_result, contig_dict, 'nt')
    utils.log('Blastn search complete.')


def predict_functions_with_cdd(mag_name, contig_dict):
    cds_aa_file = spawn_file('cds.faa', mag_name)
    cdd_file = spawn_file('cdd', mag_name)

    utils.log(f'Performing homology searches with rpsblast/cdd ...')
    if not cdd_file.exists() or FORCE:
        utils.run_external(f'rpsblast -db {Path(databases.DBDIR, "cdd", "Cdd")} -query {cds_aa_file} -out {cdd_file} -outfmt 6 -evalue 1e-7')
    else:
        utils.log("Reusing existing result file.")

    count = 0
    with utils.TabularBlastParser(cdd_file) as handle:
        for blast_result in handle:
            deciph_feat_id = utils.decipher_metaerg_id(blast_result[0])
            target_feature = contig_dict[deciph_feat_id['contig_id']].features[deciph_feat_id['gene_number']]
            top_hit = True
            for h in blast_result[1]:
                cdd_id = int(h["hit_id"][4:])
                databases.CDD_CACHE.add(cdd_id)
                if top_hit:
                    hit_length = abs(h["hit_end"] - h["hit_start"])
                    cdd_item = databases.CDD[cdd_id]
                    cdd_descr = f'{cdd_item[0]}|{cdd_item[1]} {cdd_item[2]}'
                    utils.set_feature_qualifier(target_feature, 'cdd',
                        f'[{hit_length}/{cdd_item[3]}]@{h["percent_id"]:.1f}% [{h["query_start"]}-{h["query_end"]}] {cdd_descr}')
                    top_hit = False
            count += 1
    utils.log(f'RPSBlast CDD search complete. Found hits for {count} proteins.')


def predict_functions_with_antismash(mag_name, contig_dict):
    antismash_dir = spawn_file('antismash', mag_name)
    gbk_file = spawn_file('gbk', mag_name)

    utils.log(f'Performing homology searches with antismash ...')
    if not antismash_dir.exists() or FORCE:
        if antismash_dir.exists():
            shutil.rmtree(antismash_dir)
        utils.run_external(f'antismash --genefinding-tool none --output-dir {antismash_dir} {gbk_file}')
    else:
        utils.log("Reusing existing results.")

    antismash_hit_count = 0
    for f in sorted(antismash_dir.glob("*region*.gbk")):
        with open(f) as handle:
            antismash_region_name = ''
            for gb_record in SeqIO.parse(handle, "genbank"):
                for feature in gb_record.features:
                    if 'region' == feature.type:
                        antismash_region_name = utils.get_feature_qualifier(feature, "rules")
                    elif 'CDS' in feature.type:
                        d_id = utils.decipher_metaerg_id(utils.get_feature_qualifier(feature, "locus_tag"))
                        metaerg_feature = contig_dict[d_id["contig_id"]].features[d_id["gene_number"]]
                        if antismash_region_name:
                            utils.set_feature_qualifier(metaerg_feature, 'antismash_region', antismash_region_name)
                        antismash_gene_function = utils.get_feature_qualifier(feature, "gene_functions")
                        if antismash_gene_function:
                            utils.set_feature_qualifier(metaerg_feature, 'antismash_function', antismash_gene_function)
                        antismash_gene_category =  utils.get_feature_qualifier(feature, "gene_kind")
                        if antismash_gene_category:
                            utils.set_feature_qualifier(metaerg_feature, 'antismash_category', antismash_gene_category)
                        antismash_hit_count += 1
    utils.log(f'Antismash search complete. Found hits for {antismash_hit_count} proteins (CDS).')


def predict_transmembrane_helixes(mag_name, contig_dict):
    tmhmm_file =  spawn_file('tmhmm', mag_name)
    cds_aa_file = spawn_file('cds.faa', mag_name)

    utils.log(f'Discovering transmembrane helixes with tmhmm...')
    if not tmhmm_file.exists() or FORCE:
        with open(tmhmm_file, 'w') as output, open(cds_aa_file) as input:
            utils.run_external('tmhmm', stdin=input, stdout=output)
        for file in tmhmm_file.parent.glob(f'TMHMM_*'):
            if file.is_dir():
                shutil.rmtree(file)
    else:
        utils.log("Reusing existing results.")

    count = 0
    current_feature = None
    current_txt = ""
    feature_tmh_count = 0
    with open(tmhmm_file) as tmhmm_handle:
        for line in tmhmm_handle:
            if line.startswith("#"):
                continue
            words = line.split()
            d_id = utils.decipher_metaerg_id(words[0])
            new_feature = contig_dict[d_id["contig_id"]].features[d_id["gene_number"]]
            if not current_feature:
                current_feature = new_feature
                if 'inside' == words[2]:
                    current_txt = "i,"
                elif 'outside' == words[2]:
                    current_txt = "o,"
            elif current_feature != new_feature:
                if feature_tmh_count:
                    utils.set_feature_qualifier(current_feature, 'transmembrane_helixes', feature_tmh_count)
                    utils.set_feature_qualifier(current_feature, 'tmh_topology', current_txt[:-1])
                    count += 1
                current_feature = new_feature
                feature_tmh_count = 0
                if 'inside' == words[2]:
                    current_txt = "i,"
                elif 'outside' == words[2]:
                    current_txt = "o,"
            if "TMhelix" == words[2]:
                feature_tmh_count += 1
                current_txt += f'{words[3]}-{words[4]},'
        if feature_tmh_count:
            utils.set_feature_qualifier(current_feature, 'transmembrane_helixes', feature_tmh_count)
            utils.set_feature_qualifier(current_feature, 'tmh_topology', current_txt[:-1])
            count += 1
    utils.log(f'Transmembrane helix discovery complete. Found {count} membrane proteins.')


def predict_signal_peptides(mag_name, contig_dict):
    signalp_dir = spawn_file('signalp', mag_name)
    cds_aa_file = spawn_file('cds.faa', mag_name)

    utils.log(f'Discovering signal peptides with signalp...')
    if not signalp_dir.exists() or FORCE:
        if signalp_dir.exists():
            shutil.rmtree(signalp_dir)
        utils.run_external(f'signalp6 --fastafile {cds_aa_file} --output_dir {signalp_dir} --format none --organism other')
    else:
        utils.log("Reusing existing results.")

    count = 0
    with open(Path(signalp_dir, 'prediction_results.txt')) as signalp_handle:
        for line in signalp_handle:
            if line.startswith("#"):
                continue
            words = line.split()
            if "OTHER" == words[1]:
                continue
            d_id = utils.decipher_metaerg_id(words[0])
            feature = contig_dict[d_id["contig_id"]].features[d_id["gene_number"]]
            utils.set_feature_qualifier(feature, "signal_peptide", words[1])
            count += 1
        utils.log(f'Signal peptide discovery complete. Found {count} proteins with signal peptide.')


def write_databases(mag_name, contig_dict):
    databases.save_cache(Path("."))

