import re
import shutil
from pathlib import Path

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation

from metaerg import utils
from metaerg import databases

FILTERED_CONTIGS = 0
MASKED_SEQ = 0
TOTAL_SEQ = 0
VALID_GBK_FEATURE_KEYS = set('CDS tRNA rRNA ncRNA repeat_region misc_feature'.split())
FEATURE_INFERENCES = ['minced', 'aragorn', 'cmscan', 'ltrharvest', 'tandem-repeat-finder', 'repeatscout', 'prodigal']
FEATURE_ID_TAGS = ['crispr', 'trna', 'rna', 'ltr', 'tr', 'repeat', 'cds']
FEATURE_ID_PATTERN = re.compile("(.+)_(\d{5})_(crispr|trna|rna|ltr|tr|repeat|cds)$")
NON_CODING_RNA_TYPES = {'LSU_rRNA_bacteria':'rRNA',
                        'LSU_rRNA_archaea': 'rRNA',
                        'LSU_rRNA_eukarya': 'rRNA',
                        'SSU_rRNA_bacteria': 'rRNA',
                        'SSU_rRNA_archaea': 'rRNA',
                        'SSU_rRNA_eukarya': 'rRNA',
                        'SSU_rRNA_microsporidia': 'rRNA',
                        '5S_rRNA': 'rRNA',
                        '5_8S_rRNA': 'rRNA',
                        'tmRNA': 'tmRNA',
                        'tRNA': 'tRNA'}
TRANSLATION_TABLE = 11
SOURCE = 'meta'


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
        for feature in contig.features:
            utils.set_feature_qualifier(feature, 'parent', old_to_new_id_map[
                utils.get_feature_qualifier(feature, "parent")])


def decipher_metaerg_id(id):
    m = FEATURE_ID_PATTERN.match(id)
    return {'contig_id': m.group(1),
            'gene_number': int(m.group(2)),
            'gene_type': m.group(3)}


def gff_words_to_seqfeature(words: list, inference):
    seq_feature = utils.gff_words_to_seqfeature(words)
    utils.set_feature_qualifier(seq_feature, "inference", inference)

    if not words[2] in VALID_GBK_FEATURE_KEYS:
        utils.set_feature_qualifier(seq_feature, "note", words[2])
        seq_feature.type = 'misc_feature'

    for unwanted in 'gc_cont conf score cscore sscore rscore uscore tscore rpt_type rpt_family ltr_similarity ' \
                    'seq_number'.split():
        try:
            del seq_feature.qualifiers[unwanted]
        except KeyError:
            pass

    return seq_feature


def predict_crisprs_with_minced(fasta_file:Path, contig_dict):
    # Sequence tool CRISPR crispr.start crispr.end numRepeat_CRISPR
    # NODE_883_length_34292_cov_11.1505	minced:0.3.2	repeat_region	33687	33944	4	.	.	ID=CRISPR3;bin=91;rpt_family=CRISPR;rpt_type=direct;rpt_unit_seq=GTCGCACTGGGCTTCTAAAGCCCATGAGGATTGAAAC
    # NODE_883_length_34292_cov_11.1505	minced:0.3.2	repeat_unit	33687	33723	1	.	.	ID=DR.CRISPR3.1;Parent=CRISPR3;bin=91
    utils.create_masked_contig_fasta_file(fasta_file, contig_dict)
    minced_file = Path(fasta_file.stem + '.minced')
    utils.log('Predicting CRISPR arrays with minced...')
    utils.run_external(f'minced -gffFull {fasta_file} {minced_file}')
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
            feature = gff_words_to_seqfeature(words, words[1])
            contig.features.append(feature)
    utils.log(f'CRISPR prediction complete. Found {crispr_region_count} repeat regions.')


def predict_trnas_with_aragorn(fasta_file:Path, contig_dict):
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
    utils.create_masked_contig_fasta_file(fasta_file, contig_dict)
    aragorn_file = Path(fasta_file.stem + '.aragorn')
    utils.log('Predicting tRNAs with aragorn...')
    utils.run_external(f'aragorn -l -t -gc{TRANSLATION_TABLE} {fasta_file} -w -o {aragorn_file}')
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


def non_coding_rna_hmm_hit_to_seq_record(hit, contig_dict):
    contig = contig_dict[hit['query_id']]
    type = 'ncRNA'
    try:
        type = NON_CODING_RNA_TYPES[hit['hit_id']]
    except KeyError:
        if hit['hit_id'].startswith('CRISPR'):
            type = 'crispr'

    f = SeqFeature(location=FeatureLocation(hit['query_start'] - 1, hit['query_end'], strand=hit['query_strand']),
                   type=type, qualifiers={'name': [hit["hit_id"]],
                                          'inference': ['cmscan'],
                                          'profile': [hit["descr"]]})
    contig.features.append(f)


def predict_non_coding_rna_features_with_infernal(fasta_file:Path, contig_dict):
    utils.create_masked_contig_fasta_file(fasta_file, contig_dict)
    cmscan_file = Path(fasta_file.stem + '.cmscan')
    utils.log('Predicting non-coding RNAs (rRNA, tRNA, CRISPRs, etc.) with infernal (cmscan)...')
    utils.run_external(f'cmscan --tblout {cmscan_file} {Path(databases.DBDIR, "Rfam.cm")} {fasta_file}')
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
        non_coding_rna_hmm_hit_to_seq_record(hit, contig_dict)
    utils.log(f'Non-coding RNA prediction complete. Found {len(hits)} ncRNAs.')


def predict_retrotransposons_with_ltrharvest(fasta_file, contig_dict):
    utils.create_masked_contig_fasta_file(fasta_file, contig_dict)
    ltr_index = Path(fasta_file.stem + '.ltr_index')
    ltr_harvest_file = Path(fasta_file.stem + '.ltr_harvest.gff')
    utils.log('Predicting retrotransposons with genometools/ltrharvest...')
    utils.run_external(f'gt suffixerator -db {fasta_file} -indexname {ltr_index} -tis -suf -lcp -des -ssp -sds -dna')
    utils.run_external(f'gt ltrharvest -index {ltr_index} -gff3 {ltr_harvest_file} -seqids')
    retrotransposon_count = 0
    with open(ltr_harvest_file) as ltr_handle:
        for line in ltr_handle:
            if line.startswith('#'):
                continue
            words = line.split('\t')
            if len(words) < 9:
                continue
            contig: SeqRecord = contig_dict[words[0]]
            if 'repeat_region' == words[2]:
                retrotransposon_count += 1
            contig.features.append(gff_words_to_seqfeature(words, words[1]))
    utils.log(f'Retrotransposon prediction complete. Found {retrotransposon_count} repeat regions.')


def predict_tandem_repeats_with_trf(fasta_file, contig_dict):
    utils.create_masked_contig_fasta_file(fasta_file, contig_dict)
    trf_file = Path(fasta_file.stem + '.trf')
    utils.log('Predicting tandem repeats with trf...')
    if not utils.SILENT:
        with open(trf_file, 'w') as output:
            utils.run_external(f'trf {fasta_file} 2 7 7 80 10 50 500 -d -h -ngs', stdout=output)
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


def predict_remaining_repeats_with_repeatmasker(fasta_file: Path, contig_dict):
    utils.create_masked_contig_fasta_file(fasta_file, contig_dict)
    utils.log('Predicting remaining repeats with repeatmasker...')
    lmer_table = Path(fasta_file.stem + '.lmers')
    repeat_scout_raw = Path(fasta_file.stem + '.repeatscout')
    repeat_scout_filtered = Path(fasta_file.stem + '.repeatscout.filtered')
    repeat_masker_file = Path(fasta_file.stem + '.fna.out')

    utils.run_external(f'build_lmer_table -sequence {fasta_file} -freq {lmer_table}')
    utils.run_external(f'RepeatScout -sequence {fasta_file} -output {repeat_scout_raw} -freq {lmer_table}')
    with open(repeat_scout_filtered, 'w') as output, open(repeat_scout_raw) as input:
        utils.run_external('filter-stage-1.prl', stdin=input, stdout=output)
    utils.run_external(f'RepeatMasker -lib {repeat_scout_filtered} -dir . {fasta_file}')
    repeat_count = 0
    repeat_hash = dict()
    with open(repeat_masker_file) as repeatmasker_handle:
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


def predict_coding_sequences_with_prodigal(fasta_file:Path, contig_dict):
    utils.create_masked_contig_fasta_file(fasta_file, contig_dict)
    utils.log(f'Predicting coding sequences with prodigal (using -p {SOURCE}) ...')
    faa_file = Path(fasta_file.stem + '.cds.faa')
    fna_file = Path(fasta_file.stem + '.cds.fna')
    gff_file = Path(fasta_file.stem + '.cds.gff')
    c = ''
    if 'meta' != SOURCE:
        c = f'-g {TRANSLATION_TABLE} '
        utils.log(f'... (and using -g {TRANSLATION_TABLE}) ...')

    utils.run_external(f'prodigal -p {SOURCE} {c}-m -f gff -q -i {fasta_file} -a {faa_file} -d {fna_file} -o {gff_file}')
    cds_found = 0
    cds_id_re = re.compile(r'^ID=[0-9]+_[0-9]+;')

    with open(gff_file) as prodigal_handle:
        for line in prodigal_handle:
            if line.startswith("#"):
                continue
            words = line.split('\t')
            if len(words) < 9:
                continue
            contig: SeqRecord = contig_dict[words[0]]
            words[8] = cds_id_re.sub('', words[8], 1)
            feature = gff_words_to_seqfeature(words, words[1].lower())
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


def add_homology_search_results_to_feature(blast_result,  contig_dict, alphabet):
    deciph_feat_id = decipher_metaerg_id(blast_result[0])
    target_feature = contig_dict[deciph_feat_id['contig_id']].features[deciph_feat_id['gene_number']]
    gene_function = None
    identical_function_count = 1
    hit_count = 0
    for hit in blast_result[1]:
        hit_count += 1
        deciph_db_id = databases.decipher_database_id(hit['hit_id'])
        if not gene_function:
            gene_function = deciph_db_id["descr_id"]
        elif gene_function == deciph_db_id["descr_id"]:
            identical_function_count += 1
    top_hit = blast_result[1][0]
    deciph_db_id = databases.decipher_database_id(top_hit['hit_id'])
    utils.set_feature_qualifier(target_feature, "taxonomy", deciph_db_id["taxon"])
    utils.set_feature_qualifier(target_feature, "product", f'[{top_hit["aligned_length"]}/{deciph_db_id["length"]}'\
                                           f'{alphabet}@{top_hit["percent_id"]}%] [{identical_function_count}/{hit_count}]'\
                                           f' {deciph_db_id["descr"]}')


def annotate_features_by_homology_diamond(fasta_file:Path, contig_dict):
    utils.log(f'Performing homology searches with diamond...')
    cds_aa_file = Path(fasta_file.stem + '.cds.aa')
    with open(cds_aa_file, 'w') as faa_handle:
        for contig in contig_dict.values():
            for f in contig.features:
                if f.type == 'CDS':
                    feature_seq = utils.pad_seq(f.extract(contig)).translate(table=TRANSLATION_TABLE)[:-1]
                    utils.set_feature_qualifier(f, 'translation', feature_seq.seq)
                    feature_seq.id = utils.get_feature_qualifier(f, 'id')
                    feature_seq.description = ''
                    if '*' in feature_seq:
                        utils.log(f'Warning, internal stop codon(s) in CDS {utils.get_feature_qualifier(f, "id")} {feature_seq.seq}')
                    SeqIO.write(feature_seq, faa_handle, "fasta")
    diamond_file = Path(fasta_file.stem + '.diamond.tab.txt')
    utils.run_external(f'diamond blastp -d {Path(databases.DBDIR, "db_protein.faa")} -q {cds_aa_file} -o {diamond_file} -f 6')
    # load descriptions
    with utils.TabularBlastParser(diamond_file) as handle:
        for blast_result in handle:
            add_homology_search_results_to_feature(blast_result, contig_dict, 'aa')
    utils.log('Diamond search complete.')


def annotate_features_by_homology_blastn(fasta_file:Path, contig_dict):
    utils.log(f'Performing homology searches with blastn ...')
    rna_nt_file = Path(fasta_file.stem + '.rna.nt')
    with open(rna_nt_file, 'w') as fna_handle:
        for contig in contig_dict.values():
            for f in contig.features:
                if f.type in ['tRNA', 'rRNA', 'ncRNA']:
                    feature_seq = f.extract(contig)
                    feature_seq.id = utils.get_feature_qualifier(f, 'id')
                    feature_seq.description = ''
                    SeqIO.write(feature_seq, fna_handle, "fasta")
    blastn_file = Path(fasta_file.stem + '.blastn.tab.txt')
    utils.run_external(f'blastn -db {Path(databases.DBDIR, "db_rna.fna")} -query {rna_nt_file} -out {blastn_file} -max_target_seqs 25 -outfmt 6')
    with utils.TabularBlastParser(blastn_file) as handle:
        for blast_result in handle:
            add_homology_search_results_to_feature(blast_result, contig_dict, 'nt')
    utils.log('Blastn search complete.')


def annotate_features_by_homology_cdd(fasta_file: Path, contig_dict):
    utils.log(f'Performing homology searches with rbsblast/cdd ...')
    cds_aa_file = Path(fasta_file.stem + '.cds.aa')
    cdd_file = Path(fasta_file.stem + '.cdd.tab.txt')
    utils.run_external(f'rpsblast -db {Path(databases.DBDIR, "cdd", "Cdd")} -query {cds_aa_file} -out {cdd_file} -outfmt 6 -evalue 1e-7')
    count = 0
    with utils.TabularBlastParser(cdd_file) as handle:
        for blast_result in handle:
            deciph_feat_id = decipher_metaerg_id(blast_result[0])
            target_feature = contig_dict[deciph_feat_id['contig_id']].features[deciph_feat_id['gene_number']]
            cdd_hit_str = ""
            for h in blast_result[1]:
                cdd_hit_str += f'{h["hit_id"]}:{h["percent_id"]:.1f}%@{h["query_start"]}-{h["query_end"]};'
            utils.set_feature_qualifier(target_feature, 'cdd_hits', cdd_hit_str[:-1])
            count += 1
    utils.log(f'RPSBlast CDD search complete. Found hits for {count} proteins (CDS).')



def annotate_features_by_homology_antismash(fasta_file: Path, contig_dict):
    utils.log(f'Performing homology searches with antismash ...')
    antismash_dir = Path('antismash')
    if antismash_dir.exists():
       shutil.rmtree(antismash_dir)
    #with open("antismash_output", 'w') as output:
    #    utils.run_external(f'antismash --genefinding-tool none --output-dir {antismash_dir} {fasta_file.stem + ".gbk"}', stdout=output, stderr=output)
    utils.run_external(f'antismash --genefinding-tool none --output-dir {antismash_dir} {fasta_file.stem + ".gbk"}')
    antismash_hit_count = 0
    for f in sorted(antismash_dir.glob("*region*.gbk")):
        with open(f) as handle:
            utils.log(f.name)
            antismash_region_name = '[antismash]'
            for gb_record in SeqIO.parse(handle, "genbank"):
                for feature in gb_record.features:
                    print(feature)
                    if 'region' == feature.type:
                        antismash_region_name = utils.get_feature_qualifier(feature, "rules")
                        utils.log(f'antismash_region_name "{antismash_region_name}"')
                    elif 'CDS' in feature.type:
                        d_id = decipher_metaerg_id(utils.get_feature_qualifier(feature, "locus_tag"))
                        metaerg_feature = contig_dict[d_id["contig_id"]].features[d_id["gene_number"]]
                        utils.log('metaerg feature:')
                        utils.set_feature_qualifier(metaerg_feature, 'antismash_region', antismash_region_name)
                        antismash_gene_function = utils.get_feature_qualifier(feature, "gene_functions")
                        if antismash_gene_function:
                            utils.set_feature_qualifier(metaerg_feature, 'antismash_function', antismash_gene_function)
                        antismash_gene_category =  utils.get_feature_qualifier(feature, "gene_kind")
                        if antismash_gene_category:
                            utils.set_feature_qualifier(metaerg_feature, 'antismash_category', antismash_gene_category)
                        print(metaerg_feature)
                        antismash_hit_count += 1
    utils.log(f'Antismash search complete. Found hits for {antismash_hit_count} proteins (CDS).')


def discover_transmembrane_helixes(fasta_file: Path, contig_dict):
    utils.log(f'Discovering transmembrane helixes with tmhmm...')
    cds_aa_file = Path(fasta_file.stem + '.cds.aa')
    tmhmm_file =  Path(fasta_file.stem + '.tmhmm.txt')
    with open(tmhmm_file, 'w') as output, open(cds_aa_file) as input:
        utils.run_external('tmhmm', stdin=input, stdout=output)
    count = 0
    current_feature = None
    current_txt = ""
    feature_tmh_count = 0
    with open(tmhmm_file) as tmhmm_handle:
        for line in tmhmm_handle:
            if line.startswith("#"):
                continue
            words = line.split()
            d_id = decipher_metaerg_id(words[0])
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


def discover_signal_peptides(fasta_file: Path, contig_dict):
    utils.log(f'Discovering signal peptides with signalp...')
    cds_aa_file = Path(fasta_file.stem + '.cds.aa')
    signalp_dir = Path('signalp')
    if signalp_dir.exists():
        shutil.rmtree(signalp_dir)
    utils.run_external(f'signalp6 --fastafile {cds_aa_file} --output_dir {signalp_dir} --format none --organism other')
    count = 0
    with open(Path('signalp', 'prediction_results.txt')) as signalp_handle:
        for line in signalp_handle:
            if line.startswith("#"):
                continue
            words = line.split()
            if "OTHER" == words[1]:
                continue
            d_id = decipher_metaerg_id(words[0])
            feature = contig_dict[d_id["contig_id"]].features[d_id["gene_number"]]
            utils.set_feature_qualifier(feature, "signal_peptide", words[1])
            count += 1
        utils.log(f'Signal peptide discovery complete. Found {count} proteins with signal peptide.')
