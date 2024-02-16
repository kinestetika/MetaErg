import re
from statistics import median, quantiles
from pathlib import Path
from collections import deque, Counter
from metaerg.datatypes.blast import TabularBlastParser
from metaerg import context
from metaerg.datatypes.fasta import FastaParser, write_fasta
from  metaerg.datatypes import sqlite
from metaerg import registry
from metaerg.html import html_feature_table, html_feature_details


def _run_programs():
    input_folder = context.BASE_DIR / 'faa'



def merge_and_code_fasta_input(mag_faa_dir: Path, mag_faa_file_extension: str, delimiter: str, taxa_by_orf_id: list,
                               merged_fasta_file: Path) -> list:
    context.log(f'Now merging and coding fasta input as {merged_fasta_file}..')
    taxa = list()
    file_count = 0
    orf_count = 0
    with open(merged_fasta_file, 'w') as writer:
        for file in sorted(mag_faa_dir.glob(f'*{mag_faa_file_extension}')):
            unique_ids = set()
            fasta_file = mag_faa_dir / file.name
            taxa.append(file.stem)
            with FastaParser(fasta_file, cleanup_seq=False) as fasta_reader:
                for orf in fasta_reader:
                    if orf['id'] in unique_ids:
                        context.log(f'warning: duplicate id for {orf["id"]} in {file.name}: skipping.')
                        continue
                    unique_ids.add(orf['id'])
                    if orf['seq'][-1] == '*':
                        orf['seq'] = orf['seq'][:-1]
                    recoded_orf = {'id': f'{orf_count}',
                                   'seq': orf['seq'],
                                   'descr': f'{file.stem}{delimiter}{orf["id"]} {orf["descr"]}'.strip()
                                  }
                    write_fasta(writer, recoded_orf)
                    taxa_by_orf_id.append(file_count)
                    orf_count += 1
            file_count += 1
    context.log(f'Merging and coding fasta input complete; recoded and wrote {len(taxa_by_orf_id)} proteins to fasta, '
        f'{file_count} unique taxa.')
    return taxa


class Cluster:
    def __init__(self, seq_ids: list, taxa_by_orf_id: list, blast_scores: dict):
        self.id = seq_ids[0]
        try:
            self.seq_ids = sorted(seq_ids, key=lambda k: blast_scores[self.id][k][0], reverse=True)
            self.fraction_id = sum([blast_scores[self.id][id2][1] for id2 in self.seq_ids[1:]]) / (len(self.seq_ids)-1)
            self.error = False
        except KeyError:
            context.log('Warning: no alignment information for cluster - could not compute percentage id and orthologues '
                  'could not be called for cluster.')
            self.error = True
            self.seq_ids = seq_ids
            self.fraction_id = 0
        self.paralogues = {}
        self.taxa = set()
        for seq_id in self.seq_ids:
            self.paralogues[seq_id] = taxa_by_orf_id[seq_id] in self.taxa
            self.taxa.add(taxa_by_orf_id[seq_id])
        self.fraction_orthologues = len(self.taxa) / len(self.seq_ids)
        # print(f'cluster at {self.fraction_id:.1%}')
        # for i in range(len(self.seq_ids)):
        #    print(self.seq_ids[i], taxa_by_orf_id[self.seq_ids[i]], blast_scores[self.id][self.seq_ids[i]], self.paralogues[i])


class MMSeqsClusterParser:
    # This class parses a *_cluster.tsv file created by "mmseqs ..."
    def __init__(self, path, taxa_by_orf_id: list, blast_scores: dict):
        self.path = path
        self.handle = None
        self.taxa_by_orf_id = taxa_by_orf_id
        self.blast_scores = blast_scores

    def __enter__(self):
        self.handle = open(self.path)
        return self

    def __exit__(self, exception_type, exception_value, exception_traceback):
        self.handle.close()

    def __iter__(self):
        current_cluster_seq_id_list = []
        while line := self.handle.readline():
            id1, id2 = line.strip().split('\t')
            if id1 == id2:
                if len(current_cluster_seq_id_list) > 1:
                    yield Cluster(current_cluster_seq_id_list, self.taxa_by_orf_id, self.blast_scores)
                current_cluster_seq_id_list = [int(id1)]
            elif len(current_cluster_seq_id_list):
                current_cluster_seq_id_list.append(int(id2))
        if len(current_cluster_seq_id_list) > 1:
            yield Cluster(current_cluster_seq_id_list, self.taxa_by_orf_id, self.blast_scores)


def cluster(taxa_by_orf_id: list, input_fasta_file:Path, tmp_dir: Path, fasta_output_dir: Path,
            fraction_id:float=0.5, fraction_overlap:float=0.8,
            min_fraction_orthologues:float=0.0, min_fraction_of_taxa_represented:float=0.1,
            min_taxa_represented:float=3, min_fraction_of_genes_per_taxon:float=0.1,
            include_paralogues_in_fasta_output=True):
    unique_taxa = {taxon for taxon in taxa_by_orf_id}
    context.log('Now clustering CDS genes...')
    context.log(f'Detected {len(unique_taxa)} unique taxa.')
    # prep files
    tmp_dir.mkdir(exist_ok=True)
    cluster_file_base = input_fasta_file.parent / f'{input_fasta_file.stem}.clustering'
    cluster_file_tsv = input_fasta_file.parent / f'{input_fasta_file.stem}.clustering.tsv'
    cluster_file_align = input_fasta_file.parent / f'{input_fasta_file.stem}.clustering.align'
    cluster_file_blast = input_fasta_file.parent / f'{input_fasta_file.stem}.clustering.blast'
    input_fasta_db = input_fasta_file.parent / (input_fasta_file.stem + '.mmseqdb')
    # run programs
    mmseqs_path = ''  #'/bio/bin/mmseqs/bin/'  # for testing
    context.run_external(f'{mmseqs_path}mmseqs createdb {input_fasta_file} {input_fasta_db}')
    context.run_external(f'{mmseqs_path}mmseqs cluster -c {fraction_overlap} --cov-mode 0 --min-seq-id {fraction_id} '
                 f'{input_fasta_db} {cluster_file_base} {tmp_dir}')
    context.run_external(f'{mmseqs_path}mmseqs createtsv {input_fasta_db} {input_fasta_db} {cluster_file_base} {cluster_file_tsv}')
    context.run_external(f'{mmseqs_path}mmseqs align {input_fasta_db} {input_fasta_db} {cluster_file_base} '
                 f'{cluster_file_align} -a')
    context.run_external(f'{mmseqs_path}mmseqs convertalis {input_fasta_db} {input_fasta_db} {cluster_file_align} '
                 f'{cluster_file_blast}')
    # parse results
    blast_scores = {}
    with TabularBlastParser(cluster_file_blast, 'BLAST') as reader:
        for r in reader:
            if len(r.hits) > 1:
                blast_scores[int(r.hits[0].query)] = {int(h.hit): (h.score, h.percent_id) for h in r.hits}
    seq_id2cluster = {}
    seq_id2cluster_rejected = {}
    rejected_cluster_count = 0
    error_count = 0
    cluster_count = 0
    percent_id = 0
    taxon_representaton = Counter()
    with MMSeqsClusterParser(cluster_file_tsv, taxa_by_orf_id, blast_scores) as reader:
        for cluster in reader:
            if cluster.fraction_orthologues < min_fraction_orthologues \
                    or len(cluster.taxa) < max(min_taxa_represented, min_fraction_of_taxa_represented * len(unique_taxa))\
                    or cluster.error:
                rejected_cluster_count += 1
                error_count += cluster.error
                for seq_id in cluster.seq_ids:
                    seq_id2cluster_rejected[seq_id] = cluster
                cluster.id = f'rejected_{rejected_cluster_count}'
            else:
                cluster.id = cluster_count + 1  # we do  not want a cluster with id zero
                cluster_count += 1
                percent_id += cluster.fraction_id
                for seq_id in cluster.seq_ids:
                    seq_id2cluster[seq_id] = cluster
                taxon_representaton.update(cluster.taxa)

    context.log(f'{len(seq_id2cluster)}/{len(taxa_by_orf_id)} seqs clustered ({len(seq_id2cluster)/len(taxa_by_orf_id):.1%}).')
    context.log(f'Accepted {cluster_count} clusters, rejected {rejected_cluster_count}, {error_count} due to errors.')
    context.log(f'Estimated average percent id: {percent_id/cluster_count:.1%}')
    # filter out taxa with poor representaton
    poorly_represented_taxa = {t for t in unique_taxa if taxon_representaton[t] < min_fraction_of_genes_per_taxon * cluster_count}
    for t in poorly_represented_taxa:
        context.log(f'Rejected taxon {t}; represented in only {taxon_representaton[t]/cluster_count:.1%} of clusters.')

    # write output
    fasta_output_dir.mkdir(exist_ok=True)
    rejected_fasta_output_dir = fasta_output_dir / 'rejected'
    rejected_fasta_output_dir.mkdir(exist_ok=True)
    for file in rejected_fasta_output_dir.glob('*'):
        if not file.is_dir():
            file.unlink()
    for file in fasta_output_dir.glob('*'):
        if not file.is_dir():
            file.unlink()
    i = 0
    with FastaParser(input_fasta_file, cleanup_seq=False) as fasta_reader:
        for orf in fasta_reader:
            if taxa_by_orf_id[int(orf['id'])] in poorly_represented_taxa:
                continue
            if int(orf['id']) in seq_id2cluster.keys():
                cluster = seq_id2cluster[int(orf['id'])]
                target_fasta_file = fasta_output_dir / f'{cluster.id}.faa'
            elif int(orf['id']) in seq_id2cluster_rejected.keys():
                cluster = seq_id2cluster_rejected[int(orf['id'])]
                target_fasta_file = rejected_fasta_output_dir / f'{cluster.id}.faa'
            else:
                continue
            if include_paralogues_in_fasta_output or not cluster.paralogues[int(orf['id'])]:
                paralogue_text = '(PARALOGUE) ' if cluster.paralogues[int(orf['id'])] else ''
                center_text = '(CENTER) ' if cluster.seq_ids[0] == int(orf['id']) else ''
                with open(target_fasta_file, 'a') as writer:
                    try:
                        descr_space_index = orf['descr'].index(' ')
                        orf['id'] = orf['descr'][:descr_space_index]
                        orf['descr'] = center_text + paralogue_text + orf['descr'][descr_space_index+1:]
                    except ValueError:
                        orf['id'] = orf['descr']
                        orf['descr'] = center_text + paralogue_text
                    write_fasta(writer, orf)
                    i += 1
    context.log(f'wrote {i} seqs to files in {fasta_output_dir}')


def parse_clusters(fasta_custered_dir: Path) -> dict:
    context.log(f"Now parsing clusters of homologous genes from dir '{fasta_custered_dir}'...")
    cluster_hash = {}
    taxonomy_pattern = re.compile(r'\s\(Bacteria;.+?\)$|\s\(Archaea;.+?\)$|\s\(Viruses;.+?\)$|\s\(Eukaryota;.+?\)$') # |(\s\(Archaea;.+?\)$)|(\s\(Viruses;.+?\)$)
    for cluster_fasta_file in sorted(fasta_custered_dir.glob('*.faa')):
        cluster = {'id':         int(cluster_fasta_file.stem),
                   'annotation': '',
                   'members': [],  # list of tuples (id, taxon, is_paralogue]
                   'taxa': set(),
                   'length': 0,
                   'links': [],  # links to other clusters based on gene-context, list of tuples (other cluster, score)
                   'fraction_id': 0,  # how close the members of the cluster are to each other (median of orthologous members)
                   'codon_bias_quantile': -1,  # quantile of median codon usage bias (orthologous members only), 0 = low expression, 9 = high expression
                   'selection_omega_quantile': -1  # quantile of median w = ka / ks (orthologous members only), 0 = purifying, 9 = diversifying
                   }
        lengths = []
        with FastaParser(cluster_fasta_file, cleanup_seq=False) as fasta_reader:
            for seq in fasta_reader:
                (taxon, feature_id) = seq['id'].split('~')
                is_paralogue = '(PARALOGUE)' in seq['descr']
                cluster['members'].append((feature_id, taxon, is_paralogue))
                cluster['taxa'].add(taxon)
                if '(CENTER)' in seq['descr']:
                    cluster['annotation'] = seq['descr'][9:]
                cluster['annotation'] = re.sub(taxonomy_pattern, '', cluster['annotation'])
                lengths.append(len(seq['seq']))
        cluster_hash[cluster['id']] = cluster
        cluster['length'] = median(lengths)
    context.log(f"Parsed {len(cluster_hash)} clusters of homologous genes.")
    return cluster_hash


def save_clusters_as_nt_fasta_and_compute_codon_bias(clusterhash: dict, all_taxa: list, nt_fasta_output_dir, delimiter='~'):
    # also calculates median codon usage biases for orthologues
    context.log(f"Now writing clusters of homologous genes as fasta nt to '{nt_fasta_output_dir}':")
    ids_by_taxon = {t: [] for t in all_taxa}
    # reorganize data for quicker loading
    for cluster in clusterhash.values():
        for id, taxon, is_paralogue in cluster['members']:
            if not is_paralogue:
                ids_by_taxon[taxon].append(id)
    nt_seq_hash = {}
    codon_usage_sequence_hash = {}
    # retrieve nt seqs from sql db
    context.log(f"First reading nt sequences from database...")
    for i, taxon in zip(range(len(ids_by_taxon)), ids_by_taxon.keys()):
        context.log(f'Extracting data for ({i+1}/{len(all_taxa)}) "{taxon}" from sql...')
        db_file = context.BASE_DIR / 'annotations.sqlite' / taxon
        db_connection = sqlite.connect_to_db(db_file)
        for id in ids_by_taxon[taxon]:
            feature = sqlite.read_feature_by_id(db_connection, id)
            nt_seq_hash[taxon + delimiter + id] = {'id': taxon + '~' + id,
                                                   'descr': feature.descr,
                                                   'seq': feature.nt_seq
                                                  }
            if feature.codon_bias > 0:
                codon_usage_sequence_hash[taxon + delimiter + id] = feature.codon_bias
    # write files
    context.log(f"Now writing files...")
    nt_fasta_output_dir.mkdir(exist_ok=True)
    codon_bias_distribution = []
    for cluster in clusterhash.values():
        codon_usage_biases = []
        target_fasta_file = nt_fasta_output_dir / f'{cluster["id"]}.fna'
        with open(target_fasta_file, 'w') as writer:
            for id, taxon, is_paralogue in cluster['members']:
                if not is_paralogue:
                    write_fasta(writer, nt_seq_hash[taxon + delimiter + id])
                    try:
                        codon_usage_biases.append(codon_usage_sequence_hash[taxon + delimiter + id])
                    except KeyError:
                        pass
        if len(codon_usage_biases):
            cluster['codon_bias'] = median(codon_usage_biases)
            codon_bias_distribution.append(cluster['codon_bias'])
    # aggregate codon biases, set the quantile for each cluster
    codon_bias_quantiles = [round(q, 3) for q in quantiles(codon_bias_distribution, n=10)]
    for cluster in clusterhash.values():
        try:
            for i, q in zip(range(len(codon_bias_quantiles)), codon_bias_quantiles):
                if 0 < cluster['codon_bias'] < q:
                    cluster['codon_bias_quantile'] = i
                    #print(cluster['codon_bias'], cluster['codon_bias_quantile'])
                    break
        except:
            pass



def align_aa_clusters(fasta_custered_dir: Path, alignment_dir: Path):
    context.log(f"Aligning clusters of homologous genes to '{alignment_dir}'...")
    alignment_dir.mkdir(exist_ok=True)
    famsa_path = '' #'/bio/bin/'   # only for testing
    for f in sorted(fasta_custered_dir.glob('*.faa')):
        aligned_f = alignment_dir / f.name
        context.run_external(f'{famsa_path}famsa {f} {aligned_f}')


def estimate_selection_pressure(cluster_hash: dict, aa_align_dir: Path, nt_fasta_dir: Path, delimiter = '~', minimum_mutations = 25):
    # also calculates average percent id (aminoacids) among orthologous cluster-members
    context.log(f"Estimating selective pressure for each cluster (w = Ka/Ks)... ")
    omega_distribution = []
    fractions_id = []
    for cluster in cluster_hash.values():
        # read aa alignment
        aa_align_file = aa_align_dir / (str(cluster['id']) + '.faa')
        aa_seq_hash = {}
        with FastaParser(aa_align_file, cleanup_seq=False) as fasta_reader:
            for fasta_seq in fasta_reader:
                aa_seq_hash[fasta_seq["id"]] = fasta_seq
        # read associated nt seqs
        nt_file = nt_fasta_dir / (str(cluster['id']) + '.fna')
        nt_seq_hash = {}
        with FastaParser(nt_file, cleanup_seq=False) as fasta_reader:
            for fasta_seq in fasta_reader:
                nt_seq_hash[fasta_seq["id"]] = fasta_seq
        # now analyze
        synonymous_mutations = []  # ks
        non_synonymous_mutations = []  # ka  w = ka / ks
        for id1, taxon1, is_paralogue1 in cluster["members"]:
            if is_paralogue1:
                continue
            seq_id1 = taxon1 + delimiter + id1
            aa_seq1 = aa_seq_hash[seq_id1]['seq']
            nt_seq1 = nt_seq_hash[seq_id1]['seq']
            for id2, taxon2, is_paralogue2 in cluster["members"]:
                if is_paralogue2:
                    continue
                if id1 == id2:
                    continue
                seq_id2 = taxon2 + delimiter + id2
                aa_seq2 = aa_seq_hash[seq_id2]['seq']
                nt_seq2 = nt_seq_hash[seq_id2]['seq']
                aa_seq_pos1 = 0
                aa_seq_pos2 = 0
                sm = 0
                nsm = 0
                for aligned_pos in range(len(aa_seq1)):
                    aa1 = aa_seq1[aligned_pos]
                    aa2 = aa_seq2[aligned_pos]
                    if aa1 == 'X' or aa2 == 'X':
                        aa_seq_pos1 += 1
                        aa_seq_pos2 += 1
                    elif aa1 == '-' and aa2 == '-':
                        pass
                    elif aa1 == '-':
                        aa_seq_pos2 += 1
                        nsm += 1
                    elif aa2 == '-':
                        aa_seq_pos1 += 1
                        nsm += 1
                    else:  # no gaps, no X
                        codon1 = nt_seq1[aa_seq_pos1 * 3:aa_seq_pos1 * 3 + 3]
                        codon2 = nt_seq2[aa_seq_pos2 * 3:aa_seq_pos2 * 3 + 3]
                        if aa1 != aa2:
                            nsm += 1
                        elif codon1 != codon2:
                            sm += 1
                        aa_seq_pos1 += 1
                        aa_seq_pos2 += 1
                synonymous_mutations.append(sm)
                non_synonymous_mutations.append(nsm)
                fractions_id.append((aa_seq_pos1 - nsm) / aa_seq_pos1)
        smm = median(synonymous_mutations)
        nsmm = median(non_synonymous_mutations)
        if smm and smm + nsmm > minimum_mutations:
            cluster['synonymous_mutations_median'] = smm
            cluster['non_synonymous_mutations_median'] = nsmm
            cluster['selection_omega'] = nsmm / smm
            omega_distribution.append(nsmm / smm)
        cluster['fraction_id'] = sum(fractions_id) / len(fractions_id)  # average
    # aggregate
    omega_quantiles = [round(q, 3) for q in quantiles(omega_distribution, n=10)]
    for cluster in sorted([c for c in cluster_hash.values()], key=lambda i:i.get('selection_omega', 100)):
        try:
            for i, q in zip(range(len(omega_quantiles)), omega_quantiles):
                if 0 < cluster['selection_omega'] < q:
                    cluster['selection_omega_quantile'] = i
                    break
    #        print(f'{cluster["selection_omega"]}\t{len(cluster["taxa"])}\t{len(cluster["members"])}\t{cluster["annotation"]}')
        except:
            pass


def write_overview_tsv_file(dir: Path, cluster_hash: dict, all_taxa: list):
    context.log("Now writing homology overview info to files 'homologues.tsv'...")
    # each cluster on one line, we write id<tab>annotation<tab>#taxa<tab>#seqs
    # <tab>gene-context associations to other clusters with scores
    # then (separated by <tab>) for each taxon the number of reps,
    # then (separated by <tab>) for each taxon the ids of all reps (separated by commas)
    with open(dir / 'homologues.tsv', 'w') as tsv_writer:
        tsv_writer.write('cluster id\tannotation\tlength\tselective pressure\tcodon bias\trepresentation\tcount\t% id\tlinks\t' + '\t'.join(all_taxa) + '\t' + '\t'.join(all_taxa) + '\n')
        for cluster in sorted(cluster_hash.values(), key=lambda c:len(c['taxa']), reverse=True):  # sort by representation
            tsv_writer.write(f'{cluster["id"]}\t'
                             f'{cluster["annotation"]}\t'
                             f'{cluster["length"]:.0f}\t'
                             f'{cluster["selection_omega_quantile"] if cluster["selection_omega_quantile"] >= 0 else ""}\t'
                             f'{cluster["codon_bias_quantile"] if cluster["codon_bias_quantile"] >= 0 else ""}\t'
                             f'{len(cluster["taxa"])}\t'
                             f'{len(cluster["members"])}\t'
                             f'{cluster["fraction_id"]*100:.1f}\t' +
                              ', '.join([f'{l[0]["id"]}:{l[1]:.1f}' for l in cluster["links"]]) + '\t' +
                              '\t'.join([str(sum([1 for m in cluster['members'] if m[1] == t])) for t in all_taxa]) + '\t' +
                              '\t'.join([','.join([m[0] for m in cluster['members'] if m[1] == t]) for t in all_taxa]) + '\n')


def write_homologue_info_to_sql(cluster_hash: dict, all_taxa: list):
    context.log("Now writing info about homology to each genome's feature database...")
    for i, taxon in zip(range(len(all_taxa)), all_taxa):
        targets = []
        for c in cluster_hash.values():  # round up all the features for a given taxon
            targets.extend([(m[0], m[1], m[2], c) for m in c['members'] if m[1] == taxon])
        context.log(f'Processing ({i+1}/{len(all_taxa)}) taxa: "{taxon}" has {len(targets)} protein-coding genes with homologues in other species.')
        db_file = context.BASE_DIR / 'annotations.sqlite' / taxon
        db_connection = sqlite.connect_to_db(db_file)
        for feature_id, taxon, is_paralogue, cluster in targets:
            feature = sqlite.read_feature_by_id(db_connection, feature_id)
            feature.homologous_group_feature_is_paralogue = is_paralogue
            feature.homologous_group_id = int(cluster['id'])
            feature.homologous_group_taxon_representation = len(cluster['taxa']) / len(all_taxa)
            feature.homologous_group_member_count = len(cluster['members'])
            sqlite.update_feature_in_db(db_connection, feature)


def analyse_gene_context(cluster_hash, all_taxa: list, window_size: int = 4, min_match_score: float = 0.5) -> dict:
    context.log("Now analysing gene context to identify homologous genes frequently found together...")
    context.log(f"Working with window size {window_size} and requiring {min_match_score:.1%} of pairs of homologous genes to be found close to each other.")
    profile_match_counts = Counter()
    for i, taxon in zip(range(len(all_taxa)), all_taxa):
        context.log(f'Processing ({i+1}/{len(all_taxa)}) taxa: "{taxon}"')
        window = deque(maxlen=window_size)
        db_file = context.BASE_DIR / 'annotations.sqlite' / taxon
        db_connection = sqlite.connect_to_db(db_file, target='Features')
        previous_contig = None
        for feature in sqlite.read_all_features(db_connection, type=('CDS')):
            # make sure subsequent features are actually on the same contig
            if previous_contig and feature.contig != previous_contig:
                window.clear()
            elif not feature.homologous_group_id:
                previous_contig = feature.contig
                window.append(feature)
                continue
            # here comes the logic
            for other_feature in window:
                if not other_feature.homologous_group_id or \
                        other_feature.homologous_group_id == feature.homologous_group_id:
                    continue
                profile_match_counts.update([make_match_key(feature.homologous_group_id,
                                                           other_feature.homologous_group_id)])
            # post-processing...
            previous_contig = feature.contig
            window.append(feature)
    # aggregate and update clusters
    cluster_context_links = dict()
    for match, match_count in profile_match_counts.most_common():
        (cluster_id_1, cluster_id_2) = match.split()
        cluster_id_1 = int(cluster_id_1)
        cluster_id_2 = int(cluster_id_2)
        match_score = match_count / min(len(cluster_hash[cluster_id_1]['members']),
                                        len(cluster_hash[cluster_id_2]['members']))
        if match_score > min_match_score:
            cluster_context_links[match] = match_score
            cluster_hash[cluster_id_1]['links'].append((cluster_hash[cluster_id_2], match_score))
            cluster_hash[cluster_id_2]['links'].append((cluster_hash[cluster_id_1], match_score))
    context.log(f'Identified {len(cluster_context_links)}/{len(profile_match_counts)} gene-context associations '
                f'stronger than {min_match_score:.2}')
    return cluster_context_links


def write_remaining_results_to_sql(cluster_hash: dict, cluster_context_links: dict, all_taxa: list, window_size: int = 4):
    context.log("Now computing and saving gene context to each genome's feature database...")

    for taxon in all_taxa:
        db_file = context.BASE_DIR / 'annotations.sqlite' / taxon
        db_connection = sqlite.connect_to_db(db_file)
        window = deque(maxlen=window_size)
        # only for testing
        #for feature in sqlite.read_all_features(db_connection, type=('region')):
        #    sqlite.drop_feature(db_connection, feature)
        existing_regions = list()
        for feature in sqlite.read_all_features(db_connection, type=('CDS')):
            try:
                cluster = cluster_hash[feature.homologous_group_id]
                feature.homologous_group_codon_usage_bias_quantile = cluster["codon_bias_quantile"]
                feature.homologous_group_selective_pressure_quantile = cluster["selection_omega_quantile"]
                sqlite.update_feature_in_db(db_connection, feature)
            except KeyError:
                pass
            #print(f'@{feature.id}, HGID: {feature.homologous_group_id}')
            if not len(window):
                window.append(feature)
            elif feature.contig != window[-1].contig:  # new contig -> reset
                window.clear()
                window.append(feature)
            else:
                linked_feature = None
                for prev_feature in window:
                    if cluster_context_links.get(make_match_key(feature.homologous_group_id,
                                                            prev_feature.homologous_group_id), 0):
                        linked_feature = prev_feature
                        break
                if linked_feature:
                    #print(f'[{linked_feature.id} ~ {feature.id}')
                    target_region = None
                    for existing_region in existing_regions:
                        if linked_feature in existing_region:
                            target_region = existing_region
                            break
                    if not target_region:
                        target_region = [window[pos] for pos in range(window.index(linked_feature), len(window))]
                        existing_regions.append(target_region)
                    target_region.append(feature)
                window.append(feature)
        for region_id, region in zip(range(len(existing_regions)), existing_regions):
            region_feature = sqlite.Feature(genome=region[0].genome,
                                            contig=region[0].contig,
                                            start=max(region[0].start, 0),
                                            end=region[-1].end,
                                            strand=1,
                                            type='region',
                                            inference='metaerg',
                                            descr=f'gene cluster based on homology',
                                            id=f'region_{region_id}')
            sqlite.add_new_feature_to_db(db_connection, region_feature)
            for f in sqlite.read_all_features(db_connection, contig=region_feature.contig,
                                              start=region_feature.start, end=region_feature.end):
                if f.type != 'region':
                    f.parent = region_feature.id
                    sqlite.update_feature_in_db(db_connection, f)
        context.log(f"Wrote {len(existing_regions)} homology-informed sequence regions to sql for '{taxon}'")


def run(genome_dict: dict, cluster_window_size: int = 4, min_match_score: float = 0.5):
    context.log('Mode "comparative_genomics": Clustering and analysing homologous genes across genomes...')
    # prep/clear dirs
    input_dir = context.BASE_DIR / 'faa'
    comparative_genomics_dir = context.BASE_DIR / 'comparative_genomics'
    comparative_genomics_dir.mkdir(exist_ok=True)
    for f in comparative_genomics_dir.glob('*'):
        if f.is_file():
            f.unlink()
    db_dir = comparative_genomics_dir / 'mmseqs_db'
    db_dir.mkdir(exist_ok=True)
    for f in db_dir.glob('*'):
        if f.is_file():
            f.unlink()
    merged_fasta_file = db_dir / 'all_cds_from_all_genomes_coded.faa'
    fasta_custered_dir = comparative_genomics_dir / 'clusters.faa'
    fasta_custered_dir.mkdir(exist_ok=True)
    for f in fasta_custered_dir.glob('*'):
        if f.is_file():
            f.unlink()
    nt_fasta_custered_dir = comparative_genomics_dir / 'clusters.fna'
    nt_fasta_custered_dir.mkdir(exist_ok=True)
    for f in nt_fasta_custered_dir.glob('*'):
        if f.is_file():
            f.unlink()
    alignment_dir = comparative_genomics_dir / 'clusters.faa.align'
    alignment_dir.mkdir(exist_ok=True)
    for f in alignment_dir.glob('*'):
        if f.is_file():
            f.unlink()
    #html_dir = comparative_genomics_dir / 'html'

    taxa_by_orf_id = []
    taxa = merge_and_code_fasta_input(input_dir, mag_faa_file_extension='', delimiter='~', taxa_by_orf_id=taxa_by_orf_id,
                                      merged_fasta_file=merged_fasta_file)
    cluster(taxa_by_orf_id=taxa_by_orf_id,
            input_fasta_file=merged_fasta_file,
            tmp_dir=comparative_genomics_dir / 'tmp',
            fasta_output_dir=fasta_custered_dir,
            fraction_id=0.5)
    cluster_hash = parse_clusters(fasta_custered_dir)
    align_aa_clusters(fasta_custered_dir, alignment_dir)
    save_clusters_as_nt_fasta_and_compute_codon_bias(cluster_hash, taxa, nt_fasta_custered_dir)
    estimate_selection_pressure(cluster_hash, alignment_dir, nt_fasta_custered_dir)
    write_homologue_info_to_sql(cluster_hash, taxa)
    cluster_context_links = analyse_gene_context(cluster_hash, taxa, cluster_window_size, min_match_score)
    write_overview_tsv_file(comparative_genomics_dir, cluster_hash, taxa)
    write_remaining_results_to_sql(cluster_hash, cluster_context_links, taxa, cluster_window_size)
    context.log('Updating html visualizations for all genomes...')
    for genome_name, genome in genome_dict.items():
        context.log(f'({genome_name}) Now writing final result as .html for visualization...')
        feature_db_connection = sqlite.connect_to_db(context.BASE_DIR / 'annotations.sqlite' / genome_name)
        for html_writer in registry.HTML_WRITER_REGISTRY:
            html_writer(genome, feature_db_connection, context.HTML_DIR)


def make_match_key(id1: int, id2: int):
    if id1 > id2:
        return f'{id1} {id2}'
    else:
        return f'{id2} {id1}'


def main():
    # for testing...
    # may need to update mmseqs_path in 'cluster'
    context.BASE_DIR = Path('/bio/fast/alkaline_origins/comparative_genomics_metaerg_test')
    context.LOG_FILE = context.BASE_DIR / 'log_test_comparative_genomics'

    context.log('Mode "Clade": Clustering and analysing homologous genes across genomes...')
    # prep/clear dirs
    input_dir = context.BASE_DIR / 'faa'
    comparative_genomics_dir = context.BASE_DIR / 'comparative_genomics'
    comparative_genomics_dir.mkdir(exist_ok=True)
    #for f in comparative_genomics_dir.glob('*'):
    #    if f.is_file():
    #        f.unlink()
    db_dir = comparative_genomics_dir / 'mmseqs_db'
    db_dir.mkdir(exist_ok=True)
    #for f in db_dir.glob('*'):
    #    if f.is_file():
    #        f.unlink()
    merged_fasta_file = db_dir / 'all_cds_from_all_genomes_coded.faa'
    fasta_custered_dir = comparative_genomics_dir / 'clusters.faa'
    fasta_custered_dir.mkdir(exist_ok=True)
    #for f in fasta_custered_dir.glob('*'):
    #    if f.is_file():
    #        f.unlink()
    html_dir = comparative_genomics_dir / 'html'
    nt_fasta_custered_dir = comparative_genomics_dir / 'clusters.fna'
    nt_fasta_custered_dir.mkdir(exist_ok=True)
    #for f in fasta_custered_dir.glob('*'):
    #    if f.is_file():
    #        f.unlink()
    html_dir = comparative_genomics_dir / 'html'
    alignment_dir = comparative_genomics_dir / 'clusters.faa.align'

    taxa_by_orf_id = []
    taxa = merge_and_code_fasta_input(input_dir, mag_faa_file_extension='', delimiter='~', taxa_by_orf_id=taxa_by_orf_id,
                                      merged_fasta_file=merged_fasta_file)
    #cluster(taxa_by_orf_id=taxa_by_orf_id,
    #        input_fasta_file=merged_fasta_file,
    #        tmp_dir=comparative_genomics_dir / 'tmp',
    #        fasta_output_dir=fasta_custered_dir,
    #        fraction_id=0.5)
    #cluster_hash = parse_clusters(fasta_custered_dir)
    #align_aa_clusters(fasta_custered_dir, alignment_dir)
    #save_clusters_as_nt_fasta_and_compute_codon_bias(cluster_hash, taxa, nt_fasta_custered_dir)
    #estimate_selection_pressure(cluster_hash, alignment_dir, nt_fasta_custered_dir)
    #write_homologue_info_to_sql(cluster_hash, taxa)
    #print('\n'.join(sorted(cluster_hash.keys())))
    #cluster_context_links = analyse_gene_context(cluster_hash, taxa)
    #write_overview_tsv_file(comparative_genomics_dir, cluster_hash, taxa)
    #taxa = taxa[:1]
    #write_gene_context_to_sql(cluster_context_links, taxa)


    context.log('Updating html visualizations for all genomes...')
    genome_db = sqlite.connect_to_db(context.BASE_DIR / 'genome_properties.sqlite', 'Genomes')
    for taxon in taxa:
        genome_name = taxon
        print(genome_name)
        genome = sqlite.read_genome_by_id(genome_db, genome_name)
        context.log(f'({genome_name}) Now writing final result as .html for visualization... {html_dir}')
        feature_db_connection = sqlite.connect_to_db(context.BASE_DIR / 'annotations.sqlite' / genome_name)
        html_feature_table.write_html(genome, feature_db_connection, html_dir)
        #html_feature_details.write_html(genome, feature_db_connection, html_dir)
        break


if __name__ == "__main__":
    main()
