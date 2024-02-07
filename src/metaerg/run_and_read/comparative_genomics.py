from pathlib import Path
from collections import deque, Counter
from metaerg.datatypes.blast import TabularBlastParser
from metaerg import context
from metaerg.datatypes.fasta import FastaParser, write_fasta
from  metaerg.datatypes.sqlite import connect_to_db, read_feature_by_id, update_feature_in_db, read_all_features


def _run_programs():
    input_folder = context.BASE_DIR / 'faa'



def merge_and_code_fasta_input(mag_faa_dir: Path, mag_faa_file_extension: str, delimiter: str, taxa_by_orf_id: list,
                               merged_fasta_file: Path):
    context.log(f'Now merging and coding fasta input as {merged_fasta_file}..')
    unique_ids = set()
    with open(merged_fasta_file, 'w') as writer:
        file_count = 0
        orf_count = 0
        for file in sorted(mag_faa_dir.glob(f'*{mag_faa_file_extension}')):
            fasta_file = mag_faa_dir / file.name
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



class Cluster:
    def __init__(self, seq_ids: list, taxa_by_orf_id: list, blast_scores: dict):
        self.id = seq_ids[0]
        try:
            self.seq_ids = sorted(seq_ids, key=lambda k: blast_scores[self.id][k][0], reverse=True)
            self.fraction_id = sum([blast_scores[self.id][id2][1] for id2 in self.seq_ids[1:]]) / (len(self.seq_ids)-1)
            self.error = False
        except KeyError:
            context.log('Warning: no alignment information for cluster - could not compute percentage id and orthologues '
                  'could not be called for this cluster.')
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
            min_fraction_orthologues:float=0.5, min_fraction_of_taxa_represented:float=0.1,
            min_taxa_represented:float=3, min_fraction_of_genes_per_taxon:float=0.1,
            include_paralogues_in_fasta_output=True):
    unique_taxa = {taxon for taxon in taxa_by_orf_id}
    context.log(f'Detected {len(unique_taxa)} unique taxa.')
    # prep files
    tmp_dir.mkdir(exist_ok=True)
    cluster_file_base = input_fasta_file.parent / f'{input_fasta_file.stem}.clustering'
    cluster_file_tsv = input_fasta_file.parent / f'{input_fasta_file.stem}.clustering.tsv'
    cluster_file_align = input_fasta_file.parent / f'{input_fasta_file.stem}.clustering.align'
    cluster_file_blast = input_fasta_file.parent / f'{input_fasta_file.stem}.clustering.blast'
    input_fasta_db = input_fasta_file.parent / (input_fasta_file.stem + '.mmseqdb')
    # run programs
    context.run_external(f'mmseqs createdb {input_fasta_file} {input_fasta_db}')
    context.run_external(f'mmseqs cluster -c {fraction_overlap} --cov-mode 0 --min-seq-id {fraction_id} '
                 f'{input_fasta_db} {cluster_file_base} {tmp_dir}')
    context.run_external(f'mmseqs createtsv {input_fasta_db} {input_fasta_db} {cluster_file_base} {cluster_file_tsv}')
    context.run_external(f'mmseqs align {input_fasta_db} {input_fasta_db} {cluster_file_base} '
                 f'{cluster_file_align} -a')
    context.run_external(f'mmseqs convertalis {input_fasta_db} {input_fasta_db} {cluster_file_align} '
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


def parse_clusters(fasta_custered_dir: Path):
    # [1] create a tsv with cluster info
    cluster_size_hash = {}
    with open(fasta_custered_dir / 'homologues.tsv', 'w') as tsv_writer:
        tsv_writer.write('cluster id\tannotation\trepresentation\tcount' + '\t'.join(taxa) + '\t' + '\t'.join(taxa) + '\n')
        for cluster_fasta_file in sorted(fasta_custered_dir.glob('.faa')):
            cluster_id = int(cluster_fasta_file.stem)
            annotation = ''
            seq_hash = {}
            taxa_represented = set()
            total = 0
            with FastaParser(cluster_fasta_file, cleanup_seq=False) as fasta_reader:
                for seq in fasta_reader:
                    (taxon, feature_id) = seq['id'].split('~')
                    is_paralogue = '(PARALOGUE)' in seq['descr']
                    try:
                        seq_hash[taxon].append((feature_id, is_paralogue))
                    except KeyError:
                        seq_hash[taxon] = [(feature_id, is_paralogue)]
                    taxa_represented.add(taxon)
                    total += 1
                    if '(CENTER)' in seq['desrc']:
                        annotation = seq['descr'][9:]
            tsv_writer.write(f'{cluster_id}\t{annotation}\t{len(taxa_represented)}\t{total}\t' +
                              '\t'.join([str(len(seq_hash.get(t, []))) for t in taxa]) +
                              '\t'.join([','.join(seq_hash.get(t, [])) for t in taxa]) + '\n')


def run(cluster_window_size = 4):
    context.log('Mode "Clade": Clustering homologous genes across genomes...')
    input_dir = context.BASE_DIR / 'faa'
    comparative_genomics_dir = context.BASE_DIR / 'clade'
    merged_fasta_file = comparative_genomics_dir / 'all_cds_from_all_genomes_coded'
    taxa_by_orf_id = []
    taxa = sorted([f. name for f in input_dir.glob('*')])
    fasta_custered_dir = comparative_genomics_dir / 'clusters.faa'
    merge_and_code_fasta_input(input_dir, mag_faa_file_extension='', delimiter='~', taxa_by_orf_id=taxa_by_orf_id,
                               merged_fasta_file=merged_fasta_file)
    cluster(taxa_by_orf_id=taxa_by_orf_id,
            input_fasta_file=merged_fasta_file,
            tmp_dir=comparative_genomics_dir / 'tmp',
            fasta_output_dir=fasta_custered_dir,
            fraction_id=0.5)

            # [2] for each gene, store "prevalence" in SQL and cluster id
            for taxon, feature_id_list in seq_hash.items():
                db_file = context.BASE_DIR / 'annotations.sqlite' / taxon
                for feature_id, is_paralogue in feature_id_list:
                    db_connection = connect_to_db(db_file)
                    feature = read_feature_by_id(db_connection, feature_id)
                    feature.homologous_group_feature_is_paralogue = is_paralogue
                    feature.homologous_group_id = int(cluster_id)
                    feature.homologous_group_taxon_representation = len(taxa_represented) / len(taxa)
                    feature.homologous_group_member_count = total
                    update_feature_in_db(db_connection, feature_id)
            cluster_size_hash[cluster_id] = total

    # 3 determine contextual links between orthologues (from 2)
    profile_match_counts = Counter()
    for taxon in taxa:
        window = deque(maxlen=cluster_window_size)
        db_file = context.BASE_DIR / 'annotations.sqlite' / taxon
        db_connection = connect_to_db(db_file, target='Features')
        previous_contig = None
        for feature in read_all_features(db_connection, type=('CDS', 'ncRNA')):
            # make sure subsequent features are actually found together
            if previous_contig and feature.contig != previous_contig:
                window.clear()
            if not feature.homologous_group_id:
                continue
            # here comes the logic
            for other_feature in window:
                if (not other_feature.homologous_group_id or
                        other_feature.homologous_group_id == feature.homologous_group_id):
                    continue
                profile_match_counts.update(make_match_key(feature.homologous_group_id,
                                                           other_feature.homologous_group_id))
            # post-processing...
            previous_contig = feature.contig
            window.append(feature)


def make_match_key(id1: int, id2: int):
    if id1 > id2:
        return f'{id1}_{id2}'
    else:
        return f'{id2}_{id1}'
