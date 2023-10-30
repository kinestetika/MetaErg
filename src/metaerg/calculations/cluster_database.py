from pathlib import Path
from collections import deque, Counter

from metaerg.datatypes import sqlite
from metaerg.datatypes.blast import BlastHit

CLUSTER_WINDOW_SIZE = 4
MAX_CDD_HITS = 5
MIN_MATCH_SCORE = 0.33

def compute_cluster_db(data_dir: Path, db_file: Path, data_file_extension: str = '.sqlite',
                       cluster_window_size: int = CLUSTER_WINDOW_SIZE, max_cdd_hits: int = MAX_CDD_HITS,
                       min_match_score:float = MIN_MATCH_SCORE):
    if cluster_window_size < 2:
        raise Exception(f'{cluster_window_size} < 2, should be >=2')

    profile_counts = Counter()
    profile_match_counts = Counter()
    if db_file.exists():
        with open(db_file) as prev_clustering_results:
            for line in prev_clustering_results:
                words = line.split('\t')
                (accession1, accession2) = words[0].split(' ')
                profile_counts[accession1] = int(words[3])
                profile_counts[accession2] = int(words[4])
                profile_match_counts[words[0]] = int(words[2])

    contig_files = [f for f in sorted(data_dir.glob(f'*{data_file_extension}')) if f.is_file()]
    for data_db_file in contig_files:
        window = deque(maxlen=cluster_window_size)
        db_connection = sqlite.connect_to_db(data_db_file, target='Features')
        previous_contig = None
        for feature in sqlite.read_all_features(db_connection, type=('CDS', 'ncRNA')):
            # make sure the features are actually found together
            if previous_contig and feature.contig != previous_contig:
                window.clear()
            # here comes the logic
            profile_counts.update((h.hit.accession for h in feature.cdd[:max_cdd_hits]))
            for cdd_hit_1 in feature.cdd[:max_cdd_hits]:  # a list of blast hits with a DBEntry in the 'hit' field
                profiles_done = set()  # we only want to create a match for each profile once
                for other_feature in window:
                    for cdd_hit_2 in other_feature.cdd[:max_cdd_hits]:
                        if cdd_hit_2.hit.accession in profiles_done:
                            continue
                        profile_match_counts.update((get_match_key(cdd_hit_1, cdd_hit_2),))
                        profiles_done.add(cdd_hit_2.hit.accession)
            # post-processing...
            previous_contig = feature.contig
            window.append(feature)
    # aggregate and write results
    count = 0
    with open(db_file, 'w') as writer:
        for match, match_count in profile_match_counts.most_common():
            (accession1, accession2) = match.split()
            match_score = max(match_count / profile_counts[accession1], match_count / profile_counts[accession2])
            if match_score > min_match_score:
            # (we compute the maximum score, associated with the rarest profile)
                writer.write(f'{match}\t{match_score:.3}\t{match_count}\t{profile_counts[accession1]}\t{profile_counts[accession2]}\n')
                count += 1
    print(f'Wrote {count}/{len(profile_match_counts)} gene-context associations stronger than {min_match_score:.2}')


def get_match_key(cdd_hit_1: BlastHit, cdd_hit_2: BlastHit) -> str:
    if cdd_hit_1.hit.accession <= cdd_hit_2.hit.accession:
        return f'{cdd_hit_1.hit.accession} {cdd_hit_2.hit.accession}'
    else:
        return f'{cdd_hit_2.hit.accession} {cdd_hit_1.hit.accession}'
