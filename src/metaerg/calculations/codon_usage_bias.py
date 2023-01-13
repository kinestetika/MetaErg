import math

import re
from time import sleep
from statistics import median, mean
from math import log, log10
from pathlib import Path
from collections import namedtuple

from tqdm import tqdm

from metaerg import context
from metaerg.datatypes.fasta import FastaParser
from metaerg.datatypes import ncbi_ftp
from metaerg.datatypes import sqlite

CUBData = namedtuple('CUBData', ['genome_id', 'nHE', 'CUBHE', 'consistency', 'CPB', 'filtered', 'd'])

def compute_feature_codon_counts_by_aa(feature, codon_counts_by_aa: dict|int = 0):
    if not isinstance(codon_counts_by_aa, dict):
        codon_counts_by_aa = {}

    # to deal with start codon, always translated as M:
    aa_seq = feature.aa_seq[0:1].lower() + feature.aa_seq[1:].upper()
    nt_seq = feature.nt_seq[0:3].lower() + feature.nt_seq[3:].upper()

    for i in range(1, len(aa_seq)):
        aa = aa_seq[i]
        codon = nt_seq[i * 3:i * 3 + 3]
        if 'N' in codon:
            continue
        if 'X' == aa:
            continue
        try:
            codons_for_aa = codon_counts_by_aa[aa]
        except KeyError:
            codons_for_aa = {}
            codon_counts_by_aa[aa] = codons_for_aa
        try:
            prev_count = codons_for_aa[codon]
        except:
            prev_count = 0
            codons_for_aa[codon] = 0
        codons_for_aa[codon] = prev_count + 1

    return codon_counts_by_aa


def convert_counts_to_frequencies(codon_counts_by_aa):
    codon_freq_by_aa = {}
    for aa, codons_for_aa in codon_counts_by_aa.items():
        codon_freq_by_aa[aa] = {}
        total_counts_for_aa = sum(codons_for_aa.values())
        for codon, count in codons_for_aa.items():
            codon_freq_by_aa[aa][codon] = count / total_counts_for_aa
    return codon_freq_by_aa


def compute_codon_frequencies_for_feature_data(db_connection, additional_sql: str = '') -> dict:
    genome_wide_codon_counts_by_aa = {}
    for feature in sqlite.read_all_features(db_connection, type='CDS', additional_sql=additional_sql):
        compute_feature_codon_counts_by_aa(feature, genome_wide_codon_counts_by_aa)
    genome_wide_codon_freq_by_aa = convert_counts_to_frequencies(genome_wide_codon_counts_by_aa)
    return genome_wide_codon_freq_by_aa


def compute_codon_usage_bias_for_feature(feature, genome_wide_codon_frequencies_by_aa) -> float:
    # see https://link.springer.com/article/10.1186/1471-2105-6-182, https://github.com/BioinfoHR/coRdon/blob/master/R/codonUsage.R, lines 218-255
    cds_codon_counts_by_aa = compute_feature_codon_counts_by_aa(feature)
    cds_codon_freq_by_aa = convert_counts_to_frequencies(cds_codon_counts_by_aa)
    milc = 0
    correction_factor = 0
    for aa, codons_for_aa in cds_codon_freq_by_aa.items():
        correction_factor += len(genome_wide_codon_frequencies_by_aa[aa]) - 1
        for codon, freq in codons_for_aa.items():
            milc += 2 * cds_codon_counts_by_aa[aa][codon] * log(freq / genome_wide_codon_frequencies_by_aa [aa][codon])
    milc /= len(feature.aa_seq) - 1
    correction_factor = correction_factor / (len(feature.aa_seq) - 1) - 0.5  # -0.5 as in coRdon, not +0.5 as in paper
    return milc - correction_factor


def compute_codon_pair_scores_for_feature_data(genome_id: str, db_connection, additional_sql: str = '') -> dict:
    # see https://www.science.org/doi/full/10.1126/science.1155761, supplement, Fig. 1
    codon_pair_counts = {}
    aa_pair_counts = {}
    aa_counts = {}
    codon_counts = {}
    codon2aa = {}

    ambiguous_codon_warnings = 0
    for feature in sqlite.read_all_features(db_connection, type='CDS', additional_sql=additional_sql):
        # to deal with start codon, always translated as M:
        aa_seq = feature.aa_seq[0:1].lower() + feature.aa_seq[1:].upper()
        nt_seq = feature.nt_seq[0:3].lower() + feature.nt_seq[3:].upper()

        for i in range(1, len(aa_seq)):
            x = aa_seq[i]
            a = nt_seq[i * 3:i * 3 + 3]
            if 'N' in a:
                continue
            if 'X' == x:
                continue
            try:
                prev_x = codon2aa[a]
                if x != prev_x:  # this very rarely happens
                    ambiguous_codon_warnings += 1
                    break
            except KeyError:
                pass
            codon2aa[a] = x
            try:
                aa_counts[x] += 1
            except KeyError:
                aa_counts[x] = 1
            try:
                codon_counts[a] += 1
            except KeyError:
                codon_counts[a] = 1
            try:
                y = aa_seq[i+1]
                b = nt_seq[(i + 1) * 3:(i + 1) * 3 + 3]
                if 'N' in b:
                    continue
                if 'X' == y:
                    continue
            except IndexError:
                break
            try:
                codon_pair_counts[a + b] += 1
            except KeyError:
                codon_pair_counts[a + b] = 1
            try:
                aa_pair_counts[x + y] += 1
            except KeyError:
                aa_pair_counts[x + y] = 1
    if ambiguous_codon_warnings > 0:
        print(
            f'Warning, while computing codon pair scores, skipped {ambiguous_codon_warnings}/{len(feature_data)} '
            f'orf with identical codon translated to different AAs in {genome_id}')

    codon_pair_scores = {}
    for ab, ab_counts in codon_pair_counts.items():
        a = ab[:3]
        b = ab[3:]
        x = codon2aa[a]
        y = codon2aa[b]
        codon_pair_scores[ab] = log(ab_counts / (codon_counts[a] * codon_counts[b] / (aa_counts[x] * aa_counts[y]) * aa_pair_counts[x+y]))
    return codon_pair_scores


def compute_codon_pair_bias_for_cds(feature, codon_pair_scores) -> float:
    # see https://www.science.org/doi/full/10.1126/science.1155761, supplement, Fig. 1
    cpb = 0
    nt_seq = feature.nt_seq[0:3].lower() + feature.nt_seq[3:].upper()
    for i in range(1, len(feature.aa_seq)):
        codon_pair = nt_seq[i * 3: i * 3 + 6]
        if 'N' in codon_pair:
            continue
        try:
            cpb += codon_pair_scores[codon_pair]
        except KeyError:
            pass  # can happen, a codon pair that was not seen before
    return cpb / len(feature.aa_seq)


def compute_codon_usage_bias_for_genome(genome_id, db_connection) -> CUBData:
    # codon usage bias
    background_frequencies = compute_codon_frequencies_for_feature_data(db_connection, additional_sql='end - start >= 240')
    cub_ribosomal = [compute_codon_usage_bias_for_feature(rp, background_frequencies)
                               for rp in sqlite.read_all_features(db_connection, type='CDS',
                               additional_sql='end - start >= 240 AND subsystems LIKE "%ribosomal protein%"')]
    if len(cub_ribosomal):
        codon_usage_bias = median(cub_ribosomal)
        # consistency in codon usage bias
        background_frequencies = compute_codon_frequencies_for_feature_data(db_connection,
                                 additional_sql='end - start >= 240 AND subsystems LIKE "%ribosomal protein%"')
        consistency = mean([compute_codon_usage_bias_for_feature(rp, background_frequencies)
                               for rp in sqlite.read_all_features(db_connection, type='CDS',
                               additional_sql='end - start >= 240 AND subsystems LIKE "%ribosomal protein%"')])
        # codon pair bias
        codon_pair_scores = compute_codon_pair_scores_for_feature_data(genome_id, db_connection, additional_sql='end - start >= 240')
        # codon_pair_bias = sum(codon_pair_scores.values()) / (len(codon_pair_scores) - 1)
        codon_pair_bias = median([compute_codon_pair_bias_for_cds(rp, codon_pair_scores)
                                  for rp in sqlite.read_all_features(type='CDS', additional_sql='end - start >= 240 AND subsystems LIKE "%ribosomal protein%"')])
        # keep track of # too short proteins filtered out and ribosomal proteins
        filtered_out = sum(1 for feature in sqlite.read_all_features(db_connection, type='CDS', additional_sql='end - start < 240'))
        return CUBData(genome_id, len(cub_ribosomal), codon_usage_bias, consistency, codon_pair_bias, filtered_out, 0)
    else:
        return CUBData(genome_id, 0, 0, 0, 0, 0, 0)


def compute_codon_bias_estimate_doubling_time(db_connection):
    background_frequencies = compute_codon_frequencies_for_feature_data(db_connection, additional_sql='end - start >= 240')
    cub_ribosomal = [compute_codon_usage_bias_for_feature(rp, background_frequencies)
                               for rp in sqlite.read_all_features(db_connection, type='CDS',
                               additional_sql='end - start >= 240 AND subsystems LIKE "%ribosomal protein%"')]
    if len(cub_ribosomal):
        codon_usage_bias = median(cub_ribosomal)
        doubling_time = math.pow(10, codon_usage_bias * -2.41937374 + 2.00361692)
        return codon_usage_bias, doubling_time
    else:
        context.log('Warning: No ribosomomal proteins for genome.')
        return 0, 0


def parse_training_data(list_file) -> list[CUBData]:
    training_data = []
    refseq_acc_pattern = re.compile(r'(GCF_\d+)')
    with open(list_file) as reader:
        for line in reader:
            words = line.split('\t')
            genome_id = words[-1]
            if m := refseq_acc_pattern.match(genome_id):
                genome_id = m.group(1)
                training_data.append(CUBData(genome_id, int(words[1]), float(words[2]), float(words[3]), float(words[4]),
                                        int(words[5]), float(words[6])))
    return training_data


def download_training_data(training_data: list[CUBData], destination_dir: Path):
    for entry in tqdm(training_data):
        accession = entry.genome_id
        targets = []
        extension = '_protein.faa.gz'
        aa_file = destination_dir / (accession + extension)
        if not aa_file.exists():
            targets.append((extension, aa_file))
        extension = '_cds_from_genomic.fna.gz'
        fna_file = destination_dir / (accession + extension)
        if not fna_file.exists():
            targets.append((extension, fna_file))
        if targets:
            ncbi_ftp.fetch(accession, targets)
            sleep(2)


def load_feature_data(dir: Path, entry: CUBData):
    nt_file = dir / (entry.genome_id + '_cds_from_genomic.fna.gz')
    aa_file = dir / (entry.genome_id + '_protein.faa.gz')
    aa_seq_dict = {}
    feature_data = []
    try:
        with FastaParser(aa_file, cleanup_seq=False) as reader:
            for seq in reader:
                aa_seq_dict[seq['id']] = seq

        pattern_prot_id = re.compile(r'\[protein_id=(.+?)]')
        with FastaParser(nt_file, cleanup_seq=False) as reader:
            for seq in reader:
                if m:= pattern_prot_id.search(seq['descr']):
                    prot_id = m.group(1)
                    feature_data.append({'id': prot_id,
                                         'subsystems': aa_seq_dict[prot_id]['descr'],
                                         'nt_seq': seq['seq'],
                                         'aa_seq': aa_seq_dict[prot_id]['seq'],
                                         'type': 'CDS'})
        if not feature_data:
            raise Exception(f'No data for {entry.genome_id}.')
        return pd.DataFrame(feature_data)
    except Exception:
        aa_file.unlink(missing_ok=True)
        nt_file.unlink(missing_ok=True)
        return None


def main():
    dir = Path('/bio/databases/eggo-data/')
    training_data = parse_training_data(dir / 'CodonStatistics_Training.csv')
    # download_training_data(training_data, Path('/bio/databases/eggo-data/genomes'))
    count = 0
    with open(dir / 'training_plus_metaerg.csv', 'w') as output:
        for entry in training_data:
            df = load_feature_data(dir / 'genomes', entry)
            if not df is None:
                metaerg_entry = compute_codon_usage_bias_for_genome(entry.genome_id, df)
                print(f'{count}/{len(training_data)} {metaerg_entry.nHE}/{entry.nHE}, '
                      f'{metaerg_entry.CUBHE:.3f}/{entry.CUBHE:.3f} '
                      f'{metaerg_entry.consistency:.3f}/{entry.consistency:.3f} '
                      f'{metaerg_entry.CPB:.3f}/{entry.CPB:.3f} '
                      f'{entry.genome_id}')
                count += 1
                output.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(entry.genome_id,
                                                                          entry.nHE,
                                                                          entry.filtered,
                                                                          entry.CUBHE,
                                                                          entry.consistency,
                                                                          entry.CPB,
                                                                          entry.d,
                                                                          metaerg_entry.nHE,
                                                                          metaerg_entry.CUBHE,
                                                                          metaerg_entry.consistency,
                                                                          metaerg_entry.CPB,
                                                                          ))


if __name__ == "__main__":
    main()