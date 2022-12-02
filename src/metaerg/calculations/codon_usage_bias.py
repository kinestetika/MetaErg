import pandas as pd
import re
from time import sleep
from statistics import median, mean
from math import log
from pathlib import Path
from collections import namedtuple

from tqdm import tqdm

from metaerg.datatypes.fasta import FastaParser
from metaerg.datatypes import ncbi_ftp

CUBData = namedtuple('CUBData', ['genome_id', 'nHE', 'CUBHE', 'consistency', 'CPB', 'filtered', 'd'])

def compute_feature_codon_counts_by_aa(feature, codon_counts_by_aa = 0):
    if not isinstance(codon_counts_by_aa, dict):
        codon_counts_by_aa = {}
    for i in range(1, len(feature.aa_seq)):
        aa = feature.aa_seq[i].upper()
        codon = feature.nt_seq[i * 3:i * 3 + 3].upper()
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


def compute_codon_frequencies_for_feature_data(feature_data: pd.DataFrame) -> dict:
    genome_wide_codon_counts_by_aa = {}
    for feature in feature_data.itertuples():
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
            milc += cds_codon_counts_by_aa[aa][codon] * log(freq / genome_wide_codon_frequencies_by_aa [aa][codon])
    milc /= (len(feature.aa_seq) - 1)
    correction_factor = correction_factor / (len(feature.aa_seq) - 1) - 0.5  # -0.5 as in coRdon, not +0.5 as in paper
    return milc - correction_factor


def compute_codon_pair_scores_for_feature_data(feature_data: pd.DataFrame) -> dict:
    # see https://www.science.org/doi/full/10.1126/science.1155761, supplement, Fig. 1
    codon_pair_counts = {}
    aa_pair_counts = {}
    aa_counts = {}
    codon_counts = {}
    codon2aa = {}
    for feature in feature_data.itertuples():
        for i in range(1, len(feature.aa_seq)):
            x = feature.aa_seq[i].upper()
            a = feature.nt_seq[i * 3:i * 3 + 3].upper()
            if 'N' in a:
                continue
            if 'X' == x:
                continue
            try:
                prev_x = codon2aa[a]
                if x != prev_x:
                    raise Exception(f'Codon {a} conflicted between {prev_x} and {x}')
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
                y = feature.aa_seq[i+1].upper()
                b = feature.nt_seq[(i + 1) * 3:(i + 1) * 3 + 3].upper()
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
    codon_pair_scores = {}
    for ab, ab_counts in codon_pair_counts.items():
        a = ab[:3]
        b = ab[3:]
        x = codon2aa[a]
        y = codon2aa[b]
        codon_pair_scores[ab] = log(ab_counts / (codon_counts[a] * codon_counts[b] / aa_counts[x] / aa_counts[y] * aa_pair_counts[x+y]))
    return codon_pair_scores


def compute_codon_pair_bias_for_cds(feature, codon_pair_scores) -> float:
    # see https://www.science.org/doi/full/10.1126/science.1155761, supplement, Fig. 1
    cpb = 0
    for i in range(1, len(feature.aa_seq)-1):
        codon_pair = feature.nt_seq[i * 3: i * 3 + 6]
        if 'N' in codon_pair:
            continue
        cpb += codon_pair_scores[codon_pair]
    return cpb / len(feature.aa_seq)-1


def compute_codon_usage_bias_for_genome(genome_id, feature_data: pd.DataFrame) -> CUBData:
    feature_data = feature_data[feature_data['type'] == 'CDS']
    filtered_feature_data = feature_data[feature_data['aa_seq'].str.len() >= 80]
    ribosomal_proteins = feature_data[feature_data['subsystems'].str.contains('ribosomal protein')]
    # codon usage bias
    background_frequencies = compute_codon_frequencies_for_feature_data(filtered_feature_data)
    codon_usage_bias = median([compute_codon_usage_bias_for_feature(rp, background_frequencies)
                               for rp in ribosomal_proteins.itertuples()])
    # consistency in codon usage bias
    background_frequencies = compute_codon_frequencies_for_feature_data(ribosomal_proteins)
    consistency = mean([compute_codon_usage_bias_for_feature(rp, background_frequencies)
                        for rp in ribosomal_proteins.itertuples()])
    # codon pair bias
    codon_pair_scores = compute_codon_pair_scores_for_feature_data(filtered_feature_data)
    codon_pair_bias = median([compute_codon_pair_bias_for_cds(rp, codon_pair_scores)
                              for rp in ribosomal_proteins.itertuples()])
    # keep track of # too short proteins filtered out and ribosomal proteins
    filtered_out = len(feature_data[feature_data['aa_seq'].str.len() < 80])

    return CUBData(genome_id, len(ribosomal_proteins.index), codon_usage_bias,
                   consistency, codon_pair_bias, filtered_out, 0)


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
    for entry in training_data:
        df = load_feature_data(dir / 'genomes', entry)
        if not df is None:
            metaerg_entry = compute_codon_usage_bias_for_genome(entry.genome_id, df)
            print(f' {metaerg_entry.nHE}/{entry.nHE}, {metaerg_entry.CUBHE:.3f}/{entry.CUBHE:.3f}')
            count += 1
    print(count)


if __name__ == "__main__":
    main()