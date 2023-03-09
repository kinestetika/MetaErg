from pathlib import Path
from collections import Counter

import openpyxl
from openpyxl.styles import Font

from metaerg import context
from metaerg.datatypes import sqlite
from metaerg.datatypes import functional_genes
from metaerg.calculations.codon_usage_bias import compute_codon_bias_estimate_doubling_time


GENOME_PROPERTY_FORMATS = {'genome name': '<',
                           'input fasta file': '<',
                           '# contigs': ',',
                           'size': ',',
                           '% GC': '.2%',
                           'N50 contig length': ',',
                           '# proteins': ',',
                           '% coding': '.1%',
                           'mean protein length (aa)': '.0f',
                           '# ribosomal RNA': ',',
                           '# transfer RNA': ',',
                           '# non-coding RNA': ',',
                           '# retrotransposons': ',',
                           '# CRISPR repeats': ',',
                           '# other repeats': ',',
                           '# total features': ',',
                           '% repeats': '.1%',
                           'classification (top taxon)': '<',
                           '% of CDS classified to top taxon': '.1%',
                           '% of CDS that could be classified': '.1%',
                           'codon usage bias': '.3f',
                           'doubling_time (days)': '.1f'
                           }


def compute_genome_properties(genome_name: str, input_fasta_file: Path, contig_dict: dict[str, dict], db_connection) -> dict:
    properties = {'genome name':              genome_name,
                  'input fasta file':         input_fasta_file.name,
                  '# total features':         0,
                  '# contigs':                0,
                  'N50 contig length':        0,
                  '# proteins':               0,
                  '# ribosomal RNA':          0,
                  '# transfer RNA':           0,
                  '# non-coding RNA':         0,
                  '# retrotransposons':       0,
                  '# CRISPR repeats':         0,
                  '# other repeats':          0,
                  '% coding':                 0.0,
                  '% repeats':                0.0,
                  'mean protein length (aa)': 0.0,}
    contigs:list[dict] = list(contig_dict.values())
    contigs.sort(key=lambda c:len(c['seq']), reverse=True)
    properties['size'] = sum(len(c['seq']) for c in contigs)
    properties['# contigs'] = len(contigs)
    properties['% GC'] = sum((c['seq'].count('G') + c['seq'].count('G') for c in contigs)) / \
                               (properties['size'] - sum((c['seq'].count('N') for c in contigs)))
    cum_size = 0
    for c in contigs:
        cum_size += len(c['seq'])
        if cum_size >+ properties['size'] / 2:
            properties['N50 contig length'] = len(c['seq'])
            break

    taxon_counts = Counter()
    for feature in sqlite.read_all_features(db_connection):
        properties['# total features'] += 1
        if feature.type == 'CDS':
            properties['# proteins'] += 1
            properties['% coding'] += feature.length_nt()
            properties['mean protein length (aa)'] += feature.length_aa()
            taxon_counts.update((feature.taxon,))
        elif feature.type == 'rRNA':
            properties['# ribosomal RNA'] += 1
            taxon_counts.update((feature.taxon,))
        elif feature.type == 'tRNA':
            properties['# transfer RNA'] += 1
            taxon_counts.update((feature.taxon,))
        elif feature.type == 'ncRNA':
            properties['# non-coding RNA'] += 1
            taxon_counts.update((feature.taxon,))
        elif feature.type == 'retrotransposon':
            properties['# retrotransposons'] += 1
            properties['% repeats'] += feature.length_nt()
        elif feature.type == 'crispr_repeat':
            properties['# CRISPR repeats'] += 1
            properties['% repeats'] += feature.length_nt()
        elif feature.type == 'repeat':
            properties['# other repeats'] += 1
            properties['% repeats'] += feature.length_nt()
    properties['% coding'] /= properties['size']
    properties['mean protein length (aa)'] /= properties['# proteins']
    properties['% repeats'] /= properties['size']
    properties['classification (top taxon)'] = ''
    properties['% of CDS classified to top taxon'] = 0.0
    properties['% of CDS that could be classified'] = 1 - taxon_counts[''] / taxon_counts.total()
    del taxon_counts['']
    for k, v in taxon_counts.most_common(1):
        properties['% of CDS classified to top taxon'] = v / taxon_counts.total()
        properties['classification (top taxon)'] = k
        break
    codon_usage_bias, doubling_time = compute_codon_bias_estimate_doubling_time(db_connection)
    properties['codon usage bias'] = codon_usage_bias
    properties['doubling_time (days)'] = doubling_time
    properties['subsystems'] = functional_genes.aggregate(db_connection)
    properties['subsystem_summary'] = {}
    for subsystem, subsystem_genes in properties['subsystems'].items():
        subsystem_completeness = functional_genes.get_subsystem_completeness(subsystem, subsystem_genes)
        properties['subsystem_summary'][subsystem] = subsystem_completeness
    return properties


