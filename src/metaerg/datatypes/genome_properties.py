import re
from pathlib import Path
from collections import Counter

import openpyxl

from metaerg import context
from metaerg.datatypes import sqlite
from metaerg.datatypes import functional_genes
from metaerg.datatypes.fasta import load_contigs
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
    return properties


def write_genome_properties_to_xls(genome_property_dict: dict):
    excel_file = context.BASE_DIR / 'genome_properties.xls'
    if not genome_property_dict:
        db_files = [context.spawn_file('annotations.sqlite', genome_name, context.BASE_DIR)
                         for genome_name in context.GENOME_NAMES]
        fna_files = [context.spawn_file("fna", genome_name, context.BASE_DIR)
                     for genome_name in context.GENOME_NAMES]
        genome_name_mappings = {}
        with open(context.GENOME_NAME_MAPPING_FILE) as name_mappings:
            for line in name_mappings:
                w = line.split()
                genome_name_mappings[w[0]] = Path(w[2])
        genome_property_dict = {}
        for db_file, contig_file in zip(db_files, fna_files):
            genome_name = contig_file.stem
            contig_dict = load_contigs(genome_name, contig_file, delimiter='xxxx')
            db_connection = sqlite.connect_to_db(db_file)
            genome_property_dict[genome_name] = compute_genome_properties(genome_name,
                                                                          genome_name_mappings.get(genome_name, 'unknown'),
                                                                          contig_dict, db_connection)

    e_workbook = openpyxl.Workbook()
    e_column = 3
    rows_for_taxon = max(len(gp['classification (top taxon)'].split('; ')) for gp in genome_property_dict.values())
    for genome_name, genome_properties in genome_property_dict.items():
        e_sheet = e_workbook.active
        row = 1
        for k in GENOME_PROPERTY_FORMATS.keys():
            if 'subsystems' == k or 'classification (top taxon)' == k:
                continue
            e_sheet.cell(row=row, column=1).value = k
            value = genome_properties[k]
            if '%' in GENOME_PROPERTY_FORMATS[k]:
                value *= 100
            if m := re.search(r'\d', GENOME_PROPERTY_FORMATS[k]):
                value = round(value, int(m.group(0)))
            e_sheet.cell(row=row, column=e_column).value = value
            row += 1

        e_sheet.cell(row=row, column=1).value = 'classification (top taxon)'
        for k, w in zip('kingdom phylum class order family genus species rank rank rank rank'.split(),
                        genome_properties['classification (top taxon)'].split('; ')):
            e_sheet.cell(row=row, column=2).value = k
            e_sheet.cell(row=row, column=e_column).value = w
            row += 1
        row += max(0, rows_for_taxon - len(genome_properties['classification (top taxon)'].split('; ')))

        for subsystem, genes in genome_properties['subsystems'].items():
            e_sheet.cell(row=row, column=1).value = subsystem
            for gene, feature_ids in genes.items():
                e_sheet.cell(row=row, column=2).value = gene
                e_sheet.cell(row=row, column=e_column).value = ', '.join(feature_ids)
                row += 1
        e_column += 1
    e_workbook.save(excel_file)


