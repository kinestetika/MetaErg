import re
import gzip
import textwrap

from metaerg.datatypes import sqlite

GBK_LINEWIDTH = 80
GBK_INDENT = 21
GBK_HEADER = '''LOCUS       {:<15} {:>12} bp    DNA              UNK 01-JAN-1980
DEFINITION  {}.
ACCESSION   {}
VERSION     {}
KEYWORDS    .
SOURCE      .
  ORGANISM  .
            .
FEATURES             Location/Qualifiers
'''


def gbk_write_feature(writer, feature: sqlite.Feature):
    indent = ' '*GBK_INDENT
    if feature.strand >= 0:
        location = f'{int(feature.start) + 1}..{int(feature.end)}'
    else:
        location = f'complement({int(feature.start) + 1}..{int(feature.end)})'

    writer.write('     {:<16}{}\n'.format(feature.type, textwrap.fill(location, width = GBK_LINEWIDTH,
                                                                      initial_indent = '',
                                                                      subsequent_indent = indent)))
    gbk_keys = {'locus_tag':  feature.id,
                'inference':  feature.inference,
                'product':    feature.descr,
                'taxonomy':   feature.taxon,
                'subsystems': feature.subsystems,
                'notes':      feature.notes,
                'antismash':  feature.antismash,
                'signal_peptide': feature.signal_peptide,
                'tmh':         feature.tmh,
                'tmh-topology':feature.tmh_topology,
                'translation': feature.aa_seq if feature.type == 'CDS' else ''}
    for k, v in gbk_keys.items():
        if v:
            writer.write(textwrap.fill(f'/{k}="{v}"', width = GBK_LINEWIDTH,
                                       initial_indent = indent, subsequent_indent = indent))
            writer.write('\n')


def gbk_write_genome(writer, contig_dict: dict, db_connection):
    for contig_id, contig in contig_dict.items():
        writer.write(GBK_HEADER.format(contig_id, len(contig['seq']), contig_id, contig_id, contig_id))
        for feature in sqlite.read_all_features(db_connection, contig=contig_id):
            gbk_write_feature(writer, feature)
        writer.write('ORIGIN\n')
        lw = ((GBK_LINEWIDTH - 20) // 10) * 10
        for i in range(0, len(contig['seq']), lw):
            line_seq = contig['seq'][i:i + lw].lower()
            writer.write(f'{i+1:>9}')
            for j in range(lw // 10):
                writer.write(' ')
                writer.write(line_seq[j*10:j*10+10])
            writer.write('\n')
        writer.write('//\n')


class GbkFeatureParser:
    def __init__(self, path):
        self.path = path
        self.handle = None
        self.feature_re = re.compile(r'\s{,6}(\w+)\s+(complement\()*(\d+)\.\.(\d+)\)*')
        self.qualifier_re = re.compile(r'/(\w+?)=(.+)')

    def __enter__(self):
        if str(self.path).endswith('.gz'):
            self.handle = gzip.open(self.path, 'rt')
        else:
            self.handle = open(self.path)
        return self

    def __exit__(self, exception_type, exception_value, exception_traceback):
        self.handle.close()

    def __iter__(self):
        current_contig_name = ''
        current_feature = {}
        current_k = ''
        current_v = ''
        parsing_features = False
        while line := self.handle.readline():
            if line.startswith('ACCESSION'):
                current_contig_name = line.split()[1]
                continue
            elif line.startswith('FEATURES'):
                parsing_features = True
                continue
            elif line.startswith('ORIGIN'):
                if current_feature:
                    if current_v:
                        current_feature[current_k] = current_v.strip('"')
                    yield current_feature
                parsing_features = False
                continue
            elif not parsing_features:
                continue
            line = line.strip('\n')
            if m := self.feature_re.match(line):
                if current_feature:
                    yield current_feature
                current_feature = {'contig': current_contig_name,
                                   'start': int(m.group(3))-1,
                                   'end': int(m.group(4)),
                                   'strand': -1 if m.group(2) else 1,
                                   'type': m.group(1)
                                   }
            elif current_feature:
                line =line.strip()
                if m := self.qualifier_re.fullmatch(line):
                    if current_v:
                        current_feature[current_k] = current_v.strip('"')
                    current_k = m.group(1)
                    current_v = m.group(2)
                else:
                    current_v += f' {line}'


def pad_seq(sequence):
    """ Pad sequence to multiple of 3 with N """
    remainder = len(sequence) % 3
    return sequence if remainder == 0 else sequence + 'N' * (3 - remainder)