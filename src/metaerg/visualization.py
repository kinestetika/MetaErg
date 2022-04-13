import re
import pandas as pd
from Bio import SeqIO

DATAFRAME_COLUMNS = ['id', 'contig', 'start', 'end', 'strand', 'length', 'type', 'taxonomy', 'taxon',
                     'sequence', 'inference', 'note', 'cdd_hits', 'profile',
                     'transmembrane_helixes', 'tmh_topology', 'signal peptide',
                     'antismash_region', 'antismash_function', 'antismash_category',
                     'product', 'blast_aligned', 'blast %id', 'blast hit count', 'description']
PRODUCT_RE = re.compile('\[(\d+/\d+)\w\w@([\d,.]+)%\] (\[\d+/\d+\]) (.+)')

def get_qualifier(feature, key):
    try:
        if isinstance(feature.qualifiers[key], list):
            return feature.qualifiers[key][0]
        else:
            return feature.qualifiers[key]
    except KeyError:
        return ""

def read_genome_from_gbk(gbk_file, contig_dict):
    features_as_dict = []
    with open(gbk_file) as handle:
        for gb_record in SeqIO.parse(handle, "genbank"):
            contig_dict[gb_record.id] = gb_record
            for feature in gb_record.features:
                f = {}
                for key in DATAFRAME_COLUMNS:
                    f[key] = get_qualifier(feature, key)
                # coordinates
                f['start'] = feature.location.start
                f['end'] =  feature.location.end
                f['strand'] = feature.location.strand
                f['contig'] = gb_record.id
                # type, sequence and length
                f['type'] = feature.type
                if f['type'] == 'CDS':
                    f['sequence'] = get_qualifier(feature, 'translation')
                else:
                    f['sequence'] = str(feature.extract(gb_record).seq)
                f['length'] = len(f['sequence'])
                # taxonomy
                f['taxon'] = ''
                if "~" in f['taxonomy']:
                    for t in reversed(f['taxonomy'].split("~")):
                        if not " " in t:
                            f['taxon'] = t
                            break
                # blast hit, product
                match = re.match(PRODUCT_RE, f['product'])
                if (match):
                    f['blast_aligned'] =match.group(1)
                    f['blast %id'] = float(match.group(2))
                    f['blast hit count'] = match.group(3)
                    f['description'] = match.group(4)
                if f['profile'] and not f['description']:
                   f['description'] = f['profile']

                features_as_dict.append(f)
                if len(features_as_dict) >= 500:
                    break
    features = pd.DataFrame(features_as_dict, columns=DATAFRAME_COLUMNS)
    features.set_index('id')
    print(f'loaded {len(features)} features from file.')
    return features
