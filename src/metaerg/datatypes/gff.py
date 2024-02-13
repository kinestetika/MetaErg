from urllib.parse import quote, unquote
import re
import gzip

from metaerg.datatypes import fasta
from metaerg.datatypes import sqlite


GFF_ATTRIBUTES = ('id', 'parent', 'subsystems', 'descr', 'taxon', 'signal_peptide', 'tmh_topology')


class GffParser:
    def __init__(self, path, contig_dict:dict=None, target_feature_type_dict:dict=None, inference:str=None):
        self.path = path
        self.contig_dict = contig_dict
        self.target_feature_type_dict = target_feature_type_dict  # should be a dictionary to convert feature types
        self.inference = inference

    def __enter__(self):
        if str(self.path).endswith('.gz'):
            self.handle = gzip.open(self.path, 'rt')
        else:
            self.handle = open(self.path)
        return self

    def __exit__(self, exception_type, exception_value, exception_traceback):
        self.handle.close()

    def __iter__(self):
        while line := self.handle.readline():
            words = line.strip().split('\t')
            match(words):
                case[word, *_] if word.startswith('#'):
                    continue
                case [contig_name, inference, feature_type, start, end, score, strand, frame, qualifiers]:
                    if self.target_feature_type_dict and feature_type not in self.target_feature_type_dict.keys():
                        continue
                    start = int(start) - 1
                    end = int(end)
                    strand = -1 if '-' == strand else 1
                    seq = ''
                    if self.contig_dict:
                        try:
                            contig = self.contig_dict[contig_name]
                            seq = contig['seq'][start:end]
                            if strand < 0:
                                seq = fasta.reverse_complement(seq)
                        except KeyError:
                            continue
                    inference = self.inference if self.inference else inference
                    if self.target_feature_type_dict:
                        feature_type = self.target_feature_type_dict[feature_type]
                    qualifiers = re.split(r"[=;]", qualifiers)
                    qualifiers = [unquote(str) for str in qualifiers]
                    if len(qualifiers) % 2 != 0:
                        qualifiers = qualifiers[:-1]  # this happens for example with prodigal, ending with ";"
                    qualifiers = {qualifiers[i].lower(): qualifiers[i + 1] for i in range(0, len(qualifiers), 2)}
                    feature = sqlite.Feature(contig=contig_name,
                               start=start,
                               end=end,
                               strand=strand,
                               type=feature_type,
                               inference=inference,
                               nt_seq=seq)
                    if 'name' in qualifiers.keys():
                        feature.id = qualifiers['name']
                    if 'id' in qualifiers.keys():
                        feature.id = qualifiers['id']
                    if 'parent' in qualifiers.keys():
                        feature.parent = qualifiers['parent']
                    if 'subsystems' in qualifiers.keys():
                        feature.subsystems = qualifiers['subsystems']
                    if 'descr' in qualifiers.keys():
                        feature.descr = qualifiers['descr']
                    if 'taxon' in qualifiers.keys():
                        feature.taxon = qualifiers['taxon']
                    if 'signal_peptide' in qualifiers.keys():
                        feature.signal_peptide = qualifiers['signal_peptide']
                    if 'tmh_topology' in qualifiers.keys():
                        feature.tmh_topology = qualifiers['tmh_topology']
                    yield feature


def parse_feature_qualifiers_from_gff(qualifier_str) -> {}:
    qal = re.split(r"[=;]", qualifier_str)
    qal = [unquote(str) for str in qal]
    if len(qal) % 2 != 0:
        qal = qal[:-1]  # this happens for example with prodigal which has the qualifier column ending with ";"
    qal = {qal[i].lower(): qal[i + 1] for i in range(0, len(qal), 2)}
    return qal


def gff_write_genome(writer, contig_dict: dict, db_connection):
    for contig_id, contig in contig_dict.items():
        for feature in sqlite.read_all_features(db_connection, contig=contig_id):
            if feature.strand > 0:
                strand = '+'
            elif feature.strand < 0:
                strand = '-'
            else:
                strand = '.'
            attributes = ';'.join([f'{k.capitalize()}={quote(str(v))}' for k, v in feature if v and k in GFF_ATTRIBUTES])
            attributes = 'ID' + attributes[2:]
            writer.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(feature.contig,
                                                                       feature.inference,
                                                                       feature.type,
                                                                       feature.start + 1,
                                                                       feature.end + 1,
                                                                       '.',  # score
                                                                       strand,
                                                                       '.',  # phase
                                                                       attributes
                                                                       ))
