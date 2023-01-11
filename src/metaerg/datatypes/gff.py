from urllib.parse import unquote
import re
import gzip
from metaerg.datatypes import fasta
from metaerg.datatypes.sqlite import Feature

class GffParser:
    def __init__(self, path, contig_dict, target_feature_type_dict:dict=None, inference:str=None):
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
                    try:
                        contig = self.contig_dict[contig_name]
                    except KeyError:
                        continue
                    start = int(start) - 1
                    end = int(end)
                    strand = -1 if '-' == strand else 1
                    seq = contig['seq'][start:end]
                    if strand < 0:
                        seq = fasta.reverse_complement(seq)
                    inference = self.inference if self.inference else inference
                    qualifiers = re.split(r"[=;]", qualifiers)
                    qualifiers = [unquote(str) for str in qualifiers]
                    if len(qualifiers) % 2 != 0:
                        qualifiers = qualifiers[:-1]  # this happens for example with prodigal, ending with ";"
                    qualifiers = {qualifiers[i].lower(): qualifiers[i + 1] for i in range(0, len(qualifiers), 2)}
                    feature = Feature(contig=contig_name,
                               start=start,
                               end=end,
                               strand=strand,
                               type=self.target_feature_type_dict[feature_type],
                               inference=inference,
                               nt_seq=seq)
                    yield feature


def parse_feature_qualifiers_from_gff(qualifier_str) -> {}:
    qal = re.split(r"[=;]", qualifier_str)
    qal = [unquote(str) for str in qal]
    if len(qal) % 2 != 0:
        qal = qal[:-1]  # this happens for example with prodigal which has the qualifier column ending with ";"
    qal = {qal[i].lower(): qal[i + 1] for i in range(0, len(qal), 2)}
    return qal
