import re
import time
import subprocess
from collections import namedtuple

from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.SeqRecord import SeqRecord

SILENT = False
START_TIME = time.monotonic()
LOG_TOPICS = set()
NON_IUPAC_RE = re.compile(r'[^ACTGN]')
UNWANTED_FEATURE_QUALIFIERS = 'gc_cont conf score cscore sscore rscore uscore tscore rpt_type rpt_family ' \
                              'ltr_similarity seq_number'.split()
VALID_GBK_FEATURE_KEYS = set('CDS tRNA rRNA ncRNA repeat_region retrotransposon crispr_repeat'.split())
NON_CODING_RNA_TYPES = {'LSU_rRNA_bacteria':'rRNA',
                        'LSU_rRNA_archaea': 'rRNA',
                        'LSU_rRNA_eukarya': 'rRNA',
                        'SSU_rRNA_bacteria': 'rRNA',
                        'SSU_rRNA_archaea': 'rRNA',
                        'SSU_rRNA_eukarya': 'rRNA',
                        'SSU_rRNA_microsporidia': 'rRNA',
                        '5S_rRNA': 'rRNA',
                        '5_8S_rRNA': 'rRNA',
                        'tmRNA': 'tmRNA',
                        'tRNA': 'tRNA'}
FEATURE_ID_PATTERN = re.compile("(.+)\.(\d{5})\.(crispr|trna|rna|ltr|tr|repeat|cds)$")
TRANSLATION_TABLE = 11
ESCAPE_CHARS = {'%2C': ',',
                '%3B': ';',
                '%3D': '=',
                '%25': '%',
                '%26': '&'}

def format_runtime():
    runtime = time.monotonic() - START_TIME
    return f'[{int(runtime / 3600):02d}h:{int((runtime % 3600) / 60):02d}m:{int(runtime % 60):02d}s]'


def log(log_message, values=(), topic=''):
    if not topic or topic in LOG_TOPICS:
        if len(values):
            print(f'{format_runtime()} {log_message.format(*values)}')
        else:
            print(f'{format_runtime()} {log_message}')


def get_location_from_gff_words(words):
    strand = -1 if '+' == words[6] else 1
    return FeatureLocation(int(words[3]) - 1, int(words[4]), strand=-1 if '+' == words[6] else 1)


def get_feature_qualifier(feature: SeqFeature, key):
    if key in feature.qualifiers:
        if isinstance(feature.qualifiers[key], list):
            if len(feature.qualifiers[key]) > 0:
                return feature.qualifiers[key][0]
        else:
            return feature.qualifiers[key]
    return ""


def set_feature_qualifier(feature: SeqFeature, key, value):
    feature.qualifiers[key] = [value]


def split_fasta_file(contig_dict, base_file:Path, number_of_files:int, target):
    count = len(contig_dict)
    if 'CDS' == target:
        for contig in contig_dict.values():
            for f in contig.features:
                if f.type == 'CDS':
                    count += 1
    number_of_files = min(number_of_files, count)
    seqs_per_file = count / number_of_files
    paths = [Path(base_file.parent, f'{base_file.name}.{i}') for i in range(number_of_files)]
    filehandles = [open(paths[i], 'w') for i in range(number_of_files)]
    count = 0
    for contig in contig_dict.values():
        if 'CDS' == target:
            for f in contig.features:
                if f.type == 'CDS':
                    feature_seq = pad_seq(f.extract(contig)).translate(table=TRANSLATION_TABLE)[:-1]
                    feature_seq.id = get_feature_qualifier(f, 'id')
                    feature_seq.description = feature_seq.id
                    SeqIO.write(feature_seq, filehandles[int(count / seqs_per_file)], "fasta")
                    count += 1
        else:
            SeqIO.write(contig, filehandles[int(count / seqs_per_file)], "fasta")
            count += 1
    for f in filehandles:
        f.close()
    return paths


def filter_seq(record):
    record.seq = Seq(NON_IUPAC_RE.sub('N', str(record.seq)))
    record.seq = Seq(''.join(str(record.seq).split()).upper())
    return record


def create_filtered_contig_fasta_file(fasta_file_in:Path, fasta_file_out:Path, min_length=0):
    log('Filtering contigs for length, removing gaps, replacing non-IUPAC bases with N, capitalizing...')
    input_seq_iterator = SeqIO.parse(fasta_file_in, "fasta")
    short_seq_iterator = (filter_seq(record) for record in input_seq_iterator if len(record.seq) >= min_length)
    SeqIO.write(short_seq_iterator, fasta_file_out, "fasta")


def unescape_str(str):
    for t in ESCAPE_CHARS.items():
        str = str.replace(t[0], t[1])
    if re.search('%\S\S', str):
        log(f'Warning: potentially missed escaped string in "{str}"')
    return str.strip()


def gff_words_to_seqfeature(words: list, inference=''):
    strand = None
    if '+' == words[6]:
        strand = +1
    elif '-' == words[6]:
        strand = -1

    loc = FeatureLocation(int(words[3]) - 1, int(words[4]), strand=strand)
    qal = re.split(r"[=;]", words[8])
    qal = [unescape_str(str) for str in qal]

    if len(qal) % 2 != 0:
        qal = qal[:-1]  # this happens for example with prodigal which has the qualifier column ending with ";"
    qal = {qal[i].lower(): qal[i + 1] for i in range(0, len(qal), 2)}
    seq_feature = SeqFeature(location=loc, type=words[2], qualifiers=qal)

    if (inference): # inference is only passed by predict functions, not during database constuction
        set_feature_qualifier(seq_feature, "inference", inference)

        if not words[2] in VALID_GBK_FEATURE_KEYS:
            set_feature_qualifier(seq_feature, "note", words[2])
            seq_feature.type = 'misc_feature'

        for unwanted in UNWANTED_FEATURE_QUALIFIERS:
            try:
                del seq_feature.qualifiers[unwanted]
            except KeyError:
                pass

    return seq_feature


def pad_seq(sequence):
    """ Pad sequence to multiple of 3 with N """
    remainder = len(sequence) % 3
    return sequence if remainder == 0 else sequence + Seq('N' * (3 - remainder))


def run_external(exec, stdin=None, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, log_cmd=True):
    if log_cmd:
        log(exec)
    if SILENT:
        return
    result = subprocess.run(exec.split(), stdout=stdout, stdin=stdin, stderr=stderr)
    if result.returncode != 0:
        # print(result.stderr)
        raise Exception(f'Error while trying to run "{exec}"')


def non_coding_rna_hmm_hit_to_seq_record(hit, contig_dict):
    contig = contig_dict[hit['query_id']]
    type = 'ncRNA'
    try:
        type = NON_CODING_RNA_TYPES[hit['hit_id']]
    except KeyError:
        if hit['hit_id'].startswith('CRISPR'):
            type = 'crispr'

    f = SeqFeature(location=FeatureLocation(hit['query_start'] - 1, hit['query_end'], strand=hit['query_strand']),
                   type=type, qualifiers={'name': [hit["hit_id"]],
                                          'inference': ['cmscan'],
                                          'profile': [hit["descr"]]})
    contig.features.append(f)


def decipher_metaerg_id(id):
    m = FEATURE_ID_PATTERN.match(id)
    return {'contig_id': m.group(1),
            'gene_number': int(m.group(2)),
            'gene_type': m.group(3)}


BlastHit = namedtuple('BlastHit', ['query', 'hit', 'percent_id', 'aligned_length', 'mismatches', 'gaps',
                                   'query_start', 'query_end', 'hit_start', 'hit_end', 'evalue', 'score'])

BlastResult = namedtuple('BlastResult', ['query', 'hits'])

class TabularBlastParser:
    def __init__(self, filename):
        self.filename = filename
        self.next_hit: BlastHit
        self.current_query = None

    def __enter__(self):
        self.file = open(self.filename)
        if self.load_next_hit_from_file():
            self.current_query = self.next_hit.query
            return self
        else:
            raise StopIteration

    def __exit__(self, exception_type, exception_value, exception_traceback):
        self.file.close()

    def __iter__(self):
        return self

    def load_next_hit_from_file(self):
        self.next_hit = None
        while self.file:
            line = self.file.readline()
            if not line:
                return False
            if line.startswith("#"):
                continue
            words = line.strip().split('\t')
            if len(words) < 12:
                continue
            self.next_hit = BlastHit(words[0], words[1], float(words[2]), int(words[3]), int(words[4]), int(words[5]),
                    int(words[6]), int(words[7]), int(words[8]), int(words[9]), float(words[10]), float(words[11]))
            return True

    def __next__(self):
        all_hits = list()
        if self.next_hit:
            all_hits.append(self.next_hit)
        while self.load_next_hit_from_file():
            if self.current_query != self.next_hit.query:
                prev_query = self.current_query
                self.current_query = self.next_hit.query
                return BlastResult(prev_query, tuple(all_hits))
            all_hits.append(self.next_hit)
        if len(all_hits):
            return BlastResult(self.current_query, tuple(all_hits))
        raise StopIteration


HmmHit = namedtuple('HmmHit', ['query', 'hit', 'evalue', 'score'])


class TabularHMMSearchParser:
    def __init__(self, filename):
        self.filename = filename
        self.next_hit: BlastHit
        self.current_query = None

    def __enter__(self):
        self.file = open(self.filename)
        if self.load_next_hit_from_file():
            self.current_query = self.next_hit.query
            return self
        else:
            raise StopIteration

    def __exit__(self, exception_type, exception_value, exception_traceback):
        self.file.close()

    def __iter__(self):
        return self

    def load_next_hit_from_file(self):
        self.next_hit = None
        while self.file:
            line = self.file.readline()
            if not line:
                return False
            if line.startswith("#"):
                continue
            words = line.strip().split('\t')
            if len(words) < 18:
                continue
            self.next_hit = HmmHit(words[0], words[2], float(words[4]), float(words[5]))
            return True

    def __next__(self):
        all_hits = list()
        if self.next_hit:
            all_hits.append(self.next_hit)
        while self.load_next_hit_from_file():
            if self.current_query != self.next_hit.query:
                prev_query = self.current_query
                self.current_query = self.next_hit.query
                return BlastResult(prev_query, tuple(all_hits))
            all_hits.append(self.next_hit)
        if len(all_hits):
            return BlastResult(self.current_query, tuple(all_hits))
        raise StopIteration

