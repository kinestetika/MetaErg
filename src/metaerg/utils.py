import re
import time
import subprocess

from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation

SILENT = False
START_TIME = time.monotonic()
LOG_TOPICS = set()
NON_IUPAC_RE = re.compile(r'[^ACTGN]')

ESCAPE_CHARS = {'%2C': ',',
                '%3B': ';',
                '%3D': '=',
                '%25': '%',
                '%26': '&'}


def format_runtime():
    runtime = time.monotonic() - START_TIME
    return f'[{int(runtime / 3600):02d}h:{int((runtime % 3600) / 60):02d}m:{int(runtime % 60):02d}s]'


def log(log_message, topic=''):
    if not topic or topic in LOG_TOPICS:
        print(f'{format_runtime()} {log_message}')


def filter_seq(record):
    global FILTERED_CONTIGS
    FILTERED_CONTIGS += 1
    record.seq = Seq(NON_IUPAC_RE.sub('N', str(record.seq)))
    record.seq = Seq(''.join(str(record.seq).split()).upper())
    return record


def create_filtered_contig_fasta_file(fasta_file_in:Path, fasta_file_out:Path, min_length=0):
    log('Filtering contigs for length, removing gaps, replacing non-IUPAC bases with N, capitalizing...')
    global FILTERED_CONTIGS
    FILTERED_CONTIGS = 0
    input_seq_iterator = SeqIO.parse(fasta_file_in, "fasta")
    short_seq_iterator = (filter_seq(record) for record in input_seq_iterator if len(record.seq) >= min_length)
    SeqIO.write(short_seq_iterator, fasta_file_out, "fasta")
    log(f'Wrote {FILTERED_CONTIGS} contigs of length >{min_length} '
        f'nt to {fasta_file_out}')


def mask_seq(record:SeqRecord, exceptions=None, min_mask_length=50):
    global MASKED_SEQ, TOTAL_SEQ
    seq = str(record.seq)
    TOTAL_SEQ += len(seq)
    for f in record.features:
        if f.qualifiers['inference'].lower() in exceptions:
            continue
        if len(f.location) < min_mask_length:
            continue
        fl:FeatureLocation = f.location
        MASKED_SEQ += fl.end - fl.start
        seq = seq[:fl.start] + 'N' * (fl.end - fl.start) + seq[fl.end:]
    return SeqRecord(Seq(seq), id=record.id, description=record.description)


def create_masked_contig_fasta_file(fasta_file:Path, contig_dict:dict, exceptions=None, min_mask_length=50):
    if exceptions is None:
        exceptions = set()
    global MASKED_SEQ, TOTAL_SEQ
    (MASKED_SEQ, TOTAL_SEQ) = (0,0)
    log(f'Masking contig features with Ns in file {fasta_file}...')
    seq_iterator = (mask_seq(record, exceptions=exceptions, min_mask_length=min_mask_length)
                    for record in contig_dict.values())
    SeqIO.write(seq_iterator, fasta_file, "fasta")
    log(f'Masked {MASKED_SEQ/TOTAL_SEQ*100:.1f}% of sequence data.')


def unescape_str(str):
    for t in ESCAPE_CHARS.items():
        str = str.replace(t[0], t[1])
    if re.search('%\S\S', str):
        log(f'Warning: potentially missed escaped string in "{str}"')
    return str.strip()


def gff_words_to_seqfeature(words: list):
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
    return SeqFeature(location=loc, type=words[2], qualifiers=qal)


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
        print(result.stderr)
        raise Exception(f'Error while trying to run "{exec}"')


def words_to_hit(words):
    return {'query_id': words[0],
            'hit_id': words[1],
            'percent_id': float(words[2]),
            'aligned_length': int(words[3]),
            'mismatches': int(words[4]),
            'gaps': int(words[5]),
            'query_start': int(words[6]),
            'query_end': int(words[7]),
            'hit_start': int(words[8]),
            'hit_end': int(words[9]),
            'evalue': float(words[10]),
            'score': float(words[11]),
            }


class TabularBlastParser:
    def __init__(self, filename):
        self.filename = filename
        self.next_hit = None
        self.current_query_id = None

    def __enter__(self):
        self.file = open(self.filename)
        if self.load_next_hit_from_file():
            self.current_query_id = self.next_hit['query_id']
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
            self.next_hit = words_to_hit(words)
            return True

    def __next__(self):
        all_hits = list()
        if self.next_hit:
            all_hits.append(self.next_hit)
        while self.load_next_hit_from_file():
            if self.current_query_id != self.next_hit['query_id']:
                self.current_query_id = self.next_hit['query_id']
                return self.current_query_id, all_hits
            all_hits.append(self.next_hit)
        if len(all_hits):
            return self.current_query_id, all_hits
        raise StopIteration