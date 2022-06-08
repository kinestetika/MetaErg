import gzip
import re
import collections
from pathlib import Path
from metaerg.data_model import FeatureType, MetaergSeqFeature, MetaergSeqRecord, MetaergGenome, BlastHit, BlastResult, Masker
from metaerg import context
from urllib.parse import unquote

NON_IUPAC_RE = re.compile(r'[^ACTGN]')


def _parse_fasta_header(line:str) -> tuple:
    si = line.find(' ')
    if si > 0:
        return line[1:si], line[si+1:]
    else:
        return line[1], ''


class FastaParser:
    def __init__(self, path, cleanup_seq=True):
        self.path = path
        self.header = None
        self.handle = None
        self.cleanup_seq = cleanup_seq

    def __enter__(self):
        if str(self.path).endswith('.gz'):
            self.handle = gzip.open(self.path, 'rt')
        else:
            self.handle = open(self.path)
        return self

    def __exit__(self, exception_type, exception_value, exception_traceback):
        self.handle.close()

    def __iter__(self):
        return self

    def _cleanup(self, seq):
        if self.cleanup_seq:
            seq = NON_IUPAC_RE.sub('N', seq)
            seq = ''.join(seq.split()).upper()
        return seq

    def __next__(self):
        seq_record = None
        seq = ''
        while line := self.handle.readline():
            if line.startswith('>'):
                next_header = _parse_fasta_header(line)
                if seq:
                    seq_record = MetaergSeqRecord(id=self.header[0], descr=self.header[1], seq=self._cleanup(seq))
                self.header = next_header
                if seq_record:
                    return seq_record
            else:
                seq += line
        if seq:
            return MetaergSeqRecord(id=self.header[0], descr=self.header[1], seq=self._cleanup(seq))
        raise StopIteration


def load_genome_from_file(genome_id, filename, min_contig_length=0, delimiter='.') -> MetaergGenome:
    with FastaParser(filename) as fasta_reader:
        contigs = [c for c in fasta_reader if len(c) > min_contig_length]
    contigs.sort(key=len, reverse=True)
    contigs = {c.id: c for c in contigs}
    return MetaergGenome(genome_id, contigs=contigs, delimiter=delimiter)


def write_fasta(handle, fasta, line_length=80):
    def _wf(f):
        handle.write(f'>{f.id} {f.descr}\n')
        for i in range(0, f.seq, line_length):
            handle.write(f.seq[i:min(len(f.seq), i + line_length)])
            handle.write('\n')

    if isinstance(fasta, collections.Iterable):
        for f in fasta:
            _wf(f)
    else:
        _wf(fasta)


def write_genome_fasta_files(genome, base_file: Path, split=1, target=None, mask=True, exceptions=None, min_length=50):
    """writes features (of target FeatureType), or contigs (target = None), to one or more (split) fasta files,
    optionally masking features with N"""
    if target:
        number_of_records = sum(1 for c in genome.contigs.values() for f in c.features if f.type == target
                                or (isinstance(target, collections.Sequence) and f.type in target))
    else:
        number_of_records = len(genome.contigs)
    split = min(split, number_of_records)
    records_per_file = number_of_records / split
    paths = (Path(base_file.parent, f'{base_file.name}.{i}') for i in range(split)) if split > 1 else base_file,
    filehandles = [open(p, 'w') for p in paths]
    number_of_records = 0
    masker = Masker(mask=mask, exceptions=exceptions, min_length=min_length)
    for contig in genome.contigs.values():
        if target:
            for f in contig.features:
                if f.type == target or (isinstance(target, collections.Sequence) and f.type in target):
                    write_fasta(filehandles[int(number_of_records / records_per_file)], f)
                    number_of_records += 1
        else:
            write_fasta(filehandles[int(number_of_records / records_per_file)], masker.mask(contig))
            number_of_records += 1
    for f in filehandles:
        f.close()
    if mask:
        context.log(masker.stats())
    if split > 1:
        return paths
    else:
        return base_file


class TabularBlastParser:
    def __init__(self, filename, mode, retrieve_db_entry):
        self.filename = filename
        self.mode = mode
        self.handle = None
        self.next_hit = None
        self.current_query = None
        self.retrieve_db_entry = retrieve_db_entry

    def __enter__(self):
        if str(self.filename).endswith('.gz'):
            self.handle = gzip.open(self.filename, 'rt')
        else:
            self.handle = open(self.filename)
        if self.load_next_hit_from_file():
            self.current_query = self.next_hit.query
            return self
        else:
            raise StopIteration

    def __exit__(self, exception_type, exception_value, exception_traceback):
        self.handle.close()

    def __iter__(self):
        return self

    def load_next_hit_from_file(self):
        self.next_hit = None
        while line := self.handle.readline():
            words = line.strip().split('\t')
            match words:
                case [str(word), *_] if word.startswith('#'):
                    continue
                case [query, hit, percent_id, aligned_length, mismatches, gaps, query_start, query_end, hit_start,
                      hit_end, evalue, score] if 'BLAST' == self.mode:
                    hit_db_entry = self.retrieve_db_entry(hit)
                    self.next_hit = BlastHit(query, hit_db_entry, float(percent_id), int(aligned_length),
                                             int(mismatches), int(gaps), int(query_start), int(query_end),
                                             int(hit_start), int(hit_end), float(evalue), float(score))
                case [hit, _, query, _, evalue, score, _, _, _, _, _, _, _, _, _, _, _, _, *_] if 'HMMSCAN' == self.mode:
                    hit_db_entry = self.retrieve_db_entry(hit)
                    self.next_hit = BlastHit(query, hit_db_entry, 0, 0, 0, 0, 0, 0, 0, 0, float(evalue), float(score))
                case [query, _, hit, _, evalue, score, _, _, _, _, _, _, _, _, _, _, _, _, *_] if 'HMMSEARCH' == self.mode:
                    hit_db_entry = self.retrieve_db_entry(hit)
                    self.next_hit = BlastHit(query, hit_db_entry, 0, 0, 0, 0, 0, 0, 0, 0, float(evalue), float(score))
                case [*_]:
                    continue
            return True
        return False

    def __next__(self):
        all_hits: list[BlastHit] = []
        if self.next_hit:
            all_hits.append(self.next_hit)
        while self.load_next_hit_from_file():
            if self.current_query != self.next_hit.query:
                self.current_query = self.next_hit.query
                return BlastResult(tuple(all_hits))
            all_hits.append(self.next_hit)
        if len(all_hits):
            return BlastResult(tuple(all_hits))
        raise StopIteration


class GffParser:
    def __init__(self, path):
        self.path = path

    def __enter__(self):
        if str(self.path).endswith('.gz'):
            self.handle = gzip.open(self.path, 'rt')
        else:
            self.handle = open(self.path)
        return self

    def __exit__(self, exception_type, exception_value, exception_traceback):
        self.handle.close()

    def __iter__(self):
        return self

    def __next__(self):
        while line := self.handle.readline():
            if line.startswith('#'):
                continue
            words: list[str] = line.strip().split('\t')
            if len(words) < 9:
                continue
            qal = parse_feature_qualifiers_from_gff(words[8])
            strand = -1 if '-' == words[6] else 1
            try:
                type = FeatureType[words[2]]
            except KeyError:
                type = FeatureType.ncRNA

            return MetaergSeqFeature(int(words[3])-1, int(words[4]), strand, type, words[1], '', id=words[0])
        raise StopIteration


def parse_feature_qualifiers_from_gff(qualifier_str) -> {}:
    qal = re.split(r"[=;]", qualifier_str)
    qal = [unquote(str) for str in qal]
    if len(qal) % 2 != 0:
        qal = qal[:-1]  # this happens for example with prodigal which has the qualifier column ending with ";"
    qal = {qal[i].lower(): qal[i + 1] for i in range(0, len(qal), 2)}
    return qal
