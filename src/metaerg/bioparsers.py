import gzip
import re
import textwrap
from pathlib import Path
from metaerg.data_model import SeqFeature, SeqRecord, Genome, BlastHit, BlastResult, Masker
from metaerg import context
from urllib.parse import unquote

NON_IUPAC_RE_NT = re.compile(r'[^ACTGN]')
NON_IUPAC_RE_AA = re.compile(r'[^RHKDESTNQCUGPAVILMFYW]')


def _parse_fasta_header(line:str) -> SeqRecord:
    si = line.find(' ')
    if si > 0:
        return SeqRecord(id=line[1:si], seq='', descr=line[si+1:].strip())
    else:
        return SeqRecord(id=line[1], seq='')


class FastaParser:
    def __init__(self, path, cleanup_seq=True):
        self.path = path
        self.handle = None
        self.cleanup_seq = cleanup_seq
        self.alphabet = None
        self.unknown_char = ''

    def __enter__(self):
        if str(self.path).endswith('.gz'):
            self.handle = gzip.open(self.path, 'rt')
        else:
            self.handle = open(self.path)
        return self

    def __exit__(self, exception_type, exception_value, exception_traceback):
        self.handle.close()

    def __iter__(self):
        seq_rec = None
        while line := self.handle.readline():
            if line.startswith('>'):
                if seq_rec:
                    seq_rec.seq =  self._cleanup(seq_rec.seq)
                    yield seq_rec
                seq_rec = _parse_fasta_header(line)
            elif seq_rec:
                seq_rec.seq += line.strip()
                if len(seq_rec.seq) > 10 and not self.alphabet and self.cleanup_seq:
                    seq_nt_errors = sum(1 for match in NON_IUPAC_RE_NT.finditer(seq_rec.seq))
                    seq_aa_errors = sum(1 for match in NON_IUPAC_RE_AA.finditer(seq_rec.seq))
                    if seq_nt_errors <= seq_aa_errors:
                        self.alphabet = NON_IUPAC_RE_NT
                        self.unknown_char = 'N'
                    else:
                        self.alphabet = NON_IUPAC_RE_AA
                        self.unknown_char = 'X'
        if seq_rec:
            seq_rec.seq = self._cleanup(seq_rec.seq)
            yield seq_rec

    def _cleanup(self, seq) -> str:
        if self.cleanup_seq:
            seq = seq.upper()
            seq = ''.join(seq.split())
            if self.alphabet:
                seq = self.alphabet.sub(self.unknown_char, seq)
        return seq

    # def __next__(self):
    #     seq_record = None
    #     seq = ''
    #     while line := self.handle.readline():
    #         if line.startswith('>'):
    #             next_header = _parse_fasta_header(line)
    #             if seq:
    #                 seq_record = SeqRecord(id=self.header[0], descr=self.header[1], seq=self._cleanup(seq))
    #             self.header = next_header
    #             if seq_record:
    #                 return seq_record
    #         else:
    #             seq += line.strip()
    #             if len(seq) > 10 and not self.alphabet and self.cleanup_seq:
    #                 seq_nt_errors = sum(1 for match in NON_IUPAC_RE_NT.finditer(seq))
    #                 seq_aa_errors = sum(1 for match in NON_IUPAC_RE_AA.finditer(seq))
    #                 if seq_nt_errors <= seq_aa_errors:
    #                     self.alphabet = NON_IUPAC_RE_NT
    #                     self.unknown_char = 'N'
    #                 else:
    #                     self.alphabet = NON_IUPAC_RE_AA
    #                     self.unknown_char = 'X'
    #     if seq:
    #         return SeqRecord(id=self.header[0], descr=self.header[1], seq=self._cleanup(seq))
    #     raise StopIteration


def init_genome_from_fasta_file(genome_id, filename, min_contig_length=0, delimiter='.') -> Genome:
    with FastaParser(filename) as fasta_reader:
        contigs = [c for c in fasta_reader if len(c) > min_contig_length]
    contigs.sort(key=len, reverse=True)
    contigs = {c.id: c for c in contigs}
    return Genome(genome_id, contigs=contigs, delimiter=delimiter)


def write_fasta(handle, fasta, line_length=80):
    def _wf(f):
        assert len(f.seq), 'Attempt to write zero-lentgh sequence to fasta.'
        handle.write(f'>{f.id} {f.descr}\n')
        wrapper = textwrap.TextWrapper(break_on_hyphens=False, width=line_length)
        seq_lines = wrapper.wrap(text=f.seq)
        for l in seq_lines:
            handle.write(l)
            handle.write('\n')

    try:
        for f in fasta:
            _wf(f)
    except TypeError:
        _wf(fasta)


def write_genome_to_fasta_files(genome, base_file: Path, split=1, targets: tuple = (), mask=True, exceptions=None, min_length=50):
    """writes features (of target FeatureType), or contigs (target = None), to one or more (split) fasta files,
    optionally masking features with N"""
    if targets:
        number_of_records = sum(1 for c in genome.contigs.values() for f in c.features if f.type in targets)
    else:
        number_of_records = len(genome.contigs)
    split = min(split, number_of_records)
    records_per_file = number_of_records / split
    paths = (Path(base_file.parent, f'{base_file.name}.{i}') for i in range(split)) if split > 1 else base_file,
    filehandles = [open(p, 'w') for p in paths]
    number_of_records = 0
    masker = Masker(mask=mask, exceptions=exceptions, min_length=min_length)
    for contig in genome.contigs.values():
        if targets:
            for f in contig.features:
                if f.type in targets:
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
        self.retrieve_db_entry = retrieve_db_entry

    def __enter__(self):
        if str(self.filename).endswith('.gz'):
            self.handle = gzip.open(self.filename, 'rt')
        else:
            self.handle = open(self.filename)
        return self

    def __exit__(self, exception_type, exception_value, exception_traceback):
        self.handle.close()

    def __iter__(self):
        all_hits: list[BlastHit] = []
        while next_hit := self.load_next_hit_from_file():
            if not len(all_hits) or all_hits[-1].query == next_hit.query:
                all_hits.append(next_hit)
            else:
                yield BlastResult(tuple(all_hits))
                all_hits = [next_hit]
        if len(all_hits):
            yield BlastResult(tuple(all_hits))

    def load_next_hit_from_file(self) -> BlastHit | None:
        while line := self.handle.readline():
            words = line.strip().split('\t')
            match words:
                case [str(word), *_] if word.startswith('#'):
                    continue
                case [query, hit, percent_id, aligned_length, mismatches, gaps, query_start, query_end, hit_start,
                      hit_end, evalue, score] if 'BLAST' == self.mode:
                    hit_db_entry = self.retrieve_db_entry(hit)
                    return BlastHit(query, hit_db_entry, float(percent_id), int(aligned_length),
                                    int(mismatches), int(gaps), int(query_start), int(query_end),
                                    int(hit_start), int(hit_end), float(evalue), float(score))
                case [hit, _, query, _, evalue, score, _, _, _, _, _, _, _, _, _, _, _, _, *_] if 'HMMSCAN' == self.mode:
                    hit_db_entry = self.retrieve_db_entry(hit)
                    return BlastHit(query, hit_db_entry, 0, 0, 0, 0, 0, 0, 0, 0, float(evalue), float(score))
                case [query, _, hit, _, evalue, score, _, _, _, _, _, _, _, _, _, _, _, _, *_] if 'HMMSEARCH' == self.mode:
                    hit_db_entry = self.retrieve_db_entry(hit)
                    return BlastHit(query, hit_db_entry, 0, 0, 0, 0, 0, 0, 0, 0, float(evalue), float(score))
                case [*_]:
                    continue
        return None

    # def __next__(self):
    #     all_hits: list[BlastHit] = []
    #     if self.next_hit:
    #         all_hits.append(self.next_hit)
    #     while self.load_next_hit_from_file():
    #         if self.current_query != self.next_hit.query:
    #             self.current_query = self.next_hit.query
    #             return BlastResult(tuple(all_hits))
    #         all_hits.append(self.next_hit)
    #     if len(all_hits):
    #         return BlastResult(tuple(all_hits))
    #     raise StopIteration


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
                        contig: SeqRecord = self.contig_dict[contig_name]
                    except KeyError:
                        continue
                    start = int(start) - 1
                    end = int(end)
                    strand = -1 if '-' == strand else 1
                    seq = contig.seq[start:end]
                    if strand < 0:
                        seq = reverse_complement(seq)
                    inference = self.inference if self.inference else inference
                    qualifiers = re.split(r"[=;]", qualifiers)
                    qualifiers = [unquote(str) for str in qualifiers]
                    if len(qualifiers) % 2 != 0:
                        qualifiers = qualifiers[:-1]  # this happens for example with prodigal, ending with ";"
                    qualifiers = {qualifiers[i].lower(): qualifiers[i + 1] for i in range(0, len(qualifiers), 2)}
                    feature = SeqFeature(start, end, strand, self.target_feature_type_dict[feature_type],
                                         inference=inference, seq=seq)
                    contig.features.append(feature)
                    yield feature


def parse_feature_qualifiers_from_gff(qualifier_str) -> {}:
    qal = re.split(r"[=;]", qualifier_str)
    qal = [unquote(str) for str in qal]
    if len(qal) % 2 != 0:
        qal = qal[:-1]  # this happens for example with prodigal which has the qualifier column ending with ";"
    qal = {qal[i].lower(): qal[i + 1] for i in range(0, len(qal), 2)}
    return qal


COMPLEMENT = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A',}


def reverse_complement(seq:str) -> str:
    return ''.join(COMPLEMENT.get(b, 'N') for b in reversed(seq))


def pad_seq(sequence):
    """ Pad sequence to multiple of 3 with N """
    remainder = len(sequence) % 3
    return sequence if remainder == 0 else sequence + 'N' * (3 - remainder)