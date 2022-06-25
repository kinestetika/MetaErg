import gzip
import re
import textwrap
import pandas as pd
from pathlib import Path
from metaerg.data_model import SeqFeature, SeqRecord, Genome, BlastHit, BlastResult, FeatureType
from metaerg import context
from urllib.parse import unquote

NON_IUPAC_RE_NT = re.compile(r'[^ACTGN]')
NON_IUPAC_RE_AA = re.compile(r'[^RHKDESTNQCUGPAVILMFYW]')
ALL_MASK_TARGETS = 'CDS rRNA tRNA tmRNA ncRNA repeat crispr_repeat retrotransposon'.split()

DATAFRAME_COLUMNS = 'genome contig id start end strand type inference subsystem descr taxon notes seq antismash ' \
                    'signal_peptide tmh tmh_topology blast cdd antismash'.split()
DATAFRAME_DATATYPES = str, str, str, int, int, int, "category", str, str, str, str, str, str, str, \
                      str, int, str, str, str, str
DATAFRAME_DATATYPES = {k: v for k, v in zip(DATAFRAME_COLUMNS, DATAFRAME_DATATYPES)}

def _parse_fasta_header(line:str) -> SeqRecord:
    si = line.find(' ')
    if si > 0:
        return SeqRecord(id=line[1:si], seq='', descr=line[si+1:].strip())
    else:
        return SeqRecord(id=line[1:].strip(), seq='')


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
        seq = []
        while line := self.handle.readline():
            line = line.strip()
            if line.startswith('>'):
                if len(seq) and seq_rec is  not None:
                    seq_rec.seq = self._cleanup(''.join(seq))
                    seq = []
                    yield seq_rec
                seq_rec = _parse_fasta_header(line)
            elif seq_rec is not None:
                seq.append(line)
                if len(line) > 10 and not self.alphabet and self.cleanup_seq:
                    seq_nt_errors = sum(1 for match in NON_IUPAC_RE_NT.finditer(line))
                    seq_aa_errors = sum(1 for match in NON_IUPAC_RE_AA.finditer(line))
                    if seq_nt_errors <= seq_aa_errors:
                        self.alphabet = NON_IUPAC_RE_NT
                        self.unknown_char = 'N'
                    else:
                        self.alphabet = NON_IUPAC_RE_AA
                        self.unknown_char = 'X'
        if seq_rec is not None:
            seq_rec.seq = self._cleanup(''.join(seq))
            yield seq_rec

    def _cleanup(self, seq) -> str:
        if self.cleanup_seq:
            seq = seq.upper()
            seq = ''.join(seq.split())
            if self.alphabet:
                seq = self.alphabet.sub(self.unknown_char, seq)
        return seq


class Masker:
    def __init__(self, feature_data, targets=None, min_length=50):
        if targets:
            self.feature_data = feature_data.loc[lambda df: df['type'] in targets, :]
            self.apply_mask = True
        else:
            self.feature_data = feature_data
            self.apply_mask = False
        self.min_length = min_length
        self.nt_total = 0
        self.nt_masked = 0

    def mask(self, seq_record: dict) -> dict:
        seq = seq_record['seq']
        seq_record.nt_masked = 0
        feature_data = self.feature_data.loc[lambda df: df['contig'] == seq_record['id'], :]
        if self.apply_mask:
            for index, feature in feature_data.iterrows():
                feature_length = feature['start']-feature['end']
                seq = seq[:feature['start']] + 'N' * feature_length + seq[feature['end']:]
                self.nt_masked += feature_length
        self.nt_total += len(seq)
        return {'id':    seq_record['id'],
                'descr': seq_record['descr'],
                'seq':   seq}

    def stats(self):
        return f'Masked {self.nt_masked / max(self.nt_total, 1) * 100:.1f}% of sequence data.'


def init_genome_from_fasta_file(genome_id, filename, min_contig_length=0, delimiter='.') -> Genome:
    with FastaParser(filename) as fasta_reader:
        contigs = [c for c in fasta_reader if len(c) > min_contig_length]
    contigs.sort(key=len, reverse=True)
    contigs = {c.id: c for c in contigs}
    return Genome(genome_id, contigs=contigs, delimiter=delimiter)


def write_fasta(handle, fasta, line_length=80):
    def _wf(f):
        if not f:
            raise Exception('Attempt to write zero-length sequence to fasta.')
        handle.write(f'>{f.id} {f.descr}\n')
        for i in range(0, len(f.seq), line_length):
            handle.write(f.seq[i:i+line_length])
            handle.write('\n')
    if isinstance(fasta, SeqRecord) or isinstance(fasta, SeqFeature):
        _wf(fasta)
    else:
        for f in fasta:
            _wf(f)


def write_features_to_fasta(feature_data: pd.DataFrame, base_file: Path, split=1, targets: tuple = ()):
    if targets:
        feature_data = feature_data.loc[lambda df: df['type'] in targets, :]
    number_of_records = len(feature_data.index)
    split = min(split, number_of_records)
    records_per_file = number_of_records / split
    if split > 1:
        paths = [Path(base_file.parent, f'{base_file.name}.{i}') for i in range(split)]
    else:
        paths = base_file,
    filehandles = [open(p, 'w') for p in paths]
    records_written = 0
    for index, feature in feature_data.iterrows():
        write_fasta(filehandles[int(records_written / records_per_file)], feature)
        records_written += 1


def write_contigs_to_fasta(genome_name, contig_hash: dict, feature_data: pd.DataFrame, base_file: Path, split=1,
                           mask_targets=None, mask_min_length=50):
    """writes contigs to fasta file(s), optionally masking features with N"""
    number_of_records = len(contig_hash)
    split = min(split, number_of_records)
    records_per_file = number_of_records / split
    if split > 1:
        paths = [Path(base_file.parent, f'{base_file.name}.{i}') for i in range(split)]
    else:
        paths = base_file,
    filehandles = [open(p, 'w') for p in paths]
    records_written = 0
    masker = Masker(feature_data, targets=mask_targets, min_length=mask_min_length)
    for contig in contig_hash.values():
        write_fasta(filehandles[int(records_written / records_per_file)], masker.mask(contig))
        records_written += 1
    for f in filehandles:
        f.close()
    if mask_targets:
        context.log(f'({genome_name}) {masker.stats()}')
    return paths


def write_genome_to_fasta_files(genome, base_file: Path, split=1, targets: tuple = (), mask=False, mask_exceptions=None,
                                mask_min_length=50):
    """writes features (of target FeatureType), or contigs (target = None), to one or more (split) fasta files,
    optionally masking features with N"""
    if targets:
        number_of_records = sum(1 for c in genome.contigs.values() for f in c.features if f.type in targets)
    else:
        number_of_records = len(genome.contigs)
    split = min(split, number_of_records)
    records_per_file = number_of_records / split
    if split > 1:
        paths = [Path(base_file.parent, f'{base_file.name}.{i}') for i in range(split)]
    else:
        paths = base_file,
    filehandles = [open(p, 'w') for p in paths]
    records_written = 0
    masker = Masker(mask=mask, exceptions=mask_exceptions, min_length=mask_min_length)
    for contig in genome.contigs.values():
        if targets:
            for f in contig.features:
                if f.type in targets:
                    write_fasta(filehandles[int(records_written / records_per_file)], f)
                    records_written += 1
        else:
            write_fasta(filehandles[int(records_written / records_per_file)], masker.mask(contig))
            records_written += 1
    for f in filehandles:
        f.close()
    if mask:
        context.log(f'({genome.id}) {masker.stats()}')
    return paths


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
                case [word, *_] if word.startswith('#'):
                    continue
                case [query, hit, percent_id, aligned_length, mismatches, gaps, query_start, query_end, hit_start,
                      hit_end, evalue, score] if 'BLAST' == self.mode:
                    hit_db_entry = self.retrieve_db_entry(hit)
                    # print(hit_db_entry)
                    b = BlastHit(query, hit_db_entry, float(percent_id), int(aligned_length),
                                    int(mismatches), int(gaps), int(query_start), int(query_end),
                                    int(hit_start), int(hit_end), float(evalue), float(score))
                    return(b)
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
                    seq = contig['seq'][start:end]
                    if strand < 0:
                        seq = reverse_complement(seq)
                    inference = self.inference if self.inference else inference
                    qualifiers = re.split(r"[=;]", qualifiers)
                    qualifiers = [unquote(str) for str in qualifiers]
                    if len(qualifiers) % 2 != 0:
                        qualifiers = qualifiers[:-1]  # this happens for example with prodigal, ending with ";"
                    qualifiers = {qualifiers[i].lower(): qualifiers[i + 1] for i in range(0, len(qualifiers), 2)}
                    feature = {'start': start,
                               'end': end,
                               'strand': strand,
                               'type': self.target_feature_type_dict[feature_type],
                               'inference': inference,
                               'seq': seq}
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


def gbk_write_feature(writer, f:SeqFeature, i=21, lw=80):
    indent = ' '*i
    location = f'{f.start+1}..{f.end}' if f.strand >= 0 else f'complement({f.start+1}..{f.end})'
    writer.write('     {:<16}{}\n'.format(f.type.name, textwrap.fill(location, width = lw, initial_indent = '',
                                                                     subsequent_indent = indent)))
    gbk_keys = {'locus_tag': f.id,
                'inference': f.inference,
                'product': f.descr,
                'taxonomy': f.taxon,
                'subsystem': ', '.join(f.subsystem),
                'notes': ', '.join(f.notes),
                'antismash': f.antismash,
                'signal_peptide': f.signal_peptide,
                'transmembrane_helixes': f.transmembrane_helixes,
                'translation': f.seq if f.type == FeatureType.CDS else ''}
    for k, v in gbk_keys.items():
        if v:
            writer.write(textwrap.fill(f'/{k}="{v}"', width = lw, initial_indent = indent, subsequent_indent = indent))
            writer.write('\n')

def gbk_write_record(writer, sr: SeqRecord, i=21, lw=80):
    writer.write('''LOCUS       {:<15} {:>12} bp    DNA              UNK 01-JAN-1980
DEFINITION  {}.
ACCESSION   {}
VERSION     {}
KEYWORDS    .
SOURCE      .
  ORGANISM  .
            .
FEATURES             Location/Qualifiers
'''.format(sr.id, len(sr), sr.id, sr.id, sr.id))
    for feature in sr.features:
        gbk_write_feature(writer, feature, i=i, lw=lw)
    writer.write('ORIGIN\n')
    lw = ((lw - 20) // 10) * 10
    print(lw)
    for i in range(0, len(sr.seq), lw):
        line_seq = sr.seq[i:i + lw].lower()
        writer.write(f'{i+1:>9}')
        for j in range(lw // 10):
            writer.write(' ')
            writer.write(line_seq[j*10:j*10+10])
        writer.write('\n')
    writer.write('//\n')


def gbk_write_genome(writer, genome: Genome, i=21, lw=80):
    for seq_rec in genome.contigs.values():
        gbk_write_record(writer, seq_rec, i=i, lw=lw)


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

