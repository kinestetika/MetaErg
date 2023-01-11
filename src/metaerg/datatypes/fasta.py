import gzip
import re
import numpy as np
from pathlib import Path
from metaerg import context
from metaerg.datatypes import sqlite

NON_IUPAC_RE_NT = re.compile(r'[^ACTGN]')
NON_IUPAC_RE_AA = re.compile(r'[^RHKDESTNQCUGPAVILMFYW]')
ALL_MASK_TARGETS = set('CDS rRNA tRNA tmRNA ncRNA repeat crispr_repeat retrotransposon'.split())
COMPLEMENT = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A',}


def reverse_complement(seq:str) -> str:
    return ''.join(COMPLEMENT.get(b, 'N') for b in reversed(seq))


def _parse_fasta_header(line:str) -> dict:
    si = line.find(' ')
    if si > 0:
        return {'id': line[1:si], 'seq': '', 'descr': line[si + 1:].strip()}
    else:
        return {'id': line[1:].strip(), 'seq': '', 'descr': ''}


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
            if line.startswith('#'):
                continue
            if line.startswith('>'):
                if len(seq) and seq_rec is not None:
                    seq_rec['seq'] = self._cleanup(''.join(seq))
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
            seq_rec['seq'] = self._cleanup(''.join(seq))
            yield seq_rec

    def _cleanup(self, seq) -> str:
        if self.cleanup_seq:
            seq = seq.upper()
            seq = ''.join(seq.split())
            if self.alphabet:
                seq = self.alphabet.sub(self.unknown_char, seq)
        return seq


class Masker:
    def __init__(self, db_connection, targets=None):
        self.db_connection = db_connection
        self.targets = targets
        self.nt_total = 0
        self.nt_masked = 0

    def mask(self, seq_record):
        seq = seq_record['seq']
        if self.targets:
            for feature in sqlite.read_all_features(self.db_connection,  contig=seq_record['id'], type=self.targets):
                seq = seq[:feature.start] + 'N' * feature.length_nt() + seq[feature.end:]
                self.nt_masked += feature.length_nt()

        self.nt_total += len(seq)

        masked_record = seq_record.copy()
        masked_record['seq'] = seq
        return masked_record

    def stats(self):
        return f'Masked {self.nt_masked} / {self.nt_total} nt ({self.nt_masked / max(self.nt_total, 1):.1%}).'


def write_fasta(handle, fasta, line_length=80):
    if not fasta:
        raise Exception('Attempt to write zero-length sequence to fasta.')
    handle.write(f'>{fasta["id"]}')
    try:
        if fasta["descr"] and fasta["descr"] is not np.nan:
            handle.write(f' {fasta["descr"]}')
    except KeyError:
        pass
    try:
        if fasta["subsystems"]:
            handle.write(f' ({fasta["subsystems"]})')
    except KeyError:
        pass
    try:
        if fasta["taxon"]:
            handle.write(f' ({fasta["taxon"]})')
    except KeyError:
        pass
    handle.write('\n')
    for i in range(0, len(fasta['seq']), line_length):
        handle.write(fasta['seq'][i:i+line_length])
        handle.write('\n')


def write_features_to_fasta(db_connection, seq_type: str, base_file: Path, split=1, targets = None):
    number_of_records = sum(1 for f in sqlite.read_all_features(db_connection, type=targets))
    split = min(split, number_of_records)
    records_per_file = number_of_records / split if split else number_of_records
    if split > 1:
        paths = [Path(base_file.parent, f'{base_file.name}.{i}') for i in range(split)]
    else:
        paths = base_file,
    filehandles = [open(p, 'w') for p in paths]
    records_written = 0
    for feature in sqlite.read_all_features(db_connection, type=targets):
        fasta_rec = {'id': feature.id, 'descr': feature.descr}
        if 'aa' == seq_type:
            fasta_rec['seq'] = feature.aa_seq
        elif 'nt' == seq_type:
            fasta_rec['seq'] = feature.nt_seq
        write_fasta(filehandles[int(records_written / records_per_file)], fasta_rec)
        records_written += 1
    if not number_of_records:
        context.log(f'WARNING: No [{", ".join(targets)}] features in genome, '
                    f'fasta file(s) {", ".join([str(p.name) for p in paths])} are empty!')
    return paths


def write_contigs_to_fasta(contig_dict: dict, base_file: Path, db_connection, genome_name='',
                           split=1, mask_targets=None):
    """writes contigs to fasta file(s), optionally masking features with N"""
    number_of_records = len(contig_dict)
    split = min(split, number_of_records)
    records_per_file = number_of_records / split
    if split > 1:
        paths = [Path(base_file.parent, f'{base_file.name}.{i}') for i in range(split)]
    else:
        paths = base_file,
    filehandles = [open(p, 'w') for p in paths]
    records_written = 0
    masker = Masker(db_connection, targets=mask_targets)
    for contig in contig_dict.values():
        write_fasta(filehandles[int(records_written / records_per_file)], masker.mask(contig))
        records_written += 1
    for f in filehandles:
        f.close()
    if mask_targets:
        context.log(f'({genome_name}) {masker.stats()}')
    return paths
