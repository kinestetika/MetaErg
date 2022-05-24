import re
import copy
from enum import Enum, auto
from dataclasses import dataclass, field, InitVar
from collections import Counter, Sequence
from pathlib import Path
from typing import NamedTuple

from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio import SeqUtils
from Bio import SeqIO
from BCBio import GFF

from metaerg.run_and_read import subsystems_data


class DBentry(NamedTuple):
    id: str
    gene: str
    descr: str
    taxon: str
    length: int
    pos: int

    def taxon_at_genus(self) -> str:
        for t in reversed(self.taxon.split("; ")):
            if " " not in t:
                return t
        return ''


class BlastHit(NamedTuple):
    query: str
    hit: DBentry
    percent_id: float
    aligned_length: int
    mismatches: int
    gaps: int
    query_start: int
    query_end: int
    hit_start: int
    hit_end: int
    evalue: float
    score: float


class BlastResult(NamedTuple):
    hits: tuple[BlastHit]

    def query(self):
        return self.hits[0].query

    def percent_aligned(self) -> float:
        return 100 * self.hits[0].aligned_length / self.hits[0].hit.length

    def percent_recall(self) -> float:
        return 100 * sum((1 for h in self.hits[1:] if h.hit.descr == self.hits[0].hit.descr)) / len(self.hits)

    def summary(self) -> str:
        identical_function_count = sum((1 for h in self.hits[1:] if h.hit.descr == self.hits[0].hit.descr))
        return '[{}/{}] aa@{}% [{}/{}] {}'.format(self.hits[0].aligned_length,
                                                  self.hits[0].hit.length,
                                                  self.hits[0].percent_id,
                                                  identical_function_count,
                                                  len(self.hits),
                                                  self.hits[0].hit.descr)

class TabularBlastParser:
    def __init__(self, filename, mode, retrieve_db_entry):
        self.filename = filename
        self.mode = mode
        self.next_hit = None
        self.current_query = None
        self.retrieve_db_entry = retrieve_db_entry

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


NON_IUPAC_RE = re.compile(r'[^ACTGN]')


class FeatureType(Enum):
    CDS = auto()
    rRNA = auto()
    tRNA = auto()
    tmRNA = auto()
    ncRNA = auto()
    repeat = auto()
    crispr_repeat = auto()
    retrotransposon = auto()


@dataclass(order=True)
class MetaergSeqFeature:
    """Describes a sequence feature, such as a gene."""
    start: int  # as start is the first field, it will cause features to be ordered by their position
    end: int
    strand: int
    type: FeatureType
    inference: str
    sequence: str
    id: str = ''
    product: str = ''
    taxon: str = ''
    antismash: str = ''
    transmembrane_helixes: str = ''
    signal_peptide: str = ''
    cdd: BlastResult = field(init=False)
    blast: BlastResult = field(init=False)
    subsystem: set[str] = field(default_factory=set)
    notes: set[str] = field(default_factory=set)
    # because it does not have a type, exported_keys becomes a class attribute
    exported_keys = 'id sequence inference product taxon antismash transmembrane_helixes signal_peptide'.split()

    def __len__(self):
        return self.end - self.start

    def tmh_count(self):
        try:
            return int(self.transmembrane_helixes.split()[0])
        except ValueError:
            return 0

    def taxon_at_genus(self) -> str:
        for t in reversed(self.taxon.split("; ")):
            if " " not in t:
                return t
        return ''

    def make_biopython_feature(self) -> SeqFeature:
        """Returns a BioPython SeqFeature with this content"""
        loc = FeatureLocation(self.start, self.end, self.strand)
        qal = {key: getattr(self, key) for key in MetaergSeqFeature.exported_keys if getattr(self, key)}
        return SeqFeature(location=loc, type=self.type, qualifiers=qal)

    def make_biopython_record(self) -> SeqRecord:
        return SeqRecord(Seq(self.sequence), id=self.id, description=f'{self.product} [{self.taxon}]')


@dataclass()
class SubSystem:
    id: str
    targets: list[str] = field(default_factory=list, init=False)
    hits: dict[str] = field(default_factory=dict, init=False)

    def add_hit(self, feature_id: str, target: str = 'none'):
        self.hits.setdefault(feature_id, set()).add(target)

    def get_hits(self, target):
        return (k for k, v in self.hits.items() if target in v)

    def get_stats(self):
        if self.targets:
            genes_present = len(set(self.hits.values()))
            return genes_present, len(self.targets), genes_present / len(self.targets)
        else:
            return len(self.hits), 0, 1


class SubSystems:
    def __init__(self):
        self.subsystems = {}
        self.cues = {}
        current_subsystem = None
        for line in subsystems_data.subsystem_data().split('\n'):
            line = line.strip()
            if line.startswith("#") or not len(line):
                continue
            if line.startswith(">"):
                current_subsystem = SubSystem(line[1:])
                self.subsystems[current_subsystem.id] = current_subsystem
                continue
            current_subsystem.targets.append(line)
            self.cues[line] = current_subsystem

    def match(self, feature: MetaergSeqFeature, descriptions):
        for d in descriptions:
            for cue, subsystem in self.cues.items():
                if len(d.descr) > len(cue) + 20:
                    continue
                match = re.search(r'\b' + cue + r'\b', d)
                if match and match.start() < 10:
                    subsystem.add_hit(feature.id, cue)
                    feature.subsystem.add(subsystem.id)
                    return True
        return False


@dataclass()
class MetaergSeqRecord:
    id: str
    sequence: str
    description: str = ''
    translation_table: int = 11
    features: list[MetaergSeqFeature] = field(init=False, default_factory=list)
    nt_masked: field(init=False) = 0

    def __post_init__(self):  # clean up sequence
        self.sequence = NON_IUPAC_RE.sub('N', self.sequence)
        self.sequence = ''.join(self.sequence.split()).upper()

    def __len__(self):
        return len(self.sequence)

    def make_biopython_record(self, make_features=True, masked=False, mask_excep=None, min_mask_length=50) -> SeqRecord:
        """Returns a BioPython SeqRecord with this content"""
        seq = self.sequence
        self.nt_masked = 0
        if masked:
            for f in self.features:
                if f.inference not in mask_excep and len(f) >= min_mask_length:
                    self.nt_masked += len(f)
                    seq = seq[:f.start] + 'N' * len(f) + seq[f.end:]
        record = SeqRecord(Seq(seq), id=self.id, description=self.description)
        record.annotations['molecule_type'] = 'DNA'
        if make_features:
            for f in self.features:
                record.features.append(f.make_biopython_feature())
        return record

    def spawn_feature(self, start: int, end: int, strand: int, type: FeatureType, **kwargs) -> MetaergSeqFeature:
        assert strand == 1 or strand == -1, f'invalid value {strand} for strand, needs to be 1 or -1.'
        assert end > start, f'invalid coordinates, end needs to be greater than start.'
        assert 0 <= start < len(self), f'start coordinate out of range.'
        assert 0 <= end < len(self), f'end coordinate out of range.'
        seq = Seq(self.sequence[start:end:strand])
        notes = kwargs.get("notes", set())
        if FeatureType.CDS == type:
            remainder = len(seq) % 3
            if remainder:
                seq = seq + Seq('N' * (3 - remainder))
            if '*' in seq:
                notes.add('contains internal stop codon(s).')
            seq = str(seq.translate(self.translation_table)[:-1])
        f = MetaergSeqFeature(start, end, strand, type, sequence=seq, notes=notes, **kwargs)
        self.features.append(f)
        return f

    def spawn_features(self, features):
        for f in features():
            if isinstance(f, MetaergSeqFeature):
                f_copy = copy.deepcopy(f)
                self.features.append(f_copy)
            elif isinstance(f, SeqFeature):
                self.spawn_feature(f.location.start, f.location.end, f.location.strand, FeatureType[f.type],
                                   **f.qualifiers)


@dataclass()
class MetaergGenome:
    id: str
    contig_file: InitVar[Path]
    translation_table: int = 11
    delimiter: str = '.'
    contigs: dict[str, MetaergSeqRecord] = field(default_factory=dict)
    properties: dict = field(init=False, default_factory=dict)
    subsystems: SubSystems = field(init=False)
    contig_name_mappings: dict[str, str] = field(init=False, default_factory=dict)
    rename_contigs: InitVar[bool] = True
    min_contig_length: InitVar[int] = 500

    def __post_init__(self, contig_file, rename_contigs, min_contig_length):
        """Reads contig fasta file, sorts by len, filters for min len, whitespace, and funny chars"""
        assert self.delimiter not in self.id, f'({self.id}) Genome name contains "{self.delimiter}",' \
                                              f' rename or change delimiter.'
        self.subsystems = SubSystems()
        if contig_file:
            i = 0
            c: SeqRecord
            for c in sorted([SeqIO.parse(contig_file, "fasta")], key=len, reverse=True):
                if len(c) < min_contig_length:
                    break
                if rename_contigs:
                    new_id = f'{self.id}.c{i:0>4}'
                    self.contig_name_mappings[c.id] = new_id
                else:
                    assert self.delimiter not in c.id, \
                        f'({self.id}) Contig name {c.id} contain "{self.delimiter}", rename or change delimiter.'
                    new_id = c.id
                seq = NON_IUPAC_RE.sub('N', str(c.seq))
                seq = ''.join(seq.split()).upper()

                contig = MetaergSeqRecord(new_id, seq, c.description, self.translation_table)
                self.contigs[contig.id] = contig
                contig.spawn_features(c.features)
                i += 1

    def __len__(self):
        return sum(len(c) for c in self.contigs.values())

    def generate_feature_ids(self):
        f_id = 0
        for c in self.contigs.values():
            c.features.sort()
            for f in c.features:
                f.id = self.delimiter.join((self.id, c.id, f'{f_id::05d}'))
                f_id += 1

    def get_feature(self, feature_id):
        id = feature_id.split(self.delimiter)
        return self.contigs[id[1]].features[int(id[2])]

    def write_fasta_files(self, base_file: Path, split=1, target=None, **kwargs_masking):
        """writes features (of target FeatureType), or contigs (target = None), to one or more (split) fasta files,
        optionally masking features with N"""
        if target:
            number_of_records = sum(1 for c in self.contigs.values() for f in c.features if f.type == target
                                    or (isinstance(target, Sequence) and f.type in target))
        else:
            number_of_records = len(self.contigs)
        split = min(split, number_of_records)
        records_per_file = number_of_records / split
        paths = (Path(base_file.parent, f'{base_file.name}.{i}') for i in range(split)) if split > 1 else base_file,
        filehandles = [open(p, 'w') for p in paths]
        number_of_records = 0
        nt_masked = 0

        for contig in self.contigs.values():
            if target:
                for f in contig.features:
                    if f.type == target or (isinstance(target, Sequence) and f.type in target):
                        SeqIO.write(f.make_biopython_record(), filehandles[int(number_of_records / records_per_file)],
                                    "fasta")
                        number_of_records += 1
            else:
                seq_record = contig.make_biopython_record(make_features=False, **kwargs_masking)
                SeqIO.write(seq_record, filehandles[int(number_of_records / records_per_file)], "fasta")
                number_of_records += 1
                nt_masked += contig.nt_masked
        for f in filehandles:
            f.close()
        if kwargs_masking['masked']:
            exec.log(f'Masked {nt_masked / len(self) * 100:.1f}% of sequence data.')
        if split > 1:
            return paths
        else:
            return base_file

    def write_gbk_gff(self, gbk_file=None, gff_file=None):
        record_generator = (c.make_biopython_record() for c in self.contigs.values())
        if gbk_file:
            SeqIO.write(record_generator, gbk_file, "genbank")
        if gff_file:
            with open(gff_file, "w") as gff_handle:
                GFF.write(record_generator, gff_handle)

    def write_contig_name_mappings(self, mappings_file):
        with open(mappings_file, 'w') as mapping_writer:
            for key, value in self.contig_name_mappings.items():
                mapping_writer.write(f'{key}\t{value}\n')

    def compute_properties(self):
        self.properties['size'] = len(self)
        self.properties['percent GC'] = int(sum((len(contig) * SeqUtils.GC(contig.sequence) for contig in
                                                 self.contigs.values())) / self.properties['size'] + 0.5)
        cum_size = 0
        for contig in sorted(self.contigs.values(), key=len, reverse=True):
            cum_size += len(contig)
            if cum_size >+ self.properties['size'] / 2:
                self.properties["N50"] = len(contig)
                break
        self.properties['#proteins'] = sum(1 for contig in self.contigs.values() for f in contig.features
                                           if f.type == FeatureType.CDS)
        self.properties['percent coding'] = int(sum(len(f) for contig in self.contigs.values() for f in contig.features
                                                if f.type == FeatureType.CDS) / self.properties['size'] * 100 + 0.5)
        self.properties['mean protein length (aa)'] = int(self.properties['percent coding'] * self.properties['size']
                                                          / 3 / self.properties['#proteins'])
        self.properties['#ribosomal RNA'] = sum(1 for contig in self.contigs.values() for f in contig.features
                                                if f.type == FeatureType.rRNA)
        self.properties['#transfer RNA'] = sum(1 for contig in self.contigs.values() for f in contig.features
                                               if f.type == FeatureType.tRNA)
        self.properties['#non coding RNA'] = sum(1 for contig in self.contigs.values() for f in contig.features
                                                 if f.type == FeatureType.ncRNA)
        self.properties['#retrotransposons'] = sum(1 for contig in self.contigs.values() for f in contig.features
                                                   if f.type == FeatureType.retrotransposon)
        self.properties['#CRISPR repeats'] = sum(1 for contig in self.contigs.values() for f in contig.features
                                                 if f.type == FeatureType.crispr_repeat)
        self.properties['#other repeats'] = sum(1 for contig in self.contigs.values() for f in contig.features
                                                if f.type == FeatureType.repeat)
        self.properties['percent repeats'] = int(100 * sum(len(f) for contig in self.contigs.values() for f in
                                                           contig.features if f.type in (FeatureType.repeat,
                                                           FeatureType.retrotransposon, FeatureType.crispr_repeat))
                                                 / self.properties['size'] + 0.5)
        self.properties['total # features'] = sum(len(contig.features) for contig in self.contigs.values())

        taxon_counts = Counter()
        taxon_counts.update(f.taxon for contig in self.contigs.values() for f in contig.features)
        dominant_taxon, highest_count = taxon_counts.most_common(1)[0]
        self.properties['dominant taxon'] = f'{dominant_taxon} ({highest_count/sum(taxon_counts.values()) * 100:.1f}%)'
        return self.properties


def pad_seq(sequence):
    """ Pad sequence to multiple of 3 with N """
    remainder = len(sequence) % 3
    return sequence if remainder == 0 else sequence + Seq('N' * (3 - remainder))