import re
from pathlib import Path
from enum import Enum, auto
from dataclasses import dataclass, field
from collections import Counter
from typing import NamedTuple
from metaerg.run_and_read import subsystems_data


class DBentry(NamedTuple):
    domain: str
    taxon: str
    descr: str
    ncbi: str
    gene: str
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


class FeatureType(Enum):
    CDS = auto()
    rRNA = auto()
    tRNA = auto()
    tmRNA = auto()
    ncRNA = auto()
    repeat = auto()
    crispr_repeat = auto()
    retrotransposon = auto()


RNA_FEATURES = (FeatureType.rRNA, FeatureType.tRNA, FeatureType.tmRNA, FeatureType.ncRNA, FeatureType.retrotransposon)


@dataclass(order=True)
class MetaergSeqFeature:
    """Describes a sequence feature, such as a gene."""
    start: int  # as start is the first field, it will cause features to be ordered by their position
    end: int
    strand: int
    type: FeatureType
    inference: str
    seq: str
    id: str = ''
    descr: str = ''
    taxon: str = ''
    antismash: str = ''
    transmembrane_helixes: str = ''
    signal_peptide: str = ''
    cdd: BlastResult = None
    blast: BlastResult = None
    subsystem: set[str] = field(default_factory=set)
    notes: set[str] = field(default_factory=set)
    # because it does not have a type, exported_keys becomes a class attribute
    exported_keys = 'id sequence inference product taxon antismash transmembrane_helixes signal_peptide subsystem ' \
                    'notes'.split()
    displayed_keys = 'start end strand type inference product taxon antismash transmembrane_helixes signal_peptide' \
                     'subsystem notes'.split()

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
    seq: str
    descr: str = ''
    features: list[MetaergSeqFeature] = field(init=False, default_factory=list)

    def __len__(self):
        return len(self.seq)


class Masker:
    def __init__(self, mask=True, exceptions=None, min_length=50):
        self.mask = mask
        self.exceptions = exceptions
        self.min_length = min_length
        self.nt_total = 0
        self.nt_masked = 0

    def mask(self, seq_record) -> MetaergSeqRecord:
        seq = seq_record.seq
        seq_record.nt_masked = 0
        if self.mask:
            for f in seq_record.features:
                if f.inference not in self.exceptions and len(f) >= self.min_length:
                    seq_record.nt_masked += len(f)
                    seq = seq[:f.start] + 'N' * len(f) + seq[f.end:]
                    self.nt_masked += len(f)
        self.nt_total += len(seq_record)
        return MetaergSeqRecord(seq_record.id, seq_record.descr, seq)
        # record.annotations['molecule_type'] = 'DNA'

    def stats(self):
        return f'Masked {self.nt_masked / self.nt_total * 100:.1f}% of sequence data.'


@dataclass()
class MetaergGenome:
    id: str
    contigs: dict[str, MetaergSeqRecord] = field(default_factory=dict)
    delimiter: str = '.'
    translation_table: int = 11
    properties: dict = field(init=False, default_factory=dict)
    subsystems: SubSystems = field(init=False)

    def __post_init__(self):
        self.subsystems = SubSystems()

    def validate_ids(self):
        assert self.delimiter not in self.id, f'Genome id {self.id} may not contain delimiter {self.delimiter}!'
        for c_id in self.contigs.keys():
            assert self.delimiter not in c_id, f'Contig id {c_id} may not contain delimiter {self.delimiter}!'

    def __len__(self):
        return sum(len(c) for c in self.contigs.values())

    def rename_contigs(self, mappings_file:Path):
        i = 0
        with open(mappings_file, 'w') as mapping_writer:
            for c in self.contigs.values():
                new_id = f'{self.id}.c{i:0>4}'
                mapping_writer.write(f'{c.id}\t{new_id}\n')
                c.id = new_id
                i += 1

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

    # def write_gbk_gff(self, gbk_file=None, gff_file=None):
    #     record_generator = (c.make_biopython_record() for c in self.contigs.values())
    #     if gbk_file:
    #         SeqIO.write(record_generator, gbk_file, "genbank")
    #     if gff_file:
    #         with open(gff_file, "w") as gff_handle:
    #             GFF.write(record_generator, gff_handle)

    def compute_properties(self):
        self.properties['size'] = len(self)
        self.properties['percent GC'] = int(sum((c.seq.count('G') + c.seq.count('G') for c in
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