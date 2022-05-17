from enum import Enum, auto
from dataclasses import dataclass, field, InitVar
from collections import Counter
from pathlib import Path
from collections import namedtuple

from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio import SeqUtils
from Bio import SeqIO

import subsystems
from metaerg import utils

class FeatureType(Enum):
    CDS = auto()
    rRNA = auto()
    tRNA = auto()
    tmRNA = auto()
    ncRNA = auto()
    repeat = auto()
    crispr_repeat = auto()
    retrotransposon = auto()


BlastHit = namedtuple('BlastHit', ['query', 'hit', 'percent_id', 'aligned_length', 'mismatches', 'gaps',
                                   'query_start', 'query_end', 'hit_start', 'hit_end', 'evalue', 'score', 'db_entry'])


@dataclass(order=True)
class MetaergSeqFeature:
    """Describes a sequence feature, such as a gene."""
    start: int  # as start is the first field, it will cause features to be ordered by their position
    end: int
    strand: int
    type: FeatureType
    inference: str
    id: str = ''
    description: str = ''
    cdd: list[BlastHit] = field(default_factory=list)
    blast: list[BlastHit] = field(default_factory=list)
    antismash: str = ''
    transmembrane_helixes: str = ''
    signal_peptide: str = ''
    subsystem: set[str] = field(default_factory=set)
    notes: set[str] = field(default_factory=set)
    sequence: str = None
    parent_sequence: InitVar[str] = None
    translation_table: InitVar[int] = 11

    def __post_init__(self, parent_sequence, translation_table):
        """Sets and returns self.sequence using location and parent sequence"""
        assert self.strand == 1 or self.strand == -1, f'invalid value {self.strand} for strand, needs to be 1 or -1.'
        assert self.end > self.start, f'invalid coordinates, end needs to be greater than start.'
        assert 0 <= self.start < len(parent_sequence), f'start coordinate out of range.'
        assert 0 <= self.end < len(parent_sequence), f'end coordinate out of range.'
        if parent_sequence:
            seq = Seq(parent_sequence[self.start:self.end:self.strand])
            if FeatureType.CDS == self.type:
                remainder = len(seq) % 3
                if remainder:
                    seq = seq + Seq('N' * (3 - remainder))
                if '*' in seq:
                    self.notes.add('contains internal stop codon(s).')
                seq = seq.translate(table=translation_table)[:-1]
            self.sequence = str(seq)
        else:
            assert len(self.sequence), 'No sequence or parent sequence provided to MetaergseqFeature'

    def __len__(self):
        return self.end - self.start

    def get_description(self) -> str:
        if len(self.blast):
            identical_function_count = sum((1 for h in self.blast[1:] if h.db_entry.description
                                            == self.blast[0].db_entry.description))
            return '[{}/{}] aa@{}% [{}/{}] {}'.format(self.blast[0].aligned_length,
                                                      self.blast[0].db_entry.length,
                                                      self.blast[0].percent_id,
                                                      identical_function_count,
                                                      len(self.blast),
                                                      self.blast[0].db_entry.description)
        else:
            return self.description

    def get_taxon(self) -> str:
        if len(self.blast):
            return self.blast[0].db_entry.taxon
        else:
            return ''

    def make_biopython_feature(self) -> SeqFeature:
        """Returns a BioPython SeqFeature with this content"""
        loc = FeatureLocation(self.start, self.end, self.strand)
        qal = {key: getattr(self, key) for key in
               'id sequence inference antismash transmembrane_helixes signal_peptide subsystem notes'.split()
               if getattr(self, key)}
        qal['description'] = self.get_description()
        qal['taxon'] = self.get_taxon()
        return SeqFeature(location=loc, type=self.type, qualifiers=qal)

    def make_biopython_record(self) -> SeqRecord:
        return SeqRecord(Seq(self.sequence), id=self.id, description=f'{self.get_description()} [{self.get_taxon()}]')


@dataclass()
class MetaergSeqRecord:
    id: str
    sequence: str
    description: str = ''
    translation_table: int = 11
    features: list[MetaergSeqFeature] = field(init=False, default_factory=list)

    def __len__(self):
        return len(self.sequence)

    def make_biopython_record(self, make_features=True) -> SeqRecord:
        """Returns a BioPython SeqRecord with this content"""
        record = SeqRecord(Seq(self.sequence), id=self.id, description=self.description)
        if make_features:
            for f in self.features:
                record.features.append(f.make_biopython_feature())
        return record

    def mask_seq(self, exceptions=None, min_mask_length=50) -> (SeqRecord, int):
        seq = self.sequence
        nt_masked = 0
        for f in self.features:
            if f.inference in exceptions or len(f) < min_mask_length:
                continue
            nt_masked += len(f)
            seq = seq[:f.start] + 'N' * len(f) + seq[f.end:]
        return SeqRecord(Seq(seq), id=self.id, description=self.description), nt_masked

    def spawn_feature(self, start:int, end:int, strand:int, type:FeatureType, inference:str) -> MetaergSeqFeature:
        f = MetaergSeqFeature(start, end, strand, type, inference, parent_sequence=self.sequence,
                              translation_table=self.translation_table)
        self.features.append(f)
        return f


@dataclass()
class MetaergGenome:
    id: str
    translation_table: int = 11
    contigs: dict[str, MetaergSeqRecord] = field(default_factory=dict)
    properties: dict = field(init=False, default_factory=dict)
    subsystems: list[subsystems.SubSystem] = field(init=False, default_factory=list)
    contig_dict: InitVar[dict] = None

    def __post_init__(self, contig_dict):
        self.subsystems = subsystems.init_subsystems()
        if contig_dict:
            for c in contig_dict.values():
                contig = MetaergSeqRecord(c.id, c.sequence, c.description, self.translation_table)
                self.contigs[contig.id] = contig
                for f in c.features:
                    feature = contig.spawn_feature(f.location.start, f.location.end, f.location.strand,
                                                   FeatureType[f.type], inference='')
                    for key, value in f.qualifiers.items():
                        if hasattr(feature, key):
                            setattr(feature, key, value)


    def get_feature(self, feature_id):
        id = feature_id.split('.')
        return self.contigs[id[1]].features[int(id[2])]

    def make_masked_contig_fasta_file(self, masked_fasta_file, exceptions=None, min_mask_length=50) -> Path:
        if exceptions is None:
            exceptions = set()
        nt_masked, nt_total = 0, 0
        seq_iterator = (record.mask_seq(exceptions=exceptions, min_mask_length=min_mask_length)
                        for record in self.contigs.values())
        SeqIO.write(seq_iterator, masked_fasta_file, "fasta")
        utils.log(f'Masked {nt_masked / nt_total * 100:.1f}% of sequence data.')
        return masked_fasta_file

    def make_split_fasta_files(self, base_file: Path, number_of_files, target:FeatureType = None) -> ():
        if FeatureType.CDS == target:
            number_of_records = sum(1 for c in self.contigs.values() for f in c.features if f.type == FeatureType.CDS)
        else:
            number_of_records = len(self.contigs)
        number_of_files = min(number_of_files, number_of_records)
        seqs_per_file = number_of_records / number_of_files
        paths = (Path(base_file.parent, f'{base_file.name}.{i}') for i in range(number_of_files))
        filehandles = [open(p, 'w') for p in paths]
        number_of_records = 0

        for contig in self.contigs.values():
            if FeatureType.CDS == target:
                for f in contig.features:
                    if f.type == FeatureType.CDS:
                        SeqIO.write(f.make_biopython_record(), filehandles[int(number_of_records / seqs_per_file)],
                                    "fasta")
                        number_of_records += 1
            else:
                SeqIO.write(contig.make_biopython_record(False), filehandles[int(number_of_records / seqs_per_file)],
                            "fasta")
                number_of_records += 1
        for f in filehandles:
            f.close()
        return paths

    def compute_properties(self):
        utils.log(f'({self.id}) Compiling genome stats...')
        self.properties['size'] = sum((len(contig) for contig in self.contigs.values()))
        self.properties['percent GC'] = int(sum((len(contig) * SeqUtils.GC(contig.sequence) for contig in
                                                 self.contigs.values())) / self.properties['size'] + 0.5)
        cum_size = 0
        for contig in sorted(self.contigs.values(), key=len, reverse=True):
            cum_size += len(contig)
            if cum_size > self.properties['size'] / 2:
                self.properties["N50"] = len(contig)
                break
        self.properties['#proteins'] = sum(1 for contig in self.contigs.values() for f in contig.features
                                           if f.type == FeatureType.CDS)
        self.properties['percent coding'] = sum(len(f) for contig in self.contigs.values() for f in contig.features
                                                if f.type == FeatureType.CDS) / self.properties['size']
        self.properties['mean protein length (aa)'] = int(self.properties['percent coding'] * self.properties['size']
                                                          / 3 / self.properties['#proteins'])
        self.properties['#ribosomal RNA'] = sum(1 for contig in self.contigs.values() for f in contig.features
                                                if f.type == FeatureType.rRNA)
        self.properties['#transfer RNA'] = sum(1 for contig in self.contigs.values() for f in contig.features
                                               if f.type == FeatureType.tRNA)
        self.properties['#non coding RNA'] = sum(1 for contig in self.contigs.values() for f in contig.features
                                                 if f.type == FeatureType.ncRNA)
        self.properties['#retrotransposons'] =  sum(1 for contig in self.contigs.values() for f in contig.features
                                                    if f.type == FeatureType.retrotransposon)
        self.properties['#CRISPR repeats'] =  sum(1 for contig in self.contigs.values() for f in contig.features
                                                    if f.type == FeatureType.crispr_repeat)
        self.properties['#other repeats'] =  sum(1 for contig in self.contigs.values() for f in contig.features
                                                 if f.type == FeatureType.repeat)
        self.properties['percent repeats'] = int(100 * sum(len(f) for contig in self.contigs.values() for f in
                                                           contig.features if f.type in (FeatureType.repeat,
                                                            FeatureType.retrotransposon, FeatureType.crispr_repeat))
                                                 / self.properties['size'] + 0.5)
        self.properties['total # features'] = sum(1 for contig in self.contigs.values() for f in contig.features)

        taxon_counts = Counter()
        taxon_counts.update(f.get_taxon() for contig in self.contigs.values() for f in contig.features)
        dominant_taxon, highest_count = taxon_counts.most_common(1)[0]
        self.properties['dominant taxon'] = f'{dominant_taxon} ({highest_count/sum(taxon_counts.values()) * 100:.1f}%)'

        utils.log(f'({self.id}) Compilation of stats complete...')
        return self.properties
