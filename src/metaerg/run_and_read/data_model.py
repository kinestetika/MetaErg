import re
import ast
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
    ncRNA = auto()
    repeat = auto()
    crispr_repeat = auto()


BlastHit = namedtuple('BlastHit', ['query', 'hit', 'percent_id', 'aligned_length', 'mismatches', 'gaps',
                                   'query_start', 'query_end', 'hit_start', 'hit_end', 'evalue', 'score', 'db_entry'])


@dataclass
class MetaergSeqFeature:
    """Describes a sequence feature, such as a gene."""
    type: FeatureType
    start: int
    end: int
    strand: int
    id: str
    description: str
    inference: str
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
        seq = Seq(parent_sequence[self.start:self.end:self.strand])
        if FeatureType.CDS == self.type:
            remainder = len(seq) % 3
            if remainder:
                seq = seq + Seq('N' * (3 - remainder))
            if '*' in seq:
                self.notes.add('contains internal stop codon(s).')
            seq = seq.translate(table=translation_table)[:-1]
        self.sequence = str(seq)

    def __len__(self):
        return self.start - self.end

    def get_description(self) -> str:
        if self.blast:
            identical_function_count = sum((1 for h in self.blast[1:] if h.db_entry.description
                                            == self.blast[0].db_entry.description))
            return '[{}/{}] aa@{}% [{}/{}] {}'.format(self.blast[0].aligned_length,
                                                      self.blast[0].db_entry.length,
                                                      self.blast[0].percent_id,
                                                      identical_function_count,
                                                      len(self.blast),
                                                      self.blast[0].db_entry.description)
        else:
            return ''


    def make_biopython_feature(self) -> SeqFeature:
        """Returns a BioPython SeqFeature with this content"""
        loc = FeatureLocation(self.start, self.end, self.strand)
        qal = {key: getattr(self, key) for key in
               'id sequence inference antismash transmembrane_helixes signal_peptide subsystem notes'.split()}
        qal['description'] = self.get_description()

        return SeqFeature(location=loc, type=self.type, qualifiers=qal)

    def make_biopython_record(self) -> SeqRecord:
        return SeqRecord(Seq(self.get_sequence()), id=self.id, description=self.description)

    def _clone_biopython_feature(self, seq_feature: SeqFeature):
        """Initializes the content from a BioPython SeqFeature"""
        self.location = seq_feature.location
        self.type = seq_feature.type
        for key, value in seq_feature.qualifiers.items():
            setattr(self, key, value)
        self.subsystem = ast.literal_eval(self.subsystem)


class MetaergSeqRecord:
    def __init__(self, parent, seq_record: SeqRecord = None):
        assert parent, 'Attempt to construct MetaergSeqRecord without a parent'
        self.id = ''
        self.description = ''
        self.sequence = ''
        self.parent = parent
        self.features = list()
        if seq_record:
            self._clone_biopython_record(seq_record)

    def __repr__(self):
        return f'MetaergSeqRecord(id={self.id}, {len(self.features)} features)'

    def __len__(self):
        return len(self.sequence)

    def make_biopython_record(self, make_features=True) -> SeqRecord:
        """Returns a BioPython SeqRecord with this content"""
        record = SeqRecord(Seq(self.sequence), id=self.id, description=self.description)
        if make_features:
            for f in self.features:
                record.features.append(f.make_biopython_feature())
        return record

    def _clone_biopython_record(self, seq_record: SeqRecord):
        """Initializes the content from a BioPython SeqRecord"""
        self.id = seq_record.id
        self.description = seq_record.description
        self.sequence = seq_record.seq
        for f in seq_record.features:
            self.features.append(MetaergSeqFeature(f))

    def mask_seq(self, exceptions=None, min_mask_length=50) -> (SeqRecord, int):
        seq = self.sequence
        nt_masked = 0
        for f in self.features:
            if f.inference in exceptions or len(f.location) < min_mask_length:
                continue
            fl: FeatureLocation = f.location
            nt_masked += len(fl)
            seq = seq[:fl.start] + 'N' * (fl.end - fl.start) + seq[fl.end:]
        return SeqRecord(Seq(seq), id=self.id, description=self.description), nt_masked

    def spawn_feature(self, type, location, inference) -> MetaergSeqFeature:
        f = MetaergSeqFeature(self, type, location, inference)
        self.features.append(f)
        return f


class MetaergGenome:
    def __init__(self, name, contig_dict: dict, translation_table=11):
        self.name = name
        self.translation_table = translation_table
        self.contigs = {}
        self.subsystems = subsystems.init_subsystems()
        self.properties = {}
        self.blast_results = {}
        for name, c in contig_dict:
            self.contigs[name] = MetaergSeqRecord(self, seq_record=c)

    def __repr__(self):
        return f'MetaergGenome(name={self.name}, {len(self.contigs)} contigs)'

    def get_feature(self, feature_id):
        id = feature_id.split('.')
        return self.contigs[id[1]].features[int(id[2])]

    def make_masked_contig_fasta_file(self, masked_fasta_file, exceptions=None, min_mask_length=50) -> Path:
        if exceptions is None:
            exceptions = set()
        (self.nt_masked, self.nt_total) = (0, 0)
        seq_iterator = (record.mask_seq(exceptions=exceptions, min_mask_length=min_mask_length)
                        for record in self.contigs.values())
        SeqIO.write(seq_iterator, masked_fasta_file, "fasta")
        utils.log(f'Masked {self.nt_masked / self.nt_total * 100:.1f}% of sequence data.')
        return masked_fasta_file

    def make_split_fasta_files(self, base_file: Path, number_of_files, target) -> ():
        if 'CDS' == target:
            number_of_records = sum(1 for contig in self.contigs.values() for f in contig.features if f.type == 'CDS')
        else:
            number_of_records = len(self.contigs)
        number_of_files = min(number_of_files, number_of_records)
        seqs_per_file = number_of_records / number_of_files
        paths = (Path(base_file.parent, f'{base_file.name}.{i}') for i in range(number_of_files))
        filehandles = [open(p, 'w') for p in paths]
        number_of_records = 0
        for contig in self.contigs.values():
            if 'CDS' == target:
                for f in contig.features:
                    if f.type == 'CDS':
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
        utils.log(f'({self.name}) Compiling genome stats...')
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
                                           if f.type == 'CDS')
        self.properties['percent coding'] = sum(len(f) for contig in self.contigs.values() for f in contig.features
                                                if f.type == 'CDS') / self.properties['size']
        self.properties['mean protein length (aa)'] = int(self.properties['percent coding'] * self.properties['size']
                                                          / 3 / self.properties['#proteins'])
        self.properties['#ribosomal RNA'] = sum(1 for contig in self.contigs.values() for f in contig.features
                                                if f.type == 'rRNA')
        self.properties['#transfer RNA'] = sum(1 for contig in self.contigs.values() for f in contig.features
                                               if f.type == 'tRNA')
        self.properties['#non coding RNA'] = sum(1 for contig in self.contigs.values() for f in contig.features
                                                 if f.type == 'tRNA')
        self.properties['#retrotransposons'] =  sum(1 for contig in self.contigs.values() for f in contig.features
                                                    if f.type == 'retrotransposon')
        self.properties['#CRISPR repeats'] =  sum(1 for contig in self.contigs.values() for f in contig.features
                                                    if f.type == 'crispr_repeat')
        self.properties['#other repeats'] =  sum(1 for contig in self.contigs.values() for f in contig.features
                                                 if f.type == 'repeat region')
        self.properties['percent repeats'] = int(100 * sum(len(f) for contig in self.contigs.values() for f in
                                                           contig.features if f.type
                                                           in {'repeat region', 'retrotransposon', 'crispr_repeat'}) \
                                                 / self.properties['size'] + 0.5)
        self.properties['total # features'] = sum(1 for contig in self.contigs.values() for f in contig.features)

        taxon_counts = Counter()
        taxon_counts.update(f.taxonomy for contig in self.contigs.values() for f in contig.features)
        dominant_taxon, highest_count = taxon_counts.most_common(1)[0]
        self.properties['dominant taxon'] = f'{dominant_taxon} ({highest_count/sum(taxon_counts.values()) * 100:.1f}%)'

        utils.log(f'({self.name}) Compilation of stats complete...')
        return self.properties
