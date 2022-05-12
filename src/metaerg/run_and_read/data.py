import re
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from metaerg import utils

class MetaergSeqFeature:
    def __init__(self, parent, seq_feature:SeqFeature=None, gff_line=None):
        assert parent, 'Attempt to construct MetaergSeqFeature without a parent'
        self.id = ''
        self.parent = parent
        self.sequence = ''
        self.description = ''
        self.location = None
        self.type = ''
        self.inference = ''
        self.taxonomy = ''
        self.note = ''
        self.warning = ''
        self.cdd = ''
        self.antismash_region = ''
        self.antismash_region_number = 0
        self.antismash_function = ''
        self.antismash_category = ''
        self.canthyd_function = ''
        self.transmembrane_helixes = ''
        self.tmh_topology = ''
        self.signal_peptide = ''
        self.subsystem = ''
        if seq_feature:
            self._clone_biopython_feature(seq_feature)
        elif gff_line:
            self._clone_gff_line()

    def __repr__(self):
        return f'MetaergSeqFeature(parent={self.parent}, id={self.id}, type={self.type}, location={self.location})'

    def __len__(self):
        return len(self.location)

    def make_sequence(self) -> str:
        """Sets and returns self.sequence using location and parent sequence"""
        biopython_parent = self.parent.make_biopython_record(make_features = False)
        biopython_feature = self.make_biopython_feature()
        if 'CDS' == self.type:
            seq = biopython_feature.extract(biopython_parent)
            remainder = len(seq) % 3
            if remainder:
                seq = seq + Seq('N' * (3 - remainder))
            self.sequence = str(seq.translate(table=self.parent.parent.translation_table)[:-1])
            if '*' in self.sequence:
                self.warning = 'internal stop codon(s) in CDS'
        else:
            self.sequence = str(biopython_feature.extract(biopython_parent))
        return self.sequence

    def make_biopython_feature(self) -> SeqFeature:
        """Returns a BioPython SeqFeature with this content"""
        qal = { key: getattr(self, key) for key in dir(self) if not key.startswith('_')
                and not key in ('type', 'location', 'make_biopython_feature', 'make_biopython_record')}
        return SeqFeature(location=self.location, type=self.type, qualifiers=qal)

    def make_biopython_record(self) -> SeqRecord:
        if not self.sequence:
            self.make_sequence()
        return SeqRecord(Seq(self.sequence), id=self.id, description=self.description)

    def _clone_biopython_feature(self, seq_feature:SeqFeature):
        """Initializes the content from a BioPython SeqFeature"""
        self.location = seq_feature.location
        self.type = seq_feature.type
        for key, value in seq_feature.qualifiers.items():
            setattr(self, key, value)

    def _clone_gff_line(self, line:str):
        """Initializes the content from a line of a .gff file"""
        words = line.split('\t')
        self.inference = words[1]
        self.type = words[2]
        strand = None
        if '+' == words[6]:
            strand = +1
        elif '-' == words[6]:
            strand = -1
        self.location = FeatureLocation(int(words[3]) - 1, int(words[4]), strand=strand)

        qal = re.split(r"[=;]", words[8])
        qal = [utils.unescape_str(str) for str in qal]
        if len(qal) % 2 != 0:
            qal = qal[:-1]  # this happens for example with prodigal which has the qualifier column ending with ";"
        qal = {qal[i].lower(): qal[i + 1] for i in range(0, len(qal), 2)}
        for key, value in qal:
            try:
                setattr(self, key, value)
            except AttributeError:
                continue


class MetaergSeqRecord:
    def __init__(self, parent, seq_record:SeqRecord=None):
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

    def _clone_biopython_record(self, seq_record:SeqRecord):
        """Initializes the content from a BioPython SeqRecord"""
        self.id = seq_record.id
        self.description = seq_record.description
        self.sequence = seq_record.seq
        for f in seq_record.features:
            self.features.append(MetaergSeqFeature(f))


class MetaergGenome:
    def __init__(self, name, contig_dict:dict, translation_table=11):
        self.name = name
        self.translation_table = translation_table
        self.contigs = {}
        for name, c in contig_dict:
            self.contigs[name] = MetaergSeqRecord(self, seq_record=c)
