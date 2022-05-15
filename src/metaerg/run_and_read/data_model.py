import re
import ast
from pathlib import Path
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from metaerg import utils
from Bio import SeqIO


class MetaergSeqFeature:
    def __init__(self, parent, type='', location=None, inference='', seq_feature: SeqFeature = None, gff_line=None):
        assert parent, 'Attempt to construct MetaergSeqFeature without a parent'
        self.id = ''
        self.parent = parent
        self.sequence = ''
        self.description = ''
        self.location = location
        self.type = type
        self.inference = inference
        self.taxonomy = ''
        self.note = ''
        self.warning = ''
        self.cdd = ''
        self.antismash = ''
        self.canthyd_function = ''
        self.transmembrane_helixes = ''
        self.signal_peptide = ''
        self.subsystem = []
        if seq_feature:
            self._clone_biopython_feature(seq_feature)
        elif gff_line:
            self._clone_gff_line(gff_line)

    def __repr__(self):
        return f'MetaergSeqFeature(parent={self.parent}, id={self.id}, type={self.type}, location={self.location})'

    def __len__(self):
        return len(self.location)

    def make_sequence(self) -> str:
        """Sets and returns self.sequence using location and parent sequence"""
        biopython_parent = self.parent.make_biopython_record(make_features=False)
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
        qal = {key: getattr(self, key) for key in dir(self) if not key.startswith('_')
               and key not in ('type', 'location', 'make_biopython_feature', 'make_biopython_record')}
        return SeqFeature(location=self.location, type=self.type, qualifiers=qal)

    def make_biopython_record(self) -> SeqRecord:
        if not self.sequence:
            self.make_sequence()
        return SeqRecord(Seq(self.sequence), id=self.id, description=self.description)

    def _clone_biopython_feature(self, seq_feature: SeqFeature):
        """Initializes the content from a BioPython SeqFeature"""
        self.location = seq_feature.location
        self.type = seq_feature.type
        for key, value in seq_feature.qualifiers.items():
            setattr(self, key, value)
        self.subsystem = ast.literal_eval(self.subsystem)

    def _clone_gff_line(self, line: str):
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
            if key in ('parent', 'location'):
                continue
            try:
                setattr(self, key, value)
            except AttributeError:
                continue


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
        self.subsystems = {}
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
