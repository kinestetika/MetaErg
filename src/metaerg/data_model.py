import re
from pathlib import Path
from enum import Enum, auto
from collections import Counter
from metaerg.run_and_read import subsystems_data


class DBentry:
    def __init__(self, *, domain: str, descr: str, taxon: str = '', ncbi: str='', gene: str='', length: int=0,
                 pos: int=0):
        self.domain = domain
        self.descr = descr
        self.taxon = taxon
        self.ncbi = ncbi
        self.gene = gene
        self.length = length
        self.pos = pos

    def __iter__(self):
        return ((k, v) for k, v in zip(('domain', 'descr', 'taxon', 'ncbi', 'gene', 'length', 'pos'),
                (self.domain, self.descr, self.taxon, self.ncbi, self.gene, self.length, self.pos)))

    def __repr__(self):
        return '{}({})'.format(type(self).__name__, ', '.join(f'{k}={v!r}' for (k, v) in self if v))

    def __len__(self):
        return self.length

    def taxon_at_genus(self) -> str:
        for t in reversed(self.taxon.split("; ")):
            if " " not in t:
                return t
        return ''


class BlastHit:
    def __init__(self, query: str, hit: DBentry, percent_id: float, aligned_length: int, mismatches: int, gaps: int,
                 query_start: int, query_end: int, hit_start: int, hit_end: int, evalue: float, score: float):
        self.query = query
        self.hit = hit
        self.percent_id = percent_id
        self.aligned_length = aligned_length
        self.mismatches = mismatches
        self.gaps = gaps
        self.query_start = query_start
        self.query_end = query_end
        self.hit_start = hit_start
        self.hit_end = hit_end
        self.evalue = evalue
        self.score = score

    def __iter__(self):
        return ((k, v) for k, v in zip(('query', 'hit', 'percent_id', 'aligned_length', 'mismatches', 'gaps',
                                        'query_start', 'query_end', 'hit_start', 'hit_end', 'evalue', 'score'),
                (self.query, self.hit, self.percent_id, self.aligned_length, self.mismatches, self.gaps,
                 self.query_start, self.query_end, self.hit_start, self.hit_end, self.evalue, self.score)))

    def __repr__(self):
        return '{}({})'.format(type(self).__name__, ', '.join(f'{k}={v!r}' for (k, v) in self if v))

    # def __repr__(self):
    #     return '{}({!r}, {!r}, {.1f}, {}, {}, {}, {}, {}, {}, {}, {.1e}, {.1f})'.format(type(self).__name__,
    #                                                                                     self.query,
    #                                                                                     self.hit,
    #                                                                                     self.percent_id,
    #                                                                                     self.aligned_length,
    #                                                                                     self.mismatches,
    #                                                                                     self.gaps,
    #                                                                                     self.query_start,
    #                                                                                     self.query_end,
    #                                                                                     self.hit_start,
    #                                                                                     self.hit_end,
    #                                                                                     self.evalue,
    #                                                                                     self.score)

    def __len__(self):
        return self.aligned_length


class BlastResult:
    def __init__(self, hits: tuple[BlastHit]):
        self.hits = hits
        if not len(hits):
            raise Exception('Attempt to create empty blast result.')

    def __iter__(self):
        return self.hits.__iter__()

    def __len__(self):
        return len(self.hits)

    def __repr__(self):
        return '{}({})'.format(type(self).__name__, ',\n'.join(f'{h!r}' for h in self))

    def query(self):
        return self.hits[0].query

    def percent_aligned(self) -> float:
        return 100 * len(self.hits[0]) / len(self.hits[0].hit)

    def percent_recall(self) -> float:
        return 100 * sum((1 for h in self.hits[1:] if h.hit.descr == self.hits[0].hit.descr)) / len(self)

    def summary(self) -> str:
        identical_function_count = sum((1 for h in self.hits[1:] if h.hit.descr == self.hits[0].hit.descr))
        return '[{}/{}] aa@{}% [{}/{}] {}'.format(len(self.hits[0]),
                                                  len(self.hits[0].hit),
                                                  self.hits[0].percent_id,
                                                  identical_function_count,
                                                  len(self),
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

    def __repr__(self):
        return '{}[{!r}]'.format(type(self).__name__, self.name)


RNA_FEATURES = (FeatureType.rRNA, FeatureType.tRNA, FeatureType.tmRNA, FeatureType.ncRNA, FeatureType.retrotransposon)


class SeqFeature:
    """Describes a sequence feature, such as a gene."""
    displayed_keys = 'start end strand type inference product taxon antismash transmembrane_helixes signal_peptide' \
                     'subsystem notes'.split()

    def __init__(self, start: int, end: int, strand: int, type, inference: str, seq: str, id: str = '', descr: str = '',
                 taxon: str = '', antismash: str = '', transmembrane_helixes: str = '', signal_peptide: str = '',
                 cdd: BlastResult = None, blast: BlastResult = None, subsystem = None, notes = None):
        self.start = start
        self.end = end
        self.strand = strand
        self.type = type if isinstance(type, FeatureType) else FeatureType[type]
        self.inference = inference
        self.seq = ''.join(seq.split())
        self.id = id
        self.descr = descr
        self.taxon = taxon
        self.antismash = antismash
        self.transmembrane_helixes = transmembrane_helixes
        self.signal_peptide = signal_peptide
        self.cdd = cdd
        self.blast = blast
        self.subsystem = subsystem if subsystem else set()
        self.notes = notes if notes else set()

    def __len__(self):
        return self.end - self.start

    def __iter__(self):
        return ((k, v) for k, v in zip(('id', 'type', 'start', 'end', 'strand', 'descr', 'notes', 'taxon', 'inference',
                                    'antismash', 'transmembrane_helixes', 'signal_peptide', 'subsystem', 'seq',
                                    'cdd', 'blast'),
                (self.id, self.type, self.start, self.end, self.strand, self.descr, self.notes, self.taxon,
                 self.inference, self.antismash, self.transmembrane_helixes, self.signal_peptide, self.subsystem,
                 self.seq, self.cdd, self.blast)))

    def __repr__(self):
        return '\n{}({})'.format(type(self).__name__, ',\n  '.join(f'{k}={v!r}' for k, v in self if v))

    def __lt__(self, other):
        return self.start < other.start

    def __gt__(self, other):
        return self.start > other.start

    def __eq__(self, other):
        return self.start == other.start

    def __le__(self, other):
        return self.start <= other.start

    def __ge__(self, other):
        return self.start >= other.start

    def __ne__(self, other):
        return self.start != other.start

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


class SubSystem:
    def __init__(self, id: str, targets: [str] = None, hits = None):
        self.id = id
        self.targets = targets if targets else list()
        self.hits = hits if hits else dict()

    def __repr__(self):
        return '{}({!r},{!r},{!r})'.format(type(self).__name__, self.id, self.targets, self.hits)

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
    def __init__(self, subsystems: dict[str, SubSystem] = None):
        self.subsystems = {}
        self.cues = {}
        current_subsystem = None
        for line in subsystems_data.subsystem_data().split('\n'):
            line = line.strip()
            if line.startswith("#") or not len(line):
                continue
            elif line.startswith(">"):
                current_subsystem = SubSystem(line[1:])
                self.subsystems[current_subsystem.id] = current_subsystem
            elif current_subsystem is not None:
                current_subsystem.targets.append(line)
                self.cues[line] = current_subsystem
        if subsystems:
            self.subsystems = subsystems

    def __repr__(self):
        return '{}({!r})'.format(type(self).__name__, self.subsystems)

    def match(self, feature: SeqFeature, descriptions):
        for d in descriptions:
            for cue, subsystem in self.cues.items():
                if len(d) > len(cue) + 20:
                    continue
                match = re.search(r'\b' + cue + r'\b', d)
                if match and match.start() < 10:
                    subsystem.add_hit(feature.id, cue)
                    feature.subsystem.add(subsystem.id)
                    return True
        return False


class SeqRecord:
    def __init__(self, id: str, seq: str, descr: str = '', features: list[SeqFeature] = None):
        self.id = id
        self.seq = ''.join(seq.split())
        self.descr = descr
        self.features = features if features else list()

    def __repr__(self):
        seq_lines = (self.seq[i:i+80] for i in range(0, len(self.seq), 80))
        return "{}(id={!r},descr={!r},features={!r},\nseq='''{}''')\n".format(type(self).__name__, self.id, self.descr,
                                                                              self.features, '\n'.join(seq_lines))

    def __len__(self):
        return len(self.seq)


class Masker:
    def __init__(self, mask=True, exceptions=None, min_length=50):
        self.apply_mask = mask
        self.exceptions = exceptions if exceptions else list()
        self.min_length = min_length
        self.nt_total = 0
        self.nt_masked = 0

    def mask(self, seq_record: SeqRecord) -> SeqRecord:
        seq = seq_record.seq
        seq_record.nt_masked = 0
        if self.apply_mask:
            for f in seq_record.features:
                if f.inference not in self.exceptions and len(f) >= self.min_length:
                    seq = seq[:f.start] + 'N' * len(f) + seq[f.end:]
                    self.nt_masked += len(f)
        self.nt_total += len(seq_record)
        return SeqRecord(id=seq_record.id, descr=seq_record.descr, seq=seq)
        # record.annotations['molecule_type'] = 'DNA'

    def stats(self):
        return f'Masked {self.nt_masked / max(self.nt_total, 1) * 100:.1f}% of sequence data.'


class Genome:
    def __init__(self, id: str, contigs: dict[str, SeqRecord]=None, delimiter: str = '.',
                 translation_table: int = 11, properties: dict = None, subsystems: SubSystems = None):
        self.id = id
        self.contigs = contigs if contigs else dict()
        self.delimiter = delimiter
        self.translation_table = translation_table
        self.properties = properties if properties else dict()
        self.subsystems = subsystems if subsystems else SubSystems()

    def __len__(self):
        return sum(len(c) for c in self.contigs.values())

    def __repr__(self):
        return '{}(id={!r},\ndelimiter={!r},\ntranslation_table={!r},\n' \
               'properties={!r},\nsubsystems={!r},\ncontigs={!r})\n'.format(type(self).__name__,
                                                                         self.id,
                                                                         self.delimiter,
                                                                         self.translation_table,
                                                                         self.properties,
                                                                         self.subsystems,
                                                                         self.contigs)

    def validate_ids(self):
        if self.delimiter in self.id:
            raise Exception(f'Genome id {self.id} contains {self.delimiter}; change using --delimiter')
        for c_id in self.contigs.keys():
            if self.delimiter in c_id:
                raise Exception(f'Contig id {c_id} contains {self.delimiter}; change using --delimiter')

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
                f.id = self.delimiter.join((self.id, c.id, f'{f_id:05d}'))
                f_id += 1

    def get_feature(self, feature_id):
        id = feature_id.split(self.delimiter)
        return self.contigs[id[1]].features[int(id[2])]

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