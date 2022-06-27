import gzip

class DBentry:
    def __init__(self, *, domain: str, descr: str, taxon: str = '', ncbi: str = '', gene: str = '', length: int = 0,
                 pos: int = 0):
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
        return 100 * min(1.0, len(self.hits[0]) / len(self.hits[0].hit))

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