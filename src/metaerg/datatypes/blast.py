import gzip


class DBentry:
    FIELDS = ('domain', 'descr', 'taxon', 'accession', 'gene', 'length', 'pos', 'min_score', 'min_t_score')

    def __init__(self, *, domain: str, descr: str, taxon: str = '', accession: str = '', gene: str = '',
                 length: int = 0, pos: int = 0, min_score: int = 0, min_t_score: int = 0):
        self.domain = domain
        self.descr = descr
        self.taxon = taxon
        self.accession = accession
        self.gene = gene
        self.length = length
        self.pos = pos
        self.min_score = min_score
        self.min_t_score = min_t_score

    def __iter__(self):
        return ((k, v) for k, v in zip(DBentry.FIELDS, (self.domain, self.descr, self.taxon, self.accession,
                                                        self.gene, self.length, self.pos, self.min_score,
                                                        self.min_t_score)))

    def __repr__(self):
        return '{}({})'.format(type(self).__name__, ', '.join(f'{k}={v!r}' for (k, v) in self))

    def __len__(self):
        return self.length


def taxon_at_genus(taxon) -> str:
    if type(taxon) == str:
        for t in reversed(taxon.split("; ")):
            if " " not in t:
                return t
    return ''


class BlastHit:
    FIELDS = ('query', 'hit', 'percent_id', 'aligned_length', 'mismatches', 'gaps', 'query_start',
              'query_end', 'hit_start', 'hit_end', 'evalue', 'score')

    def __init__(self, query: str, hit: DBentry, percent_id: float, aligned_length: int,
                 query_start: int, query_end: int, hit_start: int, hit_end: int, evalue: float, score: float,
                 mismatches: int = 0, gaps: int = 0):
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
        return ((k, v) for k, v in zip(BlastHit.FIELDS, (self.query, self.hit, self.percent_id, self.aligned_length,
                                                         self.mismatches, self.gaps, self.query_start, self.query_end,
                                                         self.hit_start, self.hit_end, self.evalue, self.score)))

    def __repr__(self):
        return '{}({})'.format(type(self).__name__, ', '.join(f'{k}={v!r}' for (k, v) in self))

    def __len__(self):
        return max(self.aligned_length, 1)


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
        return '{}(({},))'.format(type(self).__name__, ',\n'.join(f'{h!r}' for h in self))

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
            words = line.strip().split()  # \t works for blast, but hmmer and cmscan use spaces
            match words:
                case [word, *_] if word.startswith('#'):
                    continue
                case [query, hit, percent_id, aligned_length, mismatches, gaps, query_start, query_end, hit_start,
                      hit_end, evalue, score] if 'BLAST' == self.mode:
                    hit_db_entry = self.retrieve_db_entry(hit)
                    b = BlastHit(query, hit_db_entry, float(percent_id), int(aligned_length), int(query_start),
                                 int(query_end), int(hit_start), int(hit_end), float(evalue), float(score),
                                 int(mismatches), int(gaps))
                    return(b)
                case [hit, _, query, _, evalue, score, _, _, _, _, _, _, _, _, _, _, _, _, *_] if 'HMMSCAN' == self.mode:
                    hit_db_entry = self.retrieve_db_entry(hit)
                    return BlastHit(query, hit_db_entry, 0, hit_db_entry.length, 0, 0, 0, 0, float(evalue), float(score), 0, 0)
                case [hit, _, _, query, _, _, evalue, score, _, _, _, _, _, _, _, hit_start, hit_end,
                      query_start, query_end, *_] if 'HMMSCAN_DOM_TABLE' == self.mode:
                    hit_db_entry = self.retrieve_db_entry(hit)
                    hit_start = int(hit_start)
                    hit_end = int(hit_end)
                    return BlastHit(query, hit_db_entry, 0, hit_end - hit_start, int(query_start), int(query_end),
                                    hit_start, hit_end, float(evalue), float(score))

                case [query, _, hit, _, evalue, score, _, _, _, _, _, _, _, _, _, _, _, _, *_] if 'HMMSEARCH' == self.mode:
                    hit_db_entry = self.retrieve_db_entry(hit)
                    return BlastHit(query, hit_db_entry, 0, hit_db_entry.length, 0, 0, 0, 0, float(evalue), float(score), 0, 0)
                case [*_]:
                    continue
        return None
