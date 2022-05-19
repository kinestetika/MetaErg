from metaerg.run_and_read.data_model import BlastHit, BlastResult

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
            match(words):
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
                    self.next_hit = BlastHit(query, hit, 0, 0, 0, 0, 0, 0, 0, 0, float(evalue), float(score))
                case [query, _, hit, _, evalue, score, _, _, _, _, _, _, _, _, _, _, _, _, *_] if 'HMMSEARCH' == self.mode:
                    hit_db_entry = self.retrieve_db_entry(hit)
                    self.next_hit = BlastHit(query, hit, 0, 0, 0, 0, 0, 0, 0, 0, float(evalue), float(score))
                case [*_]:
                    continue
            return True

    def __next__(self):
        all_hits: list[BlastHit] = []
        if self.next_hit:
            all_hits.append(self.next_hit)
        while self.load_next_hit_from_file():
            if self.current_query != self.next_hit.query:
                prev_query = self.current_query
                self.current_query = self.next_hit.query
                return BlastResult(prev_query, tuple(all_hits))
            all_hits.append(self.next_hit)
        if len(all_hits):
            return BlastResult(self.current_query, tuple(all_hits))
        raise StopIteration
