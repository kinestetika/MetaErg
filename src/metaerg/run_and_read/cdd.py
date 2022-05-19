from pathlib import Path
import shutil
from collections import namedtuple
from concurrent.futures import ProcessPoolExecutor
from metaerg.run_and_read.data_model import MetaergSeqFeature
from metaerg.run_and_read import abc
from metaerg import utils
from metaerg import subsystems


class CDD(abc.Annotator):
    def __init__(self, genome, exec_env: abc.ExecutionEnvironment):
        super().__init__(genome, exec_env)
        self.cdd_file = self.spawn_file('cdd')
        self.db_cdd_index = Path(self.exec.database_dir, 'cddid.tbl')
        self.db_cdd = Path(self.exec.database_dir, "cdd", "Cdd")
        self.cdd = {}  # this is the cdd index

    def __repr__(self):
        return f'CDD({self.genome}, {self.exec})'

    def _purpose(self) -> str:
        """Should return the purpose of the tool"""
        return 'function prediction using RPSBlast and the conserved domain database'

    def _programs(self) -> tuple:
        """Should return a tuple with the programs needed"""
        return 'rpsblast',

    def _databases(self) -> tuple:
        """Should return a tuple with database files needed"""
        return self.db_cdd_index, self.db_cdd

    def _result_files(self) -> tuple:
        """Should return a tuple with the result files (Path objects) created by the programs"""
        return self.cdd_file,

    def _run_programs(self):
        """Should execute the helper programs to complete the analysis"""
        cds_aa_file = self.spawn_file('cds.faa')
        if self.exec.threads > 1:
            split_fasta_files = self.genome.write_fasta_files(cds_aa_file, self.exec.threads, target='CDS')
            split_cdd_files = [Path(self.cdd_file.parent, f'{self.cdd_file.name}.{i}')
                               for i in range(len(split_fasta_files))]
            with ProcessPoolExecutor(max_workers=self.exec.threads) as executor:
                for split_input, split_output in zip(split_fasta_files, split_cdd_files):
                    executor.submit(utils.run_external, f'rpsblast -db {self.db_cdd} -query {split_input} '
                                                        f'-out {split_output} -outfmt 6 -evalue 1e-7')

            with open(self.cdd_file, 'wb') as output:
                for split_input_file, split_output_file in zip(split_fasta_files, split_cdd_files):
                    with open(split_output_file, 'rb') as input:
                        shutil.copyfileobj(input, output)
                    split_input_file.unlink()
                    split_output_file.unlink()
        else:
            utils.run_external(f'rpsblast -db {self.db_cdd} -query {cds_aa_file} '
                               f'-out {self.cdd_file} -outfmt 6 -evalue 1e-7')

    def _read_results(self) -> int:
        """Should parse the result files and return the # of positives."""
        # load cdd index
        if self.db_cdd_index.exists():
            with open(self.db_cdd_index) as db_handle:
                for line in db_handle:
                    words = line.split("\t")
                    self.cdd[int(words[0])] = DBEntry(words[1], words[2], words[3], int(words[4]))
            utils.log(f'Parsed {len(self.cdd)} entries from conserved domain database.')
        # parse cdd results
        cdd_result_count = 0
        self.genome.blast_results['blast'] = {}  # html modules will need this later on
        subsystem_cues = {}
        for subsystem in self.genome.subsystems.values():
            subsystem_cues |= subsystem.get_cues_as_hash()
        with utils.TabularBlastParser(self.cdd_file, 'BLAST') as handle:
            for blast_result in handle:
                feature: MetaergSeqFeature = self.genome.get_feature(blast_result.query)
                self.genome.blast_results['cdd'][blast_result.query] = blast_result.hits
                for h in blast_result[1]:
                    cdd_result_count += 1
                    cdd_item: DBEntry = self.cdd[int(h.hit[4:])]
                    cdd_descr = f'{cdd_item.name}|{cdd_item.gene} {cdd_item.descr}'
                    if not feature.cdd:
                        feature.cdd = '[{}/{}]@{:.1f}% [{}-{}] {}'.format(h.aligned_length, cdd_item.length,
                                                                          h.percent_id,
                                                                          h.query_start, h.query_end,
                                                                          cdd_descr)
                    if h.aligned_length / cdd_item.length >= 0.8:
                        cue, subsystem_name = subsystems.match_subsystem(cdd_descr, subsystem_cues)
                        if subsystem_name:
                            feature.subsystem.add(subsystem_name)
                            self.genome.subsystems[subsystem_name].add_hit(feature.id, cue)
        return cdd_result_count
