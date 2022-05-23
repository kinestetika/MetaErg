from pathlib import Path
import shutil
from concurrent.futures import ProcessPoolExecutor
from metaerg.run_and_read.data_model import MetaergSeqFeature, FeatureType, TabularBlastParser, DBentry
from metaerg.run_and_read.abc import Annotator, ExecutionEnvironment, register
from metaerg import utils


@register
class CDD(Annotator):
    def __init__(self, genome, exec_env: ExecutionEnvironment):
        super().__init__(genome, exec_env)
        self.cdd_file = self.spawn_file('cdd')
        self.db_cdd_index = Path(self.exec.database_dir, 'cddid.tbl')
        self.db_cdd = Path(self.exec.database_dir, "cdd", "Cdd")
        self.cdd = {}  # this is the cdd index
        self._pipeline_position = 71
        self._purpose = 'function prediction using RPSBlast and the conserved domain database'
        self._programs = ('rpsblast',)
        self._databases = (self.db_cdd_index, self.db_cdd)
        self._result_files = (self.cdd_file,)

    def _run_programs(self):
        """Executes the helper programs to complete the analysis"""
        cds_aa_file = self.spawn_file('cds.faa')
        if self.exec.cpus_per_genome > 1:
            split_fasta_files = self.genome.write_fasta_files(cds_aa_file, self.exec.cpus_per_genome, target=FeatureType.CDS)
            split_cdd_files = [Path(self.cdd_file.parent, f'{self.cdd_file.name}.{i}')
                               for i in range(len(split_fasta_files))]
            with ProcessPoolExecutor(max_workers=self.exec.cpus_per_genome) as executor:
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

    def get_cdd_entry(self, db_id) ->DBentry:
        return self.cdd[int(db_id[4:])]

    def _read_results(self) -> int:
        """Should parse the result files and return the # of positives."""
        # load cdd index
        if self.db_cdd_index.exists():
            with open(self.db_cdd_index) as db_handle:
                for line in db_handle:
                    words = line.split("\t")
                    self.cdd[int(words[0])] = DBentry(words[0], words[1], words[2], '', int(words[4]), 0)
            utils.log(f'Parsed {len(self.cdd)} entries from conserved domain database.')
        # parse cdd results
        cdd_result_count = 0
        with TabularBlastParser(self.cdd_file, 'BLAST', self.get_cdd_entry) as handle:
            for blast_result in handle:
                feature: MetaergSeqFeature = self.genome.get_feature(blast_result.query)
                cdd_result_count += 1
                feature.cdd = blast_result
                self.genome.subsystems.match(feature, (h.hit.descr for h in blast_result.hits
                                                       if h.aligned_length / h.hit.length >= 0.8))
                top_entry = blast_result.hits[0].hit
                feature.product = f'{top_entry.id}|{top_entry.gene} {top_entry.descr}'
                if len(feature.product) > 35:
                        feature.product = feature.product[:35] + '...'
        return cdd_result_count
