from pathlib import Path
from metaerg.run_and_read.data_model import MetaergSeqFeature, BlastResult, DBentry, TabularBlastParser
from metaerg.run_and_read.abc import Annotator, ExecutionEnvironment, register
from metaerg import utils


@register
class DiamondAndBlastN(Annotator):
    def __init__(self, genome, exec_env: ExecutionEnvironment):
        super().__init__(genome, exec_env)
        self.diamond_file = self.spawn_file('diamond')
        self.blastn_file = self.spawn_file('blastn')
        self.db_descr_file = Path(self.exec.database_dir, 'db_descriptions.txt')
        self.db_taxon_file = Path(self.exec.database_dir, 'db_taxonomy.txt')
        self.taxonomy = {}  # this is associated with the blast database
        self.descriptions = {}  # this is associated with the blast database
        self._pipeline_position = 81
        self._purpose = 'function prediction and taxonomic classification of proteins and RNA genes with diamond and blastn'
        self._programs = ('diamond', 'blastn')
        self._databases = (self.db_descr_file, self.db_taxon_file)
        self._result_files = (self.diamond_file, self.blastn_file)

    def _run_programs(self):
        """Executes the helper programs to complete the analysis"""
        cds_aa_file = self.spawn_file('cds.faa')
        rna_nt_file = self.spawn_file('rna.nt')
        utils.run_external(f'diamond blastp -d {Path(self.exec.database_dir, "db_protein.faa")} -q {cds_aa_file} '
                           f'-o {self.diamond_file} -f 6  --threads {self.exec.cpus_per_genome} --max-target-seqs 10')
        utils.run_external(f'blastn -db {Path(self.exec.database_dir, "db_rna.fna")} -query {rna_nt_file} '
                           f'-out {self.blastn_file} -max_target_seqs 10 -outfmt 6')

    def get_db_entry(self, db_id):
        words = db_id.split('~')  # org_acc gene_acc [pev] gene# decr# taxon#
        return DBentry(db_id, '', self.descriptions[words[2]][int(words[4])], self.taxonomy[words[2]][int(words[5])],
                       int(words[6]), int(words[3]))

    def _process_blast_result(self, blast_result: BlastResult):
        feature: MetaergSeqFeature = self.genome.get_feature(blast_result.query())
        feature.blast = blast_result
        feature.product = blast_result.summary()
        self.genome.subsystems.match(feature, (h.hit.descr for h in blast_result.hits
                                               if h.aligned_length / h.hit.length >= 0.8))

    def _read_results(self) -> int:
        """Should parse the result files and return the # of positives."""
        # (1) load databse descriptions
        if self.db_descr_file.exists():
            with open(self.db_descr_file) as descr_handle:
                for line in descr_handle:
                    words = line.split('\t')
                    self.descriptions[words[0]][int(words[1])] = words[2].strip()
            utils.log(
                f'Parsed ({len(self.descriptions["p"])}, {len(self.descriptions["e"])}, {len(self.descriptions["v"])}) '
                f'gene descriptions from db for (prokaryotes, eukaryotes and viruses) respectively. ')
        # (2) load database taxonomy
        if self.db_taxon_file.exists():
            with open(self.db_taxon_file) as taxon_handle:
                for line in taxon_handle:
                    words = line.split('\t')
                    self.taxonomy[words[0]][int(words[1])] = words[2].strip().replace('~', '; ')
            utils.log(
                f'Parsed ({len(self.taxonomy["p"])}, {len(self.taxonomy["e"])}, {len(self.taxonomy["v"])}) '
                f'taxa from db for (prokaryotes, eukaryotes and viruses) respectively.')
        # (3) parse diamond blast results
        blast_result_count = 0
        with TabularBlastParser(self.diamond_file, 'BLAST', self.get_db_entry) as handle:
            for blast_result in handle:
                blast_result_count += 1
                self._process_blast_result(blast_result)
        with TabularBlastParser(self.blastn_file, 'BLAST', self.get_db_entry) as handle:
            for blast_result in handle:
                blast_result_count += 1
                self._process_blast_result(blast_result)
        return blast_result_count
