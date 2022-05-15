from pathlib import Path
from collections import namedtuple
from metaerg.run_and_read.data_model import MetaergSeqFeature
from metaerg.run_and_read import abc
from metaerg import utils

DBEntry = namedtuple('DiamondEntry', ['taxon', 'description', 'length', 'pos'])


class DiamondAndBlastN(abc.AbstractBaseClass):
    def __init__(self, genome, exec_env: abc.ExecutionEnvironment):
        super().__init__(genome, exec_env)
        self.diamond_file = self.spawn_file('diamond')
        self.blastn_file = self.spawn_file('blastn')
        self.db_descr_file = Path(self.exec.database_dir, 'db_descriptions.txt')
        self.db_taxon_file = Path(self.exec.database_dir, 'db_taxonomy.txt')
        self.taxonomy = {}  # this is associated with the blast database
        self.descriptions = {}  # this is associated with the blast database
        self.feature_hits = {}  # these are all the blast results

    def __repr__(self):
        return f'DiamondAndBlastN({self.genome}, {self.exec})'

    def __purpose__(self) -> str:
        """Should return the purpose of the tool"""
        return 'function prediction and taxonomic classification of proteins and RNA genes with diamond and blastn'

    def __programs__(self) -> tuple:
        """Should return a tuple with the programs needed"""
        return 'diamond', 'blastn'

    def __databases__(self) -> tuple:
        """Should return a tuple with database files needed"""
        return self.db_descr_file, self.db_taxon_file

    def __result_files__(self) -> tuple:
        """Should return a tuple with the result files (Path objects) created by the programs"""
        return self.diamond_file, self.blastn_file

    def __run_programs__(self):
        """Should execute the helper programs to complete the analysis"""
        cds_aa_file = self.spawn_file('cds.faa')
        rna_nt_file = self.spawn_file('rna.nt')
        utils.run_external(f'diamond blastp -d {Path(self.exec.database_dir, "db_protein.faa")} -q {cds_aa_file} '
                           f'-o {self.diamond_file} -f 6  --threads {self.exec.threads} --max-target-seqs 10')
        utils.run_external(f'blastn -db {Path(self.exec.database_dir, "db_rna.fna")} -query {rna_nt_file} '
                           f'-out {self.blastn_file} -max_target_seqs 10 -outfmt 6')

    def get_db_entry(self, db_id):
        words = db_id.split('~')  # org_acc gene_acc [pev] gene# decr# taxon#
        return DBEntry(self.taxonomy[words[2]][int(words[5])], self.descriptions[words[2]][int(words[4])],
                            int(words[6]), int(words[3]))

    def _process_blast_result(self, blast_result):
        feature: MetaergSeqFeature = self.genome.get_feature(blast_result.query)
        top_hit: utils.BlastHit = blast_result.hits[0]
        top_entry = self.get_db_entry(top_hit.hit)
        self.feature_hits[blast_result.query] = blast_result.hits
        identical_function_count = 1
        for hit in blast_result.hits[1:]:
            db_entry = self.get_db_entry(hit.hit)
            if db_entry.description == top_entry.description:
                identical_function_count += 1
        feature.taxonomy = top_entry.taxon
        feature.description = '[{}/{}] aa@{}% [{}/{}] {}'.format(top_hit.aligned_length,
                                                                 top_entry.length,
                                                                 top_hit.percent_id,
                                                                 identical_function_count,
                                                                 len(blast_result.hits),
                                                                 top_entry.description)

    def __read_results__(self) -> int:
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
                    self.taxonomy[words[0]][int(words[1])] = words[2].strip().replace('~', '~ ')
            utils.log(
                f'Parsed ({len(self.taxonomy["p"])}, {len(self.taxonomy["e"])}, {len(self.taxonomy["v"])}) '
                f'taxa from db for (prokaryotes, eukaryotes and viruses) respectively.')
        # (3) parse diamond blast results
        blast_result_count = 0
        with utils.TabularBlastParser(self.diamond_file, 'BLAST') as handle:
            for blast_result in handle:
                blast_result_count += 1
                self._process_blast_result(blast_result)
        with utils.TabularBlastParser(self.blastn_file) as handle:
            for blast_result in handle:
                blast_result_count += 1
                self._process_blast_result(blast_result)
        return blast_result_count
