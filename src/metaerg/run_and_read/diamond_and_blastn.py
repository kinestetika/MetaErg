from pathlib import Path
from collections import namedtuple
from metaerg.run_and_read.data_model import MetaergSeqFeature
from metaerg.run_and_read import abc
from metaerg import utils
from metaerg import subsystems

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

    def __repr__(self):
        return f'DiamondAndBlastN({self.genome}, {self.exec})'

    def _purpose(self) -> str:
        """Should return the purpose of the tool"""
        return 'function prediction and taxonomic classification of proteins and RNA genes with diamond and blastn'

    def _programs(self) -> tuple:
        """Should return a tuple with the programs needed"""
        return 'diamond', 'blastn'

    def _databases(self) -> tuple:
        """Should return a tuple with database files needed"""
        return self.db_descr_file, self.db_taxon_file

    def _result_files(self) -> tuple:
        """Should return a tuple with the result files (Path objects) created by the programs"""
        return self.diamond_file, self.blastn_file

    def _run_programs(self):
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

    def _process_blast_result(self, blast_result, subsystem_cues):
        feature: MetaergSeqFeature = self.genome.get_feature(blast_result.query)
        top_hit: utils.BlastHit = blast_result.hits[0]
        top_entry = self.get_db_entry(top_hit.hit)
        self.genome.blast_results['blast'][blast_result.query] = blast_result.hits
        identical_function_count = 1
        for hit in blast_result.hits[1:]:
            db_entry = self.get_db_entry(hit.hit)
            if db_entry.description == top_entry.description:
                identical_function_count += 1
            if hit.aligned_length / db_entry.length >= 0.8:
                cue, subsystem_name = subsystems.match_subsystem(db_entry.description, subsystem_cues)
                if subsystem_name:
                    feature.subsystem.add(subsystem_name)
                    self.genome.subsystems[subsystem_name].add_hit(feature.id, cue)
        feature.taxonomy = top_entry.taxon
        feature.description = '[{}/{}] aa@{}% [{}/{}] {}'.format(top_hit.aligned_length,
                                                                 top_entry.length,
                                                                 top_hit.percent_id,
                                                                 identical_function_count,
                                                                 len(blast_result.hits),
                                                                 top_entry.description)

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
        self.genome.blast_results['blast'] = {}
        subsystem_cues = {}
        for subsystem in self.genome.subsystems.values():
            subsystem_cues |= subsystem.get_cues_as_hash()
        with utils.TabularBlastParser(self.diamond_file, 'BLAST') as handle:
            for blast_result in handle:
                blast_result_count += 1
                self._process_blast_result(blast_result, subsystem_cues)
        with utils.TabularBlastParser(self.blastn_file, 'BLAST') as handle:
            for blast_result in handle:
                blast_result_count += 1
                self._process_blast_result(blast_result, subsystem_cues)
        return blast_result_count
