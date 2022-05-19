import shutil
from collections import namedtuple

from Bio import SeqIO

from metaerg.run_and_read.data_model import MetaergSeqFeature
from metaerg.run_and_read import abc
from metaerg import utils


DBEntry = namedtuple('CDDEntry', ['name', 'gene', 'descr', 'length'])


class Antismash(abc.Annotator):
    def __init__(self, genome, exec_env: abc.ExecutionEnvironment):
        super().__init__(genome, exec_env)
        self.antismash_file = self.spawn_file('antismash')

    def __repr__(self):
        return f'Antismash({self.genome}, {self.exec})'

    def _purpose(self) -> str:
        """Should return the purpose of the tool"""
        return 'prediction of secondary metabolite genes with antismash'

    def _programs(self) -> tuple:
        """Should return a tuple with the programs needed"""
        return 'antismash',

    def _result_files(self) -> tuple:
        """Should return a tuple with the result files (Path objects) created by the programs"""
        return self.antismash_file,

    def _run_programs(self):
        """Should execute the helper programs to complete the analysis"""
        gbk_file = self.spawn_file('gbk')
        if self.antismash_file.exists():
            shutil.rmtree(self.antismash_file)
        utils.run_external(f'antismash --genefinding-tool none --output-dir {self.antismash_file} {gbk_file}')

    def _read_results(self) -> int:
        """Should parse the result files and return the # of positives."""
        antismash_hit_count = 0
        for f in sorted(self.antismash_file.glob("*region*.gbk")):
            with open(f) as handle:
                antismash_region_name = ''
                antismash_region_number = 0
                for antismash_record in SeqIO.parse(handle, "genbank"):
                    for antismash_feature in antismash_record.features:
                        antismash_hit_count += 1
                        if 'region' == antismash_feature.type:
                            antismash_region_name = utils.get_feature_qualifier(antismash_feature, "rules")
                            antismash_region_number = int(utils.get_feature_qualifier(antismash_feature, "region_number"))
                        elif 'CDS' in antismash_feature.type:
                            id = utils.decipher_metaerg_id(utils.get_feature_qualifier(antismash_feature, "locus_tag"))
                            feature: MetaergSeqFeature = self.genome.get_feature(id)
                            antismash_gene_function = utils.get_feature_qualifier(antismash_feature, "gene_functions")
                            antismash_gene_category =  utils.get_feature_qualifier(antismash_feature, "gene_kind")
                            if antismash_region_name:
                                feature.antismash = ' '.join((f'(region {antismash_region_number})',
                                                             antismash_region_name,
                                                             antismash_gene_function,
                                                             antismash_gene_category))
                                feature.subsystem.add('[secondary-metabolites]')
                                self.genome.subsystems['[secondary-metabolites]'].add_hit(feature.id)
        if not antismash_hit_count:
            self.antismash_file.mkdir(exist_ok=True)  # to prevent re-doing fruitless searches
        return antismash_hit_count