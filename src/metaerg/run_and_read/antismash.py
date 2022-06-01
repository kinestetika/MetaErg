import shutil
from pathlib import Path

from Bio import SeqIO

from metaerg import context
from metaerg.data_model import MetaergSeqFeature, MetaergGenome


def _run_programs(genome:MetaergGenome, result_files):
    """Should execute the helper programs to complete the analysis"""
    gbk_file = context.spawn_file('gbk', genome.id)
    if result_files[0].exists():
        shutil.rmtree(result_files[0])
    context.run_external(f'antismash --genefinding-tool none --output-dir {result_files[0]} {gbk_file}')


def _read_results(genome:MetaergGenome, result_files) -> int:
    """Should parse the result files and return the # of positives."""
    antismash_hit_count = 0
    for f in sorted(result_files[0].glob('*region*.gbk')):
        with open(f) as handle:
            antismash_region_name = ''
            antismash_region_number = 0
            for antismash_record in SeqIO.parse(handle, 'genbank'):
                for antismash_feature in antismash_record.features:
                    antismash_hit_count += 1
                    if 'region' == antismash_feature.type:
                        antismash_region_name = antismash_feature.qualifiers['rules'][0]
                        antismash_region_number = int(antismash_feature.qualifiers['region_number'][0])
                    elif 'CDS' in antismash_feature.type:
                        feature: MetaergSeqFeature = genome.get_feature(antismash_feature.qualifiers['locus_tag'][0])
                        antismash_gene_function = antismash_feature.qualifiers['gene_functions'][0]
                        antismash_gene_category = antismash_feature.qualifiers['gene_kind']
                        if antismash_region_name:
                            feature.antismash = ' '.join((f'(region {antismash_region_number})',
                                                         antismash_region_name,
                                                         antismash_gene_function,
                                                         antismash_gene_category))
                            feature.subsystem.add('[secondary-metabolites]')
                            genome.subsystems.subsystems['[secondary-metabolites]'].add_hit(feature.id)
    if not antismash_hit_count:
        result_files[0].mkdir(exist_ok=True)  # to prevent re-doing fruitless searches
    return antismash_hit_count


def write_html(self, filename=None):
    """need to copy the antismash result dir to the metaerg html dir."""
    shutil.copytree(self.antismash_file, Path(self.exec.html_dir, self.genome.id, 'antismash'))


@context.register_annotator
def run_and_read_antismash():
    return ({'pipeline_position': 91,
             'purpose': 'prediction of secondary metabolite genes with antismash',
             'programs': ('antismash',),
             'result_files': ('antismash',),
             'run': _run_programs,
             'read': _read_results})
