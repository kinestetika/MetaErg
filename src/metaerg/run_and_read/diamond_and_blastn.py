from pathlib import Path

from metaerg.run_and_read.data_model import MetaergGenome, MetaergSeqFeature, BlastResult, DBentry, TabularBlastParser
from metaerg.run_and_read.context import register, spawn_file, run_external, DATABASE_DIR, CPUS_PER_GENOME, log


def _run_programs(genome:MetaergGenome, result_files):
    cds_aa_file = spawn_file('cds.faa', genome.id)
    rna_nt_file = spawn_file('rna.nt', genome.id)
    blastn_db = Path(DATABASE_DIR, "db_rna.fna")
    diamond_db = Path(DATABASE_DIR, "db_protein.faa")
    run_external(f'diamond blastp -d {diamond_db} -q {cds_aa_file} -o {result_files[0]} -f 6 '
                 f'--threads {CPUS_PER_GENOME} --max-target-seqs 10')
    run_external(f'blastn -db {blastn_db} -query {rna_nt_file} -out {result_files[1]} -max_target_seqs 10 -outfmt 6')


def _read_results(genome:MetaergGenome, result_files) -> int:
    # (1) load databse descriptions
    db_descr = Path(DATABASE_DIR, 'db_descriptions.txt')
    descriptions = {'p': {}, 'e': {}, 'v': {}}
    with open(db_descr) as descr_handle:
        for line in descr_handle:
            words = line.split('\t')
            descriptions[words[0]][int(words[1])] = words[2].strip()
    log(f'Parsed ({len(descriptions["p"])}, {len(descriptions["e"])}, {len(descriptions["v"])}) '
        f'gene descriptions from db for (prokaryotes, eukaryotes and viruses) respectively. ')
    # (2) load database taxonomy
    db_taxonomy = Path(DATABASE_DIR, 'db_taxonomy.txt')
    taxonomy = {'p': {}, 'e': {}, 'v': {}}
    with open(db_taxonomy) as taxon_handle:
        for line in taxon_handle:
            words = line.split('\t')
            taxonomy[words[0]][int(words[1])] = words[2].strip().replace('~', '; ')
    log(f'Parsed ({len(taxonomy["p"])}, {len(taxonomy["e"])}, {len(taxonomy["v"])}) '
        f'taxa from db for (prokaryotes, eukaryotes and viruses) respectively.')
    # (3) parse diamond blast results
    def get_db_entry(db_id):
        words = db_id.split('~')  # org_acc gene_acc [pev] gene# decr# taxon#
        return DBentry(db_id, '', descriptions[words[2]][int(words[4])], taxonomy[words[2]][int(words[5])],
                       int(words[6]), int(words[3]))

    def process_blast_result(blast_result: BlastResult):
        feature: MetaergSeqFeature = genome.get_feature(blast_result.query())
        feature.blast = blast_result
        feature.product = blast_result.summary()
        genome.subsystems.match(feature, (h.hit.descr for h in blast_result.hits if h.aligned_length / h.hit.length >= 0.8))

    blast_result_count = 0
    with TabularBlastParser(result_files[0], 'BLAST', get_db_entry) as handle:
        for blast_result in handle:
            blast_result_count += 1
            process_blast_result(blast_result)
    with TabularBlastParser(result_files[1], 'BLAST', get_db_entry) as handle:
        for blast_result in handle:
            blast_result_count += 1
            process_blast_result(blast_result)
    return blast_result_count


@register
def run_and_read_diamond_blastn():
    return ({'pipeline_position': 81,
             'purpose': 'function prediction and taxonomic classification of genes with diamond and blastn',
             'programs': ('diamond', 'blastn'),
             'databses': ('db_protein.faa', 'db_rna.fna', 'db_descriptions.txt', 'db_taxonomy.txt'),
             'result_files': ('diamond', 'blastn'),
             'run': _run_programs,
             'read': _read_results})
