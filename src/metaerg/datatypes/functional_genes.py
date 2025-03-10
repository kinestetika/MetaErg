from pathlib import Path
from math import log10

from metaerg import context
from metaerg.datatypes import sqlite
from metaerg.datatypes.blast import BlastResult, BlastHit

GENES = []
SUBSYSTEMS = []
CONFIG_DATA_FILE_EXTENSION = '.config.txt'

class FunctionalGene:
    def __init__(self, subsystem, gene, confidence=1.0):
        self.subsystem = subsystem
        self.gene = gene
        self.confidence = confidence
        self.pathway_positions = []
        self.cues = frozenset()
        self.anti_cues = frozenset()

    def __iter__(self):
        for k, v in self.__dict__.items():
            if k == 'cues' or k == 'pathway_positions' or k == 'anti_cues':
                continue
            yield k, v

    def __repr__(self):
        return '{}({})'.format(type(self).__name__, ', '.join(f'{k}={v!r}' for (k, v) in self))

    def __len__(self):
        return 3

    def __str__(self):
        if 0 < self.confidence < 1.0:
            return f'[{self.subsystem}|{self.gene}|{self.confidence:.1f}]'
        else:
            return f'[{self.subsystem}|{self.gene}]'

    def __key(self):
        return (self.subsystem, self.gene)

    def __hash__(self):
        return hash(self.__key())

    def __eq__(self, other):
        if isinstance(other, FunctionalGene):
            return self.__key() == other.__key()
        return NotImplemented


SECONDARY_METABOLITE_GENE = FunctionalGene('Secondary-metabolites', 'Secondary-metabolites', 1)


def format_list_of_subsystem_genes(subsystem_genes: list):
    return ' '.join(str(sg) for sg in subsystem_genes)


def parse_functional_gene_config_file(file: Path):
    current_subsystem = None
    with open(file) as handle:
        for line in handle:
            line = line.strip()
            if line.startswith("#") or not len(line):
                continue
            elif line.startswith(">"):
                current_subsystem = line[1:]
                SUBSYSTEMS.append(current_subsystem)
            elif current_subsystem is not None:
                pathway_positions = []
                si = line.find(' ')
                try:
                    pathway_positions = [int(s) for s in line[:si].split('|')]
                    line = line[si + 1:].strip()
                    si = line.find(' ')
                except ValueError:
                    pass
                new_gene_def = FunctionalGene(current_subsystem, line[si + 1:])
                try:
                    (cue_text, anti_cue_text) = line[:si].split('!')
                    new_gene_def.anti_cues = frozenset(anti_cue_text.split('|'))
                except ValueError:
                    cue_text = line[:si]
                new_gene_def.cues = frozenset(cue_text.split('|'))
                new_gene_def.pathway_positions = pathway_positions
                GENES.append(new_gene_def)
    context.log(f'Parsed {len(SUBSYSTEMS)} subsystems, {len(GENES)} genes from {file} ...')
    return len(GENES)


def init_functional_gene_config():
    context.log('Parsing functional gene config files...')
    base_config_file = Path(__file__).parent.parent / 'run_and_read' / 'data' / 'functional_gene_data'
    gene_count = parse_functional_gene_config_file(base_config_file)
    user_config_dir = context.DATABASE_DIR / 'hmm' / 'user_config'
    if user_config_dir.exists() and user_config_dir.is_dir():
        for config_file in (context.DATABASE_DIR / 'user_config').glob('*' + CONFIG_DATA_FILE_EXTENSION):
            gene_count += parse_functional_gene_config_file(user_config_dir / config_file.name)
    SUBSYSTEMS.append('Secondary-metabolites')
    SUBSYSTEMS.sort()
    GENES.append(SECONDARY_METABOLITE_GENE)
    context.log(f'Complete configuration contains {gene_count} genes.')


def match(blast_result: BlastResult, number_of_hits_considered) -> list:
    matches = []
    i = 0
    for h in blast_result.hits:
        i += 1
        if i > number_of_hits_considered:
            break
        if new_matches := match_hit(h, blast_result):
            for new_match in new_matches:  # Add only one hit to any distinct functional gene
                if not new_match in matches:
                    matches.append(new_match)
    return matches


def match_hit(blast_hit: BlastHit, blast_result: BlastResult) -> list:
    matches = []
    for gene_def in GENES:
        if blast_hit.hit.min_score:
            if blast_hit.score < blast_hit.hit.min_score:
                continue
            confidence = blast_hit.score / blast_hit.hit.min_t_score
        elif blast_hit.evalue > 0:
            confidence = min(1.0, - log10(blast_hit.evalue) / 100)
        else:
            confidence = 1.0
        if blast_hit.hit.accession in gene_def.cues and blast_hit.aligned_length >= 0.7 * blast_hit.hit.length:
            anti_ques_ok = True
            for other_blast_hit in blast_result.hits:
                if other_blast_hit.hit.accession in gene_def.anti_cues and blast_hit.evalue > other_blast_hit.evalue / 1e10:
                    # (1e10 based on inspection of evalues of in-group sequences to nir outgroup hmms... the difference was actually > ~1e50)
                    anti_ques_ok = False
                    break
            if anti_ques_ok:
                new_match = FunctionalGene(gene_def.subsystem, gene_def.gene, confidence)
                if not new_match in matches:
                    matches.append(new_match)
    return matches


def aggregate(db_connection):
    aggregated_subsystem_data = {k: {} for k in SUBSYSTEMS}
    for gene_def in GENES:
        aggregated_subsystem_data[gene_def.subsystem][gene_def.gene] = []
    for feature in sqlite.read_all_features(db_connection):
        for subsystem in feature.subsystems:
            if subsystem.confidence > 0.25:
                if subsystem.subsystem not in aggregated_subsystem_data:
                    aggregated_subsystem_data[subsystem.subsystem] = {}
                if subsystem.gene not in aggregated_subsystem_data[subsystem.subsystem]:
                    aggregated_subsystem_data[subsystem.subsystem][subsystem.gene] = []
                aggregated_subsystem_data[subsystem.subsystem][subsystem.gene].append(feature.id)
    return aggregated_subsystem_data


def get_subsystem_completeness(subsystem_name: str, genes: dict):
    denominator = 0
    positions_found = set()
    relevant_feature_ids = set()  # this is for "unsorted" systems - we only want to count each feature once
    for gene_def in GENES:
        if gene_def.subsystem == subsystem_name:
            for d in gene_def.pathway_positions:
                denominator = max(denominator, d)
    for gene_def in GENES:
        if gene_def.subsystem == subsystem_name:
            if gene_def.gene in genes.keys() and len(genes[gene_def.gene]):
                if denominator:
                    for p in gene_def.pathway_positions:
                        positions_found.add(p)
                else:
                    for feature_id in genes[gene_def.gene]:
                        relevant_feature_ids.add(feature_id)
    return len(positions_found) / denominator if denominator > 0 else len(relevant_feature_ids)


def test_alignment_cdd_and_functional_gene_data(database_dir):
    print('initing config')
    init_functional_gene_config()

    cdd_index = {}
    cdd_index_file = Path(database_dir, 'cdd', 'cddid.tbl')
    with open(cdd_index_file) as db_handle:
        for line in db_handle:
            words = line.strip().split("\t")
            cdd_index[words[1]] = f'{words[2]} {words[3]}'
    context.log(f'Parsed {len(cdd_index)} entries from conserved domain database.')

    for gene in GENES:
        for cue in gene.cues:
            if not cue in cdd_index.keys() and (cue.startswith('cl') or cue.startswith('PRK')):
                print(f'{cue} {gene.subsystem} {gene.gene[:15]}')
        for cue in gene.anti_cues:
            if not cue in cdd_index.keys() and (cue.startswith('cl') or cue.startswith('PRK')):
                print(f'{cue}')


def main ():
    test_alignment_cdd_and_functional_gene_data('/bio/data/databases/metaerg')


if __name__ == "__main__":
    main()
