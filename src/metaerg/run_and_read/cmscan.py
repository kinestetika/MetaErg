import shutil
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor
from collections import namedtuple

from metaerg.run_and_read.context import register_annotator, spawn_file, run_external, DATABASE_DIR, CPUS_PER_GENOME
from metaerg.run_and_read.data_model import FeatureType, MetaergGenome


def _run_programs(genome:MetaergGenome, result_files):
    fasta_file = genome.write_fasta_files(spawn_file('masked', genome.id), masked=True)
    rfam_database = Path(DATABASE_DIR, "Rfam.cm")
    if CPUS_PER_GENOME > 1:
        split_fasta_files = genome.write_fasta_files(fasta_file, CPUS_PER_GENOME)
        split_cmscan_files = [Path(result_files[0].parent, f'{result_files[0].name}.{i}')
                              for i in range(len(split_fasta_files))]
        with ProcessPoolExecutor(max_workers=CPUS_PER_GENOME) as executor:
            for split_input, split_output in zip(split_fasta_files, split_cmscan_files):
                executor.submit(run_external, f'cmscan --rfam --tblout {split_output} {rfam_database} {split_input}')
        with open(result_files[0], 'wb') as output:
            for split_input_file, split_output_file in zip(split_fasta_files, split_cmscan_files):
                with open(split_output_file, 'rb') as input:
                    shutil.copyfileobj(input, output)
                split_input_file.unlink()
                split_output_file.unlink()
    else:
        run_external(f'cmscan --rfam --tblout {result_files[0]} {rfam_database} {fasta_file}')


def _read_results(genome:MetaergGenome, result_files) -> int:
    NON_CODING_RNA_TYPES = {'LSU_rRNA_bacteria': FeatureType.rRNA,
                            'LSU_rRNA_archaea': FeatureType.rRNA,
                            'LSU_rRNA_eukarya': FeatureType.rRNA,
                            'SSU_rRNA_bacteria': FeatureType.rRNA,
                            'SSU_rRNA_archaea': FeatureType.rRNA,
                            'SSU_rRNA_eukarya': FeatureType.rRNA,
                            'SSU_rRNA_microsporidia': FeatureType.rRNA,
                            '5S_rRNA': FeatureType.rRNA,
                            '5_8S_rRNA': FeatureType.rRNA,
                            'tmRNA': FeatureType.tmRNA,
                            'tRNA': FeatureType.tRNA}
    hits = []
    Hit = namedtuple('Hit', ('query_id', 'hit_id', 'query_start', 'query_end',
                             'query_strand', 'score', 'descr'))
    with open(result_files[0]) as hmm_handle:
        for line in hmm_handle:
            words = line.strip().split()
            words[17] = ' '.join(words[17:])
            match  words:
                case [*_] if line.startswith('#'):
                    continue
                case [hit, _, query, _, _, _, _, start, end, '-', _, _, _, _, score, _, '!', descr]:
                    hit = Hit(query, hit, int(end), int(start), -1, float(score), descr)
                case [hit, _, query, _, _, _, _, start, end, '+', _, _, _, _, score, _, '!', descr]:
                    hit = Hit(query, hit, int(start), int(end), 1, float(score), descr)
                case [*_]:
                    continue
            overlap = None
            for prev_hit in hits:
                if hit.query_id == prev_hit.query_id and hit.query_start < prev_hit.query_end and \
                        hit.query_end > prev_hit.query_start:
                    overlap = prev_hit  # overlap detected
                    break
            if overlap:
                if hit.score > overlap.score:
                    hits.remove(overlap)
                    hits.append(hit)
            else:
                hits.append(hit)
    for hit in hits:
        if hit.hit_id in NON_CODING_RNA_TYPES.keys():
            f_type = NON_CODING_RNA_TYPES[hit.hit_id]
        elif hit.hit_id.startswith('CRISPR'):
            f_type = FeatureType.crispr_repeat
        else:
            f_type = FeatureType.ncRNA
        contig = genome.contigs[hit.query_id]
        f = contig.spawn_feature(hit.query_start - 1, hit.query_end, hit.query_strand, f_type,
                                 inference='cmscan')
        f.description = "{} {}".format(hit.hit_id, hit.descr)
    return len(hits)


@register_annotator
def run_and_read_cmscan():
    return ({'pipeline_position': 21,
             'purpose': 'noncoding (RNA) gene prediction with cmscan',
             'programs': ('cmscan',),
             'result_files': ("cmscan",),
             'databses': ('Rfam.cm',),
             'run': _run_programs,
             'read': _read_results})
