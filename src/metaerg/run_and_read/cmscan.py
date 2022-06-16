import shutil
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor
from collections import namedtuple

from metaerg import context
from metaerg import bioparsers
from metaerg.data_model import FeatureType, SeqFeature, Genome


def _run_programs(genome: Genome, result_files):
    fasta_file = context.spawn_file('masked', genome.id)
    rfam_database = Path(context.DATABASE_DIR, "Rfam.cm")
    if context.CPUS_PER_GENOME > 1:
        split_fasta_files = bioparsers.write_genome_to_fasta_files(genome, fasta_file, context.CPUS_PER_GENOME, mask=True)
        split_cmscan_files = [Path(result_files[0].parent, f'{result_files[0].name}.{i}')
                              for i in range(len(split_fasta_files))]
        with ProcessPoolExecutor(max_workers=context.CPUS_PER_GENOME) as executor:
            for split_input, split_output in zip(split_fasta_files, split_cmscan_files):
                executor.submit(context.run_external, f'cmscan --rfam --tblout {split_output} {rfam_database} '
                                                      f'{split_input}')
        with open(result_files[0], 'wb') as output:
            for split_input_file, split_output_file in zip(split_fasta_files, split_cmscan_files):
                with open(split_output_file, 'rb') as input:
                    shutil.copyfileobj(input, output)
                split_input_file.unlink()
                split_output_file.unlink()
    else:
        bioparsers.write_genome_to_fasta_files(genome, fasta_file, mask=True)
        context.run_external(f'cmscan --rfam --tblout {result_files[0]} {rfam_database} {fasta_file}')


def _read_results(genome:Genome, result_files) -> int:
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
        seq = contig.seq[hit.query_start - 1:hit.query_end]
        if hit.query_strand < 0:
            seq = bioparsers.reverse_complement(seq)
        feature = SeqFeature(hit.query_start - 1, hit.query_end, hit.query_strand, f_type, seq=seq,
                             inference='cmscan', descr = "{} {}".format(hit.hit_id, hit.descr))
        contig.features.append(feature)
    return len(hits)


@context.register_annotator
def run_and_read_cmscan():
    return ({'pipeline_position': 21,
             'purpose': 'noncoding (RNA) gene prediction with cmscan',
             'programs': ('cmscan',),
             'result_files': ("cmscan",),
             'databases': (Path(context.DATABASE_DIR, 'rfam', 'Rfam.cm'),),
             'run': _run_programs,
             'read': _read_results})


@context.register_database_installer
def install_cmscan_database():
    if 'R' not in context.CREATE_DB_TASKS:
        return
    rfam_dir = Path(context.DATABASE_DIR, 'rfam')
    rfam_dir.mkdir(exist_ok=True, parents=True)

    rfam_file = Path(rfam_dir, 'Rfam.cm')
    if context.FORCE or not rfam_file.exists() or not rfam_file.stat().st_size:
        context.log(f'Installing the RFAM database to {rfam_file}...')
        context.run_external(
            f'wget -P {rfam_dir} http://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/Rfam.cm.gz')
        context.run_external(f'gunzip {rfam_file}.gz')
    else:
        context.log(f'Keeping existing conserved domain database in {rfam_file}, use --force to overwrite.')
    if context.FORCE or not Path(context.DATABASE_DIR, "Rfam.cm.i1f").exists():
        context.log(f'Running cmpress...')
        context.run_external(f'cmpress -F {rfam_file}')
    else:
        context.log('Skipping cmpress for previously cmpressed RFAM database...')
