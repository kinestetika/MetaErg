import shutil
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor
from collections import namedtuple

from metaerg import context
from metaerg.datatypes import fasta
from metaerg.datatypes import sqlite


def _run_programs(genome_name, contig_dict, db_connection, result_files):
    fasta_file = context.spawn_file('masked', genome_name)
    rfam_database = Path(context.DATABASE_DIR, 'rfam', 'Rfam.cm')
    if context.CPUS_PER_GENOME > 1:
        split_fasta_files = fasta.write_contigs_to_fasta(contig_dict, fasta_file, db_connection, genome_name,
                                                         mask_targets=fasta.ALL_MASK_TARGETS,
                                                         split=context.CPUS_PER_GENOME)
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
        fasta.write_contigs_to_fasta(contig_dict, fasta_file, db_connection, genome_name,
                                     mask_targets=fasta.ALL_MASK_TARGETS)
        context.run_external(f'cmscan --rfam --tblout {result_files[0]} {rfam_database} {fasta_file}')


def _read_results(genome_name, contig_dict, db_connection, result_files) -> int:
    NON_CODING_RNA_TYPES = {'LSU_rRNA_bacteria': 'rRNA',
                            'LSU_rRNA_archaea': 'rRNA',
                            'LSU_rRNA_eukarya': 'rRNA',
                            'SSU_rRNA_bacteria': 'rRNA',
                            'SSU_rRNA_archaea': 'rRNA',
                            'SSU_rRNA_eukarya': 'rRNA',
                            'SSU_rRNA_microsporidia': 'rRNA',
                            '5S_rRNA': 'rRNA',
                            '5_8S_rRNA': 'rRNA',
                            'tmRNA': 'tmRNA',
                            'tRNA': 'tRNA'}
    hits = []
    Hit = namedtuple('Hit', ('query_id', 'hit_id', 'query_start', 'query_end',
                             'query_strand', 'score', 'descr'))
    with open(result_files[0]) as hmm_handle:
        for line in hmm_handle:
            words = line.strip().split()
            if line.startswith('#') or len(words) < 17:
                continue
            words[17] = ' '.join(words[17:])
            match words:
                case [hit, _, query, _, _, _, _, start, end, '-', _, _, _, _, score, _, '!', descr, *_]:
                    hit = Hit(query, hit, int(end), int(start), -1, float(score), descr)
                case [hit, _, query, _, _, _, _, start, end, '+', _, _, _, _, score, _, '!', descr, *_]:
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
    count = 0
    for hit in hits:
        if hit.hit_id in NON_CODING_RNA_TYPES.keys():
            f_type = NON_CODING_RNA_TYPES[hit.hit_id]
        elif hit.hit_id.startswith('CRISPR'):
            f_type = 'crispr_repeat'
        else:
            f_type = 'ncRNA'
        contig = contig_dict[hit.query_id]
        seq = contig['seq'][hit.query_start - 1:hit.query_end]
        if hit.query_strand < 0:
            seq = fasta.reverse_complement(seq)
        feature = sqlite.Feature(genome = genome_name,
                   contig = hit.query_id,
                   start = hit.query_start - 1,
                   end = hit.query_end,
                   strand = hit.query_strand,
                   type = f_type,
                   inference = 'cmscan',
                   nt_seq = seq,
                   descr = "{} {}".format(hit.hit_id, hit.descr))
        sqlite.add_new_feature_to_db(db_connection, feature)
        count += 1
    return count


@context.register_annotator
def run_and_read_cmscan():
    return ({'pipeline_position': 21,
             'annotator_key': 'cmscan',
             'purpose': 'noncoding (RNA) gene prediction with cmscan',
             'programs': ('cmscan',),
             'result_files': ("cmscan",),
             'databases': (Path('rfam', 'Rfam.cm'),),
             'run': _run_programs,
             'read': _read_results})


@context.register_database_installer
def install_cmscan_database():
    if 'R' not in context.TASKS:
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
