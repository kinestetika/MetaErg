import shutil
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor
from collections import namedtuple

from metaerg.run_and_read.abc import Annotator, register, ExecutionEnvironment
from metaerg.run_and_read.data_model import FeatureType
from metaerg import utils

NON_CODING_RNA_TYPES = {'LSU_rRNA_bacteria': FeatureType.rRNA,
                        'LSU_rRNA_archaea':  FeatureType.rRNA,
                        'LSU_rRNA_eukarya':  FeatureType.rRNA,
                        'SSU_rRNA_bacteria': FeatureType.rRNA,
                        'SSU_rRNA_archaea':  FeatureType.rRNA,
                        'SSU_rRNA_eukarya':  FeatureType.rRNA,
                        'SSU_rRNA_microsporidia': FeatureType.rRNA,
                        '5S_rRNA': FeatureType.rRNA,
                        '5_8S_rRNA': FeatureType.rRNA,
                        'tmRNA': FeatureType.tmRNA,
                        'tRNA': FeatureType.tRNA}


@register
class CMScan(Annotator):
    def __init__(self, genome, exec_env: ExecutionEnvironment):
        super().__init__(genome, exec_env)
        self.cmscan_file = self.spawn_file('cmscan')
        self._pipeline_position = 21
        self._purpose = 'noncoding (RNA) gene prediction with cmscan'
        self._programs = ('cmscan',)
        self._result_files = (self.cmscan_file,)

    def _run_programs(self):
        """Executes the helper programs to complete the analysis"""
        fasta_file = self.genome.write_fasta_files(self.spawn_file('masked'), masked=True)
        rfam_database = Path(self.exec.database_dir, "Rfam.cm")
        if self.exec.cpus_per_genome > 1:
            split_fasta_files = self.genome.write_fasta_files(fasta_file, self.exec.cpus_per_genome)
            split_cmscan_files = [Path(self.cmscan_file.parent, f'{self.cmscan_file.name}.{i}')
                                  for i in range(len(split_fasta_files))]
            with ProcessPoolExecutor(max_workers=self.exec.cpus_per_genome) as executor:
                for split_input, split_output in zip(split_fasta_files, split_cmscan_files):
                    executor.submit(utils.run_external, f'cmscan --rfam --tblout {split_output} '
                                                        f'{rfam_database} {split_input}')
            with open(self.cmscan_file, 'wb') as output:
                for split_input_file, split_output_file in zip(split_fasta_files, split_cmscan_files):
                    with open(split_output_file, 'rb') as input:
                        shutil.copyfileobj(input, output)
                    split_input_file.unlink()
                    split_output_file.unlink()
        else:
            utils.run_external(f'cmscan --rfam --tblout {self.cmscan_file} {rfam_database} {fasta_file}')

    def _read_results(self) -> int:
        """Should parse the result files and return the # of positives"""
        hits = []
        Hit = namedtuple('Hit', ('query_id', 'hit_id', 'query_start', 'query_end',
                                 'query_strand', 'score', 'descr'))
        with open(self.cmscan_file) as hmm_handle:
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
            contig = self.genome.contigs[hit.query_id]
            f = contig.spawn_feature(hit.query_start - 1, hit.query_end, hit.query_strand, f_type,
                                     inference='cmscan')
            f.description = "{} {}".format(hit.hit_id, hit.descr)
        return len(hits)
