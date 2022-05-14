import shutil
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor
from collections import namedtuple

from Bio.SeqFeature import FeatureLocation
from metaerg.run_and_read import abc
from metaerg import utils

NON_CODING_RNA_TYPES = {'LSU_rRNA_bacteria': 'rRNA',
                        'LSU_rRNA_archaea':  'rRNA',
                        'LSU_rRNA_eukarya':  'rRNA',
                        'SSU_rRNA_bacteria': 'rRNA',
                        'SSU_rRNA_archaea':  'rRNA',
                        'SSU_rRNA_eukarya':  'rRNA',
                        'SSU_rRNA_microsporidia': 'rRNA',
                        '5S_rRNA': 'rRNA',
                        '5_8S_rRNA': 'rRNA',
                        'tmRNA': 'tmRNA',
                        'tRNA': 'tRNA'}


class CMScan(abc.AbstractBaseClass):
    def __init__(self, genome, exec_env: abc.ExecutionEnvironment):
        super().__init__(genome, exec_env)
        self.cmscan_file = self.spawn_file('cmscan')

    def __repr__(self):
        return f'CMScan({self.genome}, {self.exec})'

    def __purpose__(self) -> str:
        """Should return the purpose of the tool"""
        return 'noncoding (RNA) gene prediction with cmscan'

    def __programs__(self) -> tuple:
        """Should return a tuple with the programs needed"""
        return 'cmscan',

    def __result_files__(self) -> tuple:
        """Should return a tuple with the result files (Path objects) created by the programs"""
        return self.cmscan_file,

    def __run_programs__(self):
        """Should execute the helper programs to complete the analysis"""
        fasta_file = self.genome.make_masked_contig_fasta_file(self.spawn_file('masked'))
        rfam_database = Path(self.exec.database_dir, "Rfam.cm")
        if self.exec.threads > 1:
            split_fasta_files = self.genome.make_split_fasta_files(fasta_file, self.exec.threads, target='contig')
            split_cmscan_files = [Path(self.cmscan_file.parent, f'{self.cmscan_file.name}.{i}')
                                  for i in range(len(split_fasta_files))]
            with ProcessPoolExecutor(max_workers=self.exec.threads) as executor:
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

    def __read_results__(self) -> int:
        """Should parse the result files and return the # of positives"""
        hits = []
        with open(self.cmscan_file) as hmm_handle:
            for line in hmm_handle:
                if line.startswith('#'):
                    continue
                words = line.split()
                if len(words) < 18 or '?' == words[16]:
                    continue
                words[17] = ' '.join(words[17:])
                Hit = namedtuple('Hit', ('hit_id', 'query_id', 'query_start', 'query_end',
                                         'query_strand', 'score', 'descr'))
                if '-' == words[9]:  # this is the strand, if it is -, reverse start and end
                    hit = Hit(hit_id=words[0], query_id=words[2], query_start=int(words[8]), query_end=int(words[7]),
                              query_strand=-1, score=float(words[14]), descr=words[17])
                else:
                    hit = Hit(hit_id=words[0], query_id=words[2], query_start=int(words[7]), query_end=int(words[8]),
                              query_strand=1, score=float(words[14]), descr=words[17])
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
                f_type = 'crispr'
            else:
                f_type = 'ncRNA'
            contig = self.genome.contigs[hit.query_id]
            location = FeatureLocation(hit.query_start - 1, hit.query_end, strand=hit.query_strand)
            f = contig.spawn_feature(f_type, location, 'cmscan')
            f.description = "{} {}".format(hit.hit_id, hit.descr)
        return len(hits)
