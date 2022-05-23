import shutil
import os
from multiprocessing import cpu_count
from pathlib import Path
from metaerg import utils
from metaerg.run_and_read.data_model import MetaergGenome

annotator_registry = []

class ExecutionEnvironment:
    def __init__(self, contig_file, database_dir, rename_contigs, rename_genomes, min_contig_length, cpus, force,
                 file_extension, translation_table, checkm_dir, gtdbtk_dir):
        utils.log('Initializing execution environment with commenad line arguments...')
        self.contig_file = Path(contig_file).absolute()
        self.base_dir = self.contig_file if self.contig_file.is_dir() else self.contig_file.parent
        self.temp_dir = Path(self.base_dir, "temp")
        self.database_dir = Path(database_dir).absolute()
        self.checkm_dir = Path(checkm_dir).absolute()
        self.gtdbtk_dir = Path(gtdbtk_dir).absolute()
        self.genome_name_mapping_file = Path(self.temp_dir, 'genome.name.mapping.txt')
        self.html_dir = Path(self.base_dir, "html")

        self.multi_mode = self.contig_file.is_dir()
        self.rename_contigs = rename_contigs,
        self.rename_genomes = rename_genomes
        self.min_contig_length = min_contig_length
        self.force = force
        self.file_extension = file_extension
        self.translation_table = translation_table
        self.contig_files = None
        self.genome_names = None

        self.cpus_per_genome = cpus
        self.cpus_available = cpu_count()

    def __repr__(self):
        return '{}({})'.format(self.__class__.__name__, ', '.join(f'{k}={v}' for k, v in self.__dict__))

    def prep_environment(self):
        if not self.contig_file.exists():
            utils.log(f'Input file "{self.contig_file}" is missing. Expecting dir or a nt fasta file.')
            exit(1)
        if self.temp_dir.exists():
            utils.log('Warning: may overwrite existing temp files...')
            if self.temp_dir.is_file():
                utils.log(f'Expected folder at {self.temp_dir}, found regular file, crash! Delete this file first')
                exit(1)
        else:
            os.mkdir(self.temp_dir)
        shutil.rmtree(self.html_dir, ignore_errors=True)
        os.mkdir(self.html_dir)
        if self.cpus_per_genome > 0:
            self.cpus_per_genome = min(self.cpus_per_genome, self.cpus_available)
            self.cpus_available = self.cpus_per_genome
        else:
            self.cpus_per_genome = self.cpus_available
        utils.log(f'Detected {self.cpus_available} available threads/cpus, will use {self.cpus_per_genome}.')

        if self.contig_file.is_dir():
            self.contig_files = [x.absolute() for x in sorted(self.contig_file.glob(f'*{self.file_extension}'))]
            if not len(self.contig_files):
                utils.log(f'Did not find any contig files with extension "{self.file_extension}" '
                          f'in dir "{self.contig_file}"')
                exit(1)
            self.cpus_per_genome = max(1, int(self.cpus_per_genome / len(self.contig_files)))
            self.parallel_annotations = min(self.cpus_available, len(self.contig_files))
        else:
            self.contig_files = [self.contig_file]
        if self.rename_genomes:
            self.genome_names = [f'g{self.contig_files.index(f):0>4}' for f in self.contig_files]
        else:
            self.genome_names = [f.name for f in self.contig_files]
        utils.log(f'writing genome names to {self.genome_name_mapping_file} ')
        with open(self.genome_name_mapping_file, 'w') as mapping_file:
            for n, o in zip(self.genome_names, self.contig_files):
                mapping_file.write(f'{n}\t{o.stem}\t{o}\n')
        utils.log(f'Ready to annotate {len(self.contig_files)} genomes in dir "{self.base_dir}" with '
                  f'{self.cpus_per_genome} threads per genome.')

    def spawn_file(self, program_name, genome_id, base_dir = None) -> Path:
        """computes a Path genome_id.program_name or, if multimode==True, program_name/genome_id"""
        target_dir = base_dir if base_dir else self.temp_dir
        if self.multi_mode:
            dir = Path(target_dir, program_name)
            if not dir.exists():
                dir.mkdir()
            elif dir.is_file():
                if self.force:
                    dir.unlink()
                    dir.mkdir()
                else:
                    raise Exception("Use force to overwrite existing results")
            return Path(dir, genome_id)
        else:
            file = Path(target_dir, f'{genome_id}.{program_name}')
            if file.exists() and file.is_dir():
                if self.force:
                    shutil.rmtree(file)
            return file


class Annotator:
    def __init__(self, genome: MetaergGenome, exec_env: ExecutionEnvironment):
        self.genome = genome
        self.exec = exec_env
        self._pipeline_position = 999
        self._purpose = 'override'
        self._programs = tuple()
        self._databases = tuple()
        self._result_files = tuple()

    def __repr__(self):
        return '{}({})'.format(self.__class__.__name__, ', '.join(f'{k}={v}' for k, v in self.__dict__))

    def _run_programs(self):
        """Executes the helper programs to complete the analysis."""
        pass

    def _read_results(self) -> int:
        """Parses the result files and returns the # of positives."""
        pass

    def __call__(self):
        self.run_and_read()

    def spawn_file(self, program_name) -> Path:
        """Comvenience wrapper for exec.spawn_file."""
        return self.exec.spawn_file(program_name, self.genome.id)

    def run_and_read(self):
        """Runs programs and reads results."""
        utils.log('({}) {} started...', (self.genome.id,
                                         self._purpose))
        # (1) First make sure that the helper programs are available:
        all_programs_in_path = True
        for p in self._programs:
            program_path = shutil.which(p, mode=os.X_OK)
            if not program_path:
                all_programs_in_path = False
                utils.log('({}) Unable to run {}, helper program "{}" not in path', (self.genome.id,
                                                                                     self._purpose,
                                                                                     p))
        # (2) Then, make sure required databases are available
        for d in self._databases:
            if not d.exists() or not d.stat().st_size:
                utils.log('({}) Unable to run {}, or parse results, database "{}" missing', (self.genome.id,
                                                                                             self._purpose,
                                                                                             d))
                return
        # (3) Then, if force or the results files are not yet there, run the programs:
        if all_programs_in_path:
            previous_results_missing = False
            for f in self._result_files:
                if not f.exists() or not f.stat().st_size:
                    previous_results_missing = True
                    break
            if self.exec.force or previous_results_missing:
                self._run_programs()
            else:
                utils.log('({}) Reusing existing results in {}.', (self.genome.id,
                                                                   self._result_files))
        # (4) If all results files are there, read the results:
        all_results_created = True
        for f in self._result_files:
            if not f.exists() or not f.stat().st_size:
                all_results_created = False
                utils.log('({}) Missing expected result file {}.', (self.genome.id, f))
        if all_results_created:
            positive_count = self._read_results()
        else:
            positive_count = 0
        utils.log('({}) {} complete. Found {}.', (self.genome.id,
                                                  self._purpose,
                                                  positive_count))
        # (5) Save (intermediate) results:
        gbk_file = self.exec.spawn_file("gbk", self.genome.id)
        gff_file = self.exec.spawn_file("gff", self.genome.id)
        self.genome.write_gbk_gff(gbk_file=gbk_file, gff_file=gff_file)


def register(annotator: Annotator):
    annotator_registry.append(annotator)
    annotator_registry.sort(key=lambda annotator: annotator.pipeline_position)
    return annotator
