import shutil
import os
from multiprocessing import cpu_count
from pathlib import Path
from metaerg import utils
from metaerg.run_and_read.data_model import MetaergGenome

annotator_registry = []

class ExecutionEnvironment:
    def __init__(self, **kwargs):
        self.database_dir = kwargs['database_dir']
        self.force = kwargs.get('force', False)
        self.rename_contigs = kwargs.get('rename_contigs', True)
        self.min_contig_length = int(kwargs.get('min_contig_length', 0))
        self.cpus_per_genome = int(kwargs.get('cpus', 0))
        cpus_available = cpu_count()
        if self.cpus_per_genome > 0:
            self.cpus_per_genome = min(self.cpus_per_genome, cpus_available)
            cpus_available = self.cpus_per_genome
        else:
            self.cpus_per_genome = cpus_available
        utils.log(f'Detected {cpus_available} available threads/cpus, will use {self.cpus_per_genome}.')
        self.contig_file = Path(kwargs['contig_file']).absolute()
        if self.contig_file.is_dir():
            self.contig_files = [x.absolute() for x in sorted(self.contig_file.glob(f'*{kwargs["file_extension"]}'))]
            if not len(self.contig_files):
                utils.log(f'Did not find any contig files with extension "{kwargs["file_extension"]}" '
                          f'in dir "{self.contig_file}"')
                exit(1)
            self.multi_mode = True
            self.cpus_per_genome = max(1, int(self.cpus_per_genome / len(self.contig_files)))
            self.parallel_annotations = min(cpus_available, len(self.contig_files))
            utils.log(f'Ready to annotate {len(self.contig_files)} genomes in dir "{self.contig_file}" with '
                      f'{self.cpus_per_genome} threads per genome.')
        else:
            self.multi_mode = False
            utils.log(f'Ready to annotate genome in nucleotide fasta file "{self.contig_file}".')

    def __repr__(self):
        return 'ExecutionEnvironment(database_dir={}, force={}, multi_mode={}, threads={})'.format(
            self.database_dir, self.force, self.multi_mode, self.cpus_per_genome)

    def spawn_file(self, program_name, genome_id) -> Path:
        if self.multi_mode:
            folder = Path(program_name)
            if not folder.exists():
                folder.mkdir()
            elif folder.is_file():
                if self.force:
                    folder.unlink()
                    folder.mkdir()
                else:
                    raise Exception("Use force to overwrite existing results")
            return Path(folder, genome_id)
        else:
            file = Path(f'{genome_id}.{program_name}')
            if file.exists() and file.is_dir():
                if self.force:
                    shutil.rmtree(file)
            return file


class Annotator:
    def __init__(self, genome: MetaergGenome, exec_env: ExecutionEnvironment):
        self.genome = genome
        self.exec = exec_env
        self.pipeline_position = 999

    def _purpose(self) -> str:
        """Returns the purpose of the annotaton tool."""
        pass

    def _programs(self) -> tuple:
        """Returns a tuple with the programs needed."""
        pass

    def _databases(self) -> tuple:
        """Returns a tuple with database files needed."""
        return ()

    def _result_files(self) -> tuple:
        """Returns a tuple with the result files (Path objects) created by the programs."""
        pass

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
                                         self._purpose()))
        # (1) First make sure that the helper programs are available:
        all_programs_in_path = True
        for p in self._programs():
            program_path = shutil.which(p, mode=os.X_OK)
            if not program_path:
                all_programs_in_path = False
                utils.log('({}) Unable to run {}, helper program "{}" not in path', (self.genome.id,
                                                                                     self._purpose(),
                                                                                     p))
        # (2) Then, make sure required databases are available
        for d in self._databases():
            if not d.exists() or not d.stat().st_size:
                utils.log('({}) Unable to run {}, or parse results, database "{}" missing', (self.genome.id,
                                                                                             self._purpose(),
                                                                                             d))
                return
        # (3) Then, if force or the results files are not yet there, run the programs:
        if all_programs_in_path:
            previous_results_missing = False
            for f in self._result_files():
                if not f.exists() or not f.stat().st_size:
                    previous_results_missing = True
                    break
            if self.exec.force or previous_results_missing:
                self._run_programs()
            else:
                utils.log('({}) Reusing existing results in {}.', (self.genome.id,
                                                                   self._result_files()))
        # (4) If all results files are there, read the results:
        all_results_created = True
        for f in self._result_files():
            if not f.exists() or not f.stat().st_size:
                all_results_created = False
                utils.log('({}) Missing expected result file {}.', (self.genome.id, f))
        if all_results_created:
            positive_count = self._read_results()
        else:
            positive_count = 0
        utils.log('({}) {} complete. Found {}.', (self.genome.id,
                                                  self._purpose(),
                                                  positive_count))
        # (5) Save (intermediate) results:
        gbk_file = self.exec.spawn_file("gbk", self.genome.id)
        gff_file = self.exec.spawn_file("gff", self.genome.id)
        self.genome.write_gbk_gff(gbk_file=gbk_file, gff_file=gff_file)


def register(annotator: Annotator):
    annotator_registry.append(annotator)
    annotator_registry.sort(key=lambda annotator: annotator.pipeline_position)
    return annotator
