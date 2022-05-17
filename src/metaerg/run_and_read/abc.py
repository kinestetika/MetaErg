import shutil
import os
from pathlib import Path

from metaerg import utils
from metaerg.run_and_read.data_model import MetaergGenome


class ExecutionEnvironment:
    def __init__(self, database_dir: Path, force=False, multi_mode=False, threads=1):
        self.database_dir = database_dir
        self.force = force
        self.multi_mode = multi_mode
        self.threads = threads

    def __repr__(self):
        return 'ExecutionEnvironment(database_dir={}, force={}, multi_mode={}, threads={})'.format(
            self.database_dir, self.force, self.multi_mode, self.threads)


class AbstractBaseClass:
    def __init__(self, genome: MetaergGenome, exec_env: ExecutionEnvironment):
        self.genome = genome
        self.exec = exec_env

    def _purpose(self) -> str:
        """Should return the purpose of the tool"""
        pass

    def _programs(self) -> tuple:
        """Should return a tuple with the programs needed"""
        pass

    def _databases(self) -> tuple:
        """Should return a tuple with database files needed"""
        return ()

    def _result_files(self) -> tuple:
        """Should return a tuple with the result files (Path objects) created by the programs"""
        pass

    def _run_programs(self):
        """Should execute the helper programs to complete the analysis"""
        pass

    def _read_results(self) -> int:
        """Should parse the result files and return the # of positives"""
        pass

    def __call__(self):
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

    def spawn_file(self, program_name):
        if self.exec.multi_mode:
            folder = Path(program_name)
            if not folder.exists():
                folder.mkdir()
            elif folder.is_file():
                if self.exec.force:
                    folder.unlink()
                    folder.mkdir()
                else:
                    raise Exception("Use force to overwrite existing results")
            return Path(folder, self.genome.id)
        else:
            file = Path(f'{self.genome.id}.{program_name}')
            if file.exists() and file.is_dir():
                if self.exec.force:
                    shutil.rmtree(file)
            return file
