import shutil
import os
from pathlib import Path

from metaerg import utils
from metaerg.run_and_read.data import MetaergSeqRecord
from metaerg.run_and_read.data import MetaergGenome

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import FeatureLocation


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
    def __init__(self, genome:MetaergGenome, exec:ExecutionEnvironment):
        self.genome = genome
        self.exec = exec

    def __purpose__(self) -> str:
        """Should return the purpose of the tool"""
        pass

    def __programs__(self) -> tuple:
        """Should return a tuple with the programs needed"""
        pass

    def __result_files__(self) -> tuple:
        """Should return a tuple with the result files (Path objects) created by the programs"""
        pass

    def __run_programs__(self):
        """Should execute the helper programs to complete the analysis"""
        pass

    def __read_results__(self) -> int:
        """Should parse the result files and return the # of positives"""
        pass

    def __call__(self):
        utils.log('({}) {} started...', (self.genome.name,
                                         self.__purpose__()))
        # (1) First make sure that the helper programs are available:
        all_programs_in_path = True
        for p in self.__programs__():
            program_path = shutil.which(p, mode=os.X_OK)
            if not program_path:
                all_programs_in_path = False
                utils.log('({}) Unable to run {}, helper program "{}" not in path', (self.genome.name,
                                                                                     self.__purpose__(),
                                                                                     p))
        # (2) Then, if force or the results files are not yet there, run the programs:
        if all_programs_in_path:
            previous_results_missing = False
            for f in self.__result_files__():
                if not f.exists() or not f.stat().st_size:
                    previous_results_missing = True
                    break
            if self.exec.force or previous_results_missing:
                self.__run_programs__()
            else:
                utils.log('({}) Reusing existing results in {}.', (self.genome.name,
                                                                   self.__result_files__()))
        # (3) If all results files are there, read the results:
        all_results_created = True
        for f in self.__result_files__():
            if not f.exists() or not f.stat().st_size:
                all_results_created = False
                utils.log('({}) Missing expected result file {}.', (self.genome.name, f))
        if all_results_created:
            positive_count = self.__read_results__()
        else:
            positive_count = 0
        utils.log('({}) {} complete. Found {}.', (self.genome.name,
                                                  self.__purpose__(),
                                                  positive_count))

    def spawn_file(self, program_name):
        if self.exec.multi_mode:
            dir =  Path(program_name)
            if not dir.exists():
                dir.mkdir()
            elif dir.is_file():
                if self.exec.force:
                    dir.unlink()
                    dir.mkdir()
                else:
                    raise Exception("Use force to overwrite existing results")
            return Path(dir, self.genome.name)
        else:
            file = Path(f'{self.genome.name}.{program_name}')
            if file.exists() and file.is_dir():
                if self.exec.force:
                    shutil.rmtree(file)
            return file

    def _mask_seq(self, record:MetaergSeqRecord, exceptions=None, min_mask_length=50) -> SeqRecord:
        seq = record.sequence
        self.nt_total += len(seq)
        for f in record.features:
            if f.inference in exceptions:
                continue
            if len(f.location) < min_mask_length:
                continue
            fl:FeatureLocation = f.location
            self.nt_masked += fl.end - fl.start
            seq[fl.start:fl.end] = 'N' * (fl.end - fl.start)
            #seq = seq[:fl.start] + 'N' * (fl.end - fl.start) + seq[fl.end:]
        return SeqRecord(Seq(seq), id=record.id, description=record.description)

    def make_masked_contig_fasta_file(self, exceptions=None, min_mask_length=50) -> Path:
        if exceptions is None:
            exceptions = set()
        (self.nt_masked, self.nt_total) = (0, 0)
        masked_fasta_file = self.spawn_file('masked')
        seq_iterator = (self._mask_seq(record, exceptions=exceptions, min_mask_length=min_mask_length)
                        for record in self.genome.contigs.values())
        SeqIO.write(seq_iterator, masked_fasta_file, "fasta")
        utils.log(f'Masked {self.nt_masked / self.nt_total * 100:.1f}% of sequence data.')
        return masked_fasta_file

    def make_split_fasta_files(self, base_file: Path, target):
            if 'CDS' == target:
                count = sum(1 for contig in self.genome.contigs.values() for f in contig.features if f.type == 'CDS')
            else:
                count = len(self.genome.contigs)
            number_of_files = min(self.exec.threads, count)
            seqs_per_file = count / number_of_files
            paths = [Path(base_file.parent, f'{base_file.name}.{i}') for i in range(number_of_files)]
            filehandles = [open(paths[i], 'w') for i in range(number_of_files)]
            count = 0
            for contig in self.genome.contigs.values():
                if 'CDS' == target:
                    for f in contig.features:
                        if f.type == 'CDS':
                            SeqIO.write(f.make_biopython_record(), filehandles[int(count / seqs_per_file)], "fasta")
                            count += 1
                else:
                    SeqIO.write(contig.make_biopython_record(False), filehandles[int(count / seqs_per_file)], "fasta")
                    count += 1
            for f in filehandles:
                f.close()
            return paths
