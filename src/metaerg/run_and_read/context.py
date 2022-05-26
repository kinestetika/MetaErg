import shutil
import os
import subprocess
import time
from multiprocessing import cpu_count
from pathlib import Path
from metaerg.run_and_read.data_model import MetaergGenome

BASE_DIR = ''
TEMP_DIR = ''
HTML_DIR = ''
DATABASE_DIR = ''
CHECKM_DIR = ''
GTDBTK_DIR = ''
GENOME_NAME_MAPPING_FILE = ''
CONTIG_FILES = []
GENOME_NAMES = []

FORCE = False
MULTI_MODE = False
RENAME_CONTIGS = False
RENAME_GENOMES = False
MIN_CONTIG_LENGTH = 0
TRANSLATION_TABLE = 0
CPUS_PER_GENOME = 0
CPUS_AVAILABLE = 0
PARAlLEL_ANNOTATIONS = 0
START_TIME = 0
LOG_TOPICS = set()
FILE_EXTENSION = ''

ANNOTATOR_REGISTRY = {}
HTML_WRITER_REGISTRY = []
DATABASE_INSTALLER_REGISTRY = []

def init(contig_file, database_dir, rename_contigs, rename_genomes, min_contig_length, cpus, force, file_extension,
         translation_table, checkm_dir, gtdbtk_dir, log_topics=''):
    global BASE_DIR, TEMP_DIR, HTML_DIR, DATABASE_DIR, CHECKM_DIR, GTDBTK_DIR, GENOME_NAME_MAPPING_FILE, MULTI_MODE,\
           RENAME_CONTIGS, RENAME_GENOMES, MIN_CONTIG_LENGTH, FORCE, FILE_EXTENSION, TRANSLATION_TABLE, CPUS_PER_GENOME, \
           CPUS_AVAILABLE, START_TIME, LOG_TOPICS, PARAlLEL_ANNOTATIONS
    contig_file = Path(contig_file).absolute()
    BASE_DIR = contig_file if contig_file.is_dir() else contig_file.parent
    TEMP_DIR = Path(BASE_DIR, "temp")
    DATABASE_DIR = Path(database_dir).absolute()
    CHECKM_DIR = Path(checkm_dir).absolute()
    GTDBTK_DIR = Path(gtdbtk_dir).absolute()
    GENOME_NAME_MAPPING_FILE = Path(TEMP_DIR, 'genome.name.mapping.txt')
    HTML_DIR = Path(BASE_DIR, "html")

    MULTI_MODE = contig_file.is_dir()
    RENAME_CONTIGS = rename_contigs,
    RENAME_GENOMES = rename_genomes
    MIN_CONTIG_LENGTH = min_contig_length
    FORCE = force
    FILE_EXTENSION = file_extension
    TRANSLATION_TABLE = translation_table

    CPUS_PER_GENOME = cpus
    CPUS_AVAILABLE = cpu_count()
    START_TIME = time.monotonic()
    LOG_TOPICS = set(log_topics.split())
    log('Initialized execution environment with commenad line arguments...')

    if not contig_file.exists():
        log(f'Input file "{contig_file}" is missing. Expecting dir or a nt fasta file.')
        exit(1)
    if TEMP_DIR.exists():
        log('Warning: may overwrite existing temp files...')
        if TEMP_DIR.is_file():
            log(f'Expected folder at {TEMP_DIR}, found regular file, crash! Delete this file first')
            exit(1)
        else:
            os.mkdir(TEMP_DIR)
        shutil.rmtree(HTML_DIR, ignore_errors=True)
        os.mkdir(HTML_DIR)
        if CPUS_PER_GENOME > 0:
            CPUS_PER_GENOME = min(CPUS_PER_GENOME, CPUS_AVAILABLE)
        else:
            CPUS_PER_GENOME = CPUS_AVAILABLE
        log(f'Detected {CPUS_AVAILABLE} available threads/cpus, will use {CPUS_PER_GENOME}.')

        if contig_file.is_dir():
            CONTIG_FILES = [x.absolute() for x in sorted(contig_file.glob(f'*{FILE_EXTENSION}'))]
            if not len(CONTIG_FILES):
                log(f'Did not find any contig files with extension "{FILE_EXTENSION}" '
                          f'in dir "{contig_file}"')
                exit(1)
            PARAlLEL_ANNOTATIONS = CPUS_PER_GENOME
            CPUS_PER_GENOME = max(1, int(CPUS_PER_GENOME / len(CONTIG_FILES)))
            PARAlLEL_ANNOTATIONS = int(PARAlLEL_ANNOTATIONS / CPUS_PER_GENOME)
        else:
            CONTIG_FILES = [contig_file]
            PARAlLEL_ANNOTATIONS = 1
        if RENAME_GENOMES:
            GENOME_NAMES = [f'g{CONTIG_FILES.index(f):0>4}' for f in CONTIG_FILES]
        else:
            GENOME_NAMES = [f.name for f in CONTIG_FILES]
        log(f'writing genome names to {GENOME_NAME_MAPPING_FILE} ')
        with open(GENOME_NAME_MAPPING_FILE, 'w') as mapping_file:
            for n, o in zip(GENOME_NAMES, CONTIG_FILES):
                mapping_file.write(f'{n}\t{o.stem}\t{o}\n')
        log(f'Ready to annotate {len(CONTIG_FILES)} genomes in dir "{BASE_DIR}" with '
                  f'{CPUS_PER_GENOME} threads per genome.')
        log_settings()


def spawn_file(program_name, genome_id, base_dir = None) -> Path:
    """computes a Path genome_id.program_name or, if multimode==True, program_name/genome_id"""
    target_dir = base_dir if base_dir else TEMP_DIR
    if MULTI_MODE:
        dir = Path(target_dir, program_name)
        if not dir.exists():
            dir.mkdir()
        elif dir.is_file():
            if FORCE:
                dir.unlink()
                dir.mkdir()
            else:
                raise Exception("Use force to overwrite existing results")
        return Path(dir, genome_id)
    else:
        file = Path(target_dir, f'{genome_id}.{program_name}')
        if file.exists() and file.is_dir():
            if FORCE:
                shutil.rmtree(file)
        return file


def log(log_message, values=(), topic=''):
    if not topic or topic in LOG_TOPICS:
        if len(values):
            print(f'{format_runtime()} {log_message.format(*values)}')
        else:
            print(f'{format_runtime()} {log_message}')


def format_runtime():
    runtime = time.monotonic() - START_TIME
    return f'[{int(runtime / 3600):02d}h:{int((runtime % 3600) / 60):02d}m:{int(runtime % 60):02d}s]'


def log_settings():
    log(',\n'.join(f'{k.lower()}={eval(k)}' for k in dir() if not k.startswith('_')))


def run_external(exec, stdin=None, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, log_cmd=True):
    if log_cmd:
        log(exec)
    result = subprocess.run(exec.split(), stdout=stdout, stdin=stdin, stderr=stderr)
    if result.returncode != 0:
        raise Exception(f'Error while trying to run "{exec}"')


def sorted_annotators():
    return (ANNOTATOR_REGISTRY[a] for a in sorted[ANNOTATOR_REGISTRY.keys()])


def register_annotator(define_annotator):
    param = define_annotator()
    def annotator(genome: MetaergGenome):
        """Runs programs and reads results."""
        log('({}) {} started...', (genome.id, param['purpose']))
        # (1) First make sure that the helper programs are available:
        all_programs_in_path = True
        for p in param.get('programs', []):
            program_path = shutil.which(p, mode=os.X_OK)
            if not program_path:
                all_programs_in_path = False
                log('({}) Unable to run {}, helper program "{}" not in path', (genome.id, param['purpose'], p))
        # (2) Then, make sure required databases are available
        for d in param.get('databases', []):
            d = Path(DATABASE_DIR, d)
            if not d.exists() or not d.stat().st_size:
                log('({}) Unable to run {}, or parse results, database "{}" missing', (genome.id,
                                                                                       param['purpose'], d))
                return
        # (3) Then, if force or the results files are not yet there, run the programs:
        result_files = [spawn_file(f, genome.id) for f in param.get('result_files', [])]
        if all_programs_in_path:
            previous_results_missing = False
            for f in result_files:
                if not f.exists() or not f.stat().st_size:
                    previous_results_missing = True
                    break
            if FORCE or previous_results_missing:
                param['run'](genome, result_files)
            else:
                log('({}) Reusing existing results in {}.'.format(genome.id, result_files))
        # (4) If all results files are there, read the results:
        all_results_created = True
        for f in result_files:
            if not f.exists() or not f.stat().st_size:
                all_results_created = False
                log('({}) Missing expected result file {}.', (genome.id, f))
        if all_results_created:
            positive_count = param['read'](genome, result_files)
        else:
            positive_count = 0
        log('({}) {} complete. Found {}.', (genome.id, param['purpose'], positive_count))
        # (5) Save (intermediate) results:
        gbk_file = spawn_file("gbk", genome.id)
        gff_file = spawn_file("gff", genome.id)
        genome.write_gbk_gff(gbk_file=gbk_file, gff_file=gff_file)

    ANNOTATOR_REGISTRY[param['pipeline_position']] = annotator
    return annotator


def register_html_writer(writer):
    HTML_WRITER_REGISTRY.append(writer)
    return writer

def register_database_installer(database_installer):
    DATABASE_INSTALLER_REGISTRY.append(database_installer)
    return database_installer