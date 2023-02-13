import re
import shutil
import subprocess
import time
import os
import traceback
from multiprocessing import cpu_count
from pathlib import Path

import httpx

from metaerg import registry

BASE_DIR = Path()
TEMP_DIR = Path()
HTML_DIR = Path()
DATABASE_DIR = Path()
CHECKM_DIR = Path()
GTDBTK_DIR = Path()
GENOME_NAME_MAPPING_FILE = ''
CONTIG_FILES = []
GENOME_NAMES = []
LOG_FILE = ''

FORCE = False
MULTI_MODE = False
RENAME_CONTIGS = False
RENAME_GENOMES = False
MIN_CONTIG_LENGTH = 0
TRANSLATION_TABLE = 0
ACTIVE_ANNOTATORS = set('antismash aragorn cdd cmscan diamond_and_blastn hmm write_genes ltr_harvest minced prodigal '
                        'signalp repeat_masker tmhmm trf'.split())
READ_ONLY = False
DELIMITER = ''
CPUS_PER_GENOME = 0
CPUS_AVAILABLE = 0
PARALLEL_ANNOTATIONS = 0
START_TIME = 0
LOG_TOPICS = set()
FILE_EXTENSION = ''
METAERG_MODE_RUN = 1
METAERG_MODE_DOWNLOAD_DATABASE = 2
METAERG_MODE_CREATE_DATABASE = 3
METAERG_MODE_INSTALL_DEPS = 4
METAERG_MODE = METAERG_MODE_RUN
TASKS = 'all'
PREFIX = 'g'
BIN_DIR =''
PATH_TO_SIGNALP = None
PATH_TO_TMHMM = None

def init(contig_file, database_dir, rename_contigs, rename_genomes, min_contig_length, cpus, force, file_extension,
         translation_table, delimiter, checkm_dir, gtdbtk_dir, prefix, create_database, download_database,
         install_deps, path_to_signalp, path_to_tmhmm, log_topics='', contig_mode = False, skip_step='', read_only=True,
         run_step=''):
    global BASE_DIR, TEMP_DIR, HTML_DIR, DATABASE_DIR, CHECKM_DIR, GTDBTK_DIR, GENOME_NAME_MAPPING_FILE, MULTI_MODE,\
           RENAME_CONTIGS, RENAME_GENOMES, MIN_CONTIG_LENGTH, FORCE, FILE_EXTENSION, TRANSLATION_TABLE, \
           CPUS_PER_GENOME, CPUS_AVAILABLE, START_TIME, LOG_TOPICS, PARALLEL_ANNOTATIONS, METAERG_MODE, \
           GENOME_NAMES, CONTIG_FILES, DELIMITER, LOG_FILE, TASKS, PREFIX, BIN_DIR, PATH_TO_SIGNALP, PATH_TO_TMHMM, \
           ACTIVE_ANNOTATORS, READ_ONLY

    START_TIME = time.monotonic()
    LOG_TOPICS = set(log_topics.split())
    LOG_FILE = Path('log.txt').absolute()
    log('Initializing execution environment with command line arguments...')

    if install_deps:
        METAERG_MODE = METAERG_MODE_INSTALL_DEPS
        BIN_DIR = Path(install_deps).absolute()
        if path_to_signalp and Path(path_to_signalp).is_file():
            PATH_TO_SIGNALP = Path(path_to_signalp).absolute()
        else:
            log('Warning: path to signalp tarbal not provided or wrong; signalp will not be installed.')
        if path_to_tmhmm and Path(path_to_tmhmm).is_file():
            PATH_TO_TMHMM = Path(path_to_tmhmm).absolute()
        else:
            log('Warning: path to tmhmm tarbal not provided or wrong; tmhmm will not be installed.')
        log(f'Ready to install helper programs at {BIN_DIR} with {PATH_TO_SIGNALP}, {PATH_TO_TMHMM}.')
        return

    DATABASE_DIR = Path(database_dir).absolute()
    if not DATABASE_DIR.is_dir():
        raise Exception("No database dir provided or database dir is not a dir.")
    if download_database:
        METAERG_MODE = METAERG_MODE_DOWNLOAD_DATABASE
        log(f'Ready to download databases.')
        return
    elif create_database:
        METAERG_MODE = METAERG_MODE_CREATE_DATABASE
        if create_database == 'all':
            TASKS = 'PVEBRCSA'
        else:
            TASKS = create_database
        log(f'Ready to create databases from scratch with tasks {TASKS}.')
        return
    else: # we're going to annotate genomes...
        # (1) collect contig files
        if ',' in contig_file:
            for f in contig_file.split(','):
                f = Path(f).absolute()
                if not f.exists():
                    log(f'Did not find contig file "{f}"')
                    continue
                CONTIG_FILES.append(f)
        else:
            contig_file = Path(contig_file).absolute()
            if contig_file.is_dir():
                CONTIG_FILES = [f.absolute() for f in sorted(contig_file.glob(f'*{FILE_EXTENSION}')) if f.is_file()]
            else:
                if contig_file.exists():
                    CONTIG_FILES = [contig_file]
        if len(CONTIG_FILES):
            BASE_DIR = CONTIG_FILES[0].parent
        else:
            log(f'No contig files to annotate, did you type the correct file or file extension?')
            exit(1)
        # (2) manage file structure / folders
        TEMP_DIR = BASE_DIR / 'temp'
        CHECKM_DIR = Path(checkm_dir).absolute()
        GTDBTK_DIR = Path(gtdbtk_dir).absolute()
        GENOME_NAME_MAPPING_FILE = TEMP_DIR / 'genome.name.mapping.txt'
        HTML_DIR = BASE_DIR / 'html'
        for folder in (TEMP_DIR, HTML_DIR):
            if folder.exists():
                print(f'Warning: may overwrite existing files in {folder}...')
                if folder.is_file():
                    print(f'Expected dir at {folder}, found regular file, crash! Delete this file first')
                    exit(1)
            else:
                folder.mkdir()
        if run_step:
            ACTIVE_ANNOTATORS = {step for step in run_step.split(',')}
        if contig_mode:
            TRANSLATION_TABLE = -1
            ACTIVE_ANNOTATORS.remove('repeat_masker')
        for skipped_step in skip_step.split(','):
            try:
                ACTIVE_ANNOTATORS.remove(skipped_step)
            except KeyError:
                pass
        READ_ONLY = read_only
        # (3) set some global variables
        RENAME_CONTIGS = rename_contigs
        RENAME_GENOMES = rename_genomes
        MIN_CONTIG_LENGTH = int(min_contig_length)
        FORCE = force
        FILE_EXTENSION = file_extension
        DELIMITER = delimiter
        PREFIX = prefix
        TRANSLATION_TABLE = translation_table
        # (4) set up multitreading
        CPUS_PER_GENOME = int(cpus)
        CPUS_AVAILABLE = cpu_count() // 2
        if CPUS_PER_GENOME > 0:
            CPUS_PER_GENOME = min(CPUS_PER_GENOME, CPUS_AVAILABLE)
        else:
            CPUS_PER_GENOME = CPUS_AVAILABLE
        PARALLEL_ANNOTATIONS = CPUS_PER_GENOME
        CPUS_PER_GENOME = max(1, CPUS_PER_GENOME // len(CONTIG_FILES))
        PARALLEL_ANNOTATIONS = int(PARALLEL_ANNOTATIONS / CPUS_PER_GENOME)
        if len(CONTIG_FILES) == 1:
            PARALLEL_ANNOTATIONS = 1
        log(f'Detected {CPUS_AVAILABLE} available threads/cpus, will use {CPUS_PER_GENOME} per genome with '
            f'{PARALLEL_ANNOTATIONS} genomes annotated in parallel.')
        MULTI_MODE = len(CONTIG_FILES) > 1
        if  (TEMP_DIR / 'metaerg_progress').exists and (TEMP_DIR / 'metaerg_progress').is_dir():
            log(f'Detected previous results - maintaining file structure for compatibility.')
            MULTI_MODE = True
        # (5) manage genome and contig naming
        if RENAME_GENOMES:
            log('Will create new names for genomes.')
            RENAME_CONTIGS = True
        else:
            log('Using original filenames as genome names.')
        if RENAME_CONTIGS:
            log('Will rename contigs.')
        else:
            log('Keeping original contig identifiers.')
        # (6) more about genome names and consistency between metaerg (re)runs
        if GENOME_NAME_MAPPING_FILE.exists():
            log(f'Detected previous results -  maintaining genome name mappings for consistency.')
            prev_name_mappings = {}
            index_for_next_genome = 0
            with open(GENOME_NAME_MAPPING_FILE) as reader:
                for line in reader:
                    words = line.split()
                    prev_name_mappings[words[1]] = words[0]
                    if match := re.search(r'\d+', words[0]):
                        index = int(match.group(0))
                        if index > index_for_next_genome:
                            index_for_next_genome = index
            GENOME_NAMES = []
            for contig_file in CONTIG_FILES:
                try:
                    GENOME_NAMES.append(prev_name_mappings[contig_file.stem])
                except KeyError:
                    if RENAME_GENOMES:
                        index_for_next_genome += 1
                        new_name = f'{PREFIX}{index_for_next_genome:0>4}'
                        GENOME_NAMES.append(new_name)
                        with open(GENOME_NAME_MAPPING_FILE, 'a') as mapping_file:
                            mapping_file.write(f'{new_name}\t{contig_file.stem}\t{contig_file}\n')
                            log(f'Detected new contig file {contig_file.name} - renamed to {new_name} and added to'
                                f'name mapping file.')
                    else:
                        GENOME_NAMES.append(contig_file.stem)
        else:
            if RENAME_GENOMES:
                GENOME_NAMES = [f'{PREFIX}{CONTIG_FILES.index(f):0>4}' for f in CONTIG_FILES]
            else:
                GENOME_NAMES = [f.stem for f in CONTIG_FILES]
            log(f'Writing genome names to {GENOME_NAME_MAPPING_FILE}...')
            with open(GENOME_NAME_MAPPING_FILE, 'w') as mapping_file:
                for n, o in zip(GENOME_NAMES, CONTIG_FILES):
                    mapping_file.write(f'{n}\t{o.stem}\t{o}\n')

        log(f'Ready to annotate {len(CONTIG_FILES)} genomes in dir "{BASE_DIR}" with '
                  f'{CPUS_PER_GENOME} threads per genome and tasks {TASKS}.')

def spawn_file(program_name, genome_id, base_dir = None) -> Path:
    """computes a Path genome_id.program_name or, if multimode==True, program_name/genome_id"""
    target_dir = Path(base_dir) if base_dir else TEMP_DIR
    if MULTI_MODE:
        dir = target_dir / program_name
        dir.mkdir(exist_ok=True)
        if dir.is_file():
            if FORCE:
                dir.unlink()
                dir.mkdir(exist_ok=True)
            else:
                raise Exception("Use force to overwrite existing results")
        return dir / genome_id
    else:
        file = target_dir / f'{genome_id}.{program_name}'
        if file.exists() and file.is_dir():
            if FORCE:
                shutil.rmtree(file)
        return file


def log(log_message, values=(), topic=''):
    if not topic or topic in LOG_TOPICS:
        if len(values):
            final_msg = f'{format_runtime()} {log_message.format(*values)}'
        else:
            final_msg = f'{format_runtime()} {log_message}'
        print(final_msg)
        try:
            with open(LOG_FILE, 'a') as log_handle:
                log_handle.write(final_msg)
                log_handle.write('\n')
        except FileNotFoundError:
            pass


def format_runtime():
    runtime = time.monotonic() - START_TIME
    return f'[{int(runtime / 3600):02d}h:{int((runtime % 3600) / 60):02d}m:{int(runtime % 60):02d}s]'


def run_external(exec, stdin=None, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, log_cmd=True):
    if log_cmd:
        log(exec)
    result = subprocess.run(exec.split(), stdout=stdout, stdin=stdin, stderr=stderr)
    if result.returncode != 0:
        raise Exception(f'WARNING: "{exec}" exited with non-zero status')


def download(url: str, file: Path):
    log(f'Started downloading {url} to {file}...')
    with httpx.stream('GET', url, timeout=6.1, follow_redirects=True) as data_stream, open(file, 'wb') as handle:
        for data in data_stream.iter_bytes():
            handle.write(data)


def sorted_annotators():
    return (registry.ANNOTATOR_REGISTRY[a] for a in sorted(list(registry.ANNOTATOR_REGISTRY.keys())))


def parse_metaerg_progress(genome_name):
    progress_file = spawn_file('metaerg_progress', genome_name)
    progress = ''
    if progress_file.exists():
        with open(progress_file) as progress_reader:
            for line in progress_reader:
                progress += line
    return progress


def write_metaerg_progress(genome_name, new_progress):
    progress_file = spawn_file('metaerg_progress', genome_name)
    with open(progress_file, 'w') as writer:
        writer.write(new_progress)


def register_annotator(define_annotator):
    param = define_annotator()

    def annotator(genome_name, contig_dict, db_connection) -> int:
        """Runs programs and reads results."""
        # (1) Read metaerg progress file, update progress and log the start of the analysis
        if not param['annotator_key'] in ACTIVE_ANNOTATORS:
            return 0
        if FORCE:
            write_metaerg_progress(genome_name, '')
        current_progress = parse_metaerg_progress(genome_name)
        if not '{}=complete'.format(param['annotator_key']) in current_progress:
            current_progress += '{}=started\n'.format(param['annotator_key'])
            write_metaerg_progress(genome_name, current_progress)
        log('({}) Starting {} ...', (genome_name, param['purpose']))

        # (2) Make sure required databases are available
        for d in param.get('databases', []):
            d = DATABASE_DIR / d
            if not d.exists() or not d.stat().st_size:
                log('({}) Unable to run {}, or parse results, database "{}" missing', (genome_name,
                                                                                       param['purpose'], d))
                return 0
        # (3) Then, if force or the results files are not yet there, run the programs:
        result_files = [spawn_file(f, genome_name) for f in param.get('result_files', [])]
        analysis_already_completed = True
        for f in result_files:
            if not f.exists() or not f.stat().st_size:
                analysis_already_completed = False
                break
        # (4) check if we tried to run this analysis before without completing...
        if '{}=started'.format(param['annotator_key']) in current_progress:
            analysis_already_completed = False
        elif '{}=complete'.format(param['annotator_key']) in current_progress:
            analysis_already_completed = True

        if not READ_ONLY and (FORCE or not analysis_already_completed):
            # (5) make sure that the helper programs are available:
            all_programs_in_path = True
            p = 'X'
            for p in param.get('programs', []):
                program_path = shutil.which(p, mode=os.X_OK)
                if not program_path:
                    all_programs_in_path = False
            if all_programs_in_path:
                try:
                    param['run'](genome_name, contig_dict, db_connection, result_files)
                except Exception as e:
                    log('({}) Error while running {}: {}', (genome_name, param['purpose'],
                                                            "".join(traceback.format_exception(e))))
                    return 0
            else:
                log('({}) Unable to run {}, helper program "{}" not in path', (genome_name, param['purpose'], p))
                return 0
        elif len(param.get('programs', [])):  # before logging this, do check if any result files will be parsed
            log('({}) Reusing existing results in {}.'. format(genome_name,
                                                              ', '.join(str(file) for file in result_files)))
        # (6) Check if all results files are there, otherwise return:
        results_complete = True
        for f in result_files:
            if not f.exists() or not f.stat().st_size:
                log('({}) Missing expected result file {}; this could mean no results were found or that a '
                    'helper program failed to run.', (genome_name, f))
                results_complete = False
        # (7) Report success, update progress
        current_progress = current_progress.replace('{}=started'.format(param['annotator_key']),
                                                    '{}=complete'.format(param['annotator_key']))
        write_metaerg_progress(genome_name, current_progress)
        positive_count = 0
        if results_complete:
            positive_count = param['read'](genome_name, contig_dict, db_connection, result_files)
        log('({}) {} complete. Found {}.', (genome_name, param['purpose'], positive_count))
        return 0

    registry.ANNOTATOR_REGISTRY[param['pipeline_position']] = annotator
    return annotator


def register_html_writer(writer):
    registry.HTML_WRITER_REGISTRY.append(writer)
    return writer


def register_database_installer(database_installer):
    registry.DATABASE_INSTALLER_REGISTRY.append(database_installer)
    # print(len(registry.DATABASE_INSTALLER_REGISTRY))
    return database_installer
