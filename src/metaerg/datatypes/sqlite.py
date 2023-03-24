import sqlite3 as sql
from metaerg.datatypes.functional_genes import FunctionalGene
from metaerg.datatypes.blast import DBentry, BlastHit, BlastResult


#FEATURE_FIELDS = tuple('id genome contig start end strand type inference subsystems descr taxon notes ' \
#                    'aa_seq nt_seq antismash signal_peptide tmh tmh_topology blast cdd hmm'.split())

RNA_TARGETS = set("rRNA tRNA tmRNA ncRNA retrotransposon".split())
SQLITE_CREATE_FEATURE_TABLE_SYNTAX = '''CREATE TABLE features(
    id TEXT,
    genome TEXT,
    contig TEXT,
    start INT,
    end INT,
    strand INT,
    type TEXT,
    inference TEXT,
    subsystems TEXT,
    descr TEXT,
    taxon TEXT,
    notes TEXT,
    aa_seq TEXT,
    nt_seq TEXT,
    antismash TEXT,
    signal_peptide TEXT,
    tmh INT,
    tmh_topology TEXT,
    blast TEXT,
    cdd TEXT,
    hmm TEXT
)'''
SQLITE_UPDATE_FEATURE_SYNTAX = '''UPDATE features SET
    id = ?,
    genome = ?,
    contig = ?,
    start = ?,
    end = ?,
    strand = ?,
    type = ?,
    inference = ?,
    subsystems = ?,
    descr = ?,
    taxon = ?,
    notes = ?,
    aa_seq = ?,
    nt_seq = ?,
    antismash = ?,
    signal_peptide = ?,
    tmh = ?,
    tmh_topology = ?,
    blast = ?,
    cdd = ?,
    hmm = ? 
WHERE rowid = ?'''

class Feature:
    def __init__(self,
                 rowid: int = 0,
                 id: str = '',
                 genome: str = '',
                 contig: str = '',
                 start: int = 0,
                 end: int = 0,
                 strand: int = 0,
                 type: str = '',
                 inference: str = '',
                 subsystems: str = '',
                 descr: str = '',
                 taxon: str = '',
                 notes: str = '',
                 aa_seq: str = '',
                 nt_seq: str = '',
                 antismash: str = '',
                 signal_peptide: str = '',
                 tmh: int = 0,
                 tmh_topology: str = '',
                 blast: str = '',
                 cdd: str = '',
                 hmm: str = ''):
        self.rowid = rowid
        self.id = id
        self.genome = genome
        self.contig = contig
        self.start = start
        self.end = end
        self.strand = strand
        self.type = type
        self.inference = inference
        self.subsystems = eval(subsystems) if subsystems else []
        self.descr = descr
        self.taxon = taxon
        self.notes = notes
        self.aa_seq = aa_seq
        self.nt_seq = nt_seq
        self.antismash = antismash
        self.signal_peptide = signal_peptide
        self.tmh = tmh
        self.tmh_topology = tmh_topology
        self.blast = eval(blast) if blast else None
        self.cdd = eval(cdd) if cdd else None
        self.hmm = eval(hmm) if hmm else None

    def __iter__(self):
        for k, v in self.__dict__.items():
            if k == 'rowid':
                continue
            yield k, v

    def __repr__(self):
        return '{}({})'.format(type(self).__name__, ', '.join(f'{k}={v!r}' for k, v in self))

    def __str__(self):
        return '{}:\n{}'.format(type(self).__name__, '\n'.join(f'  {k:14}: {v}' for k, v in self))

    def __len__(self):
        return len(self.__dict__.items())

    def length_nt(self):
        return self.end - self.start

    def length_aa(self):
        return (self.end - self.start) // 3

SQLITE_CREATE_GENOME_TABLE_SYNTAX = '''CREATE TABLE genomes(
    name TEXT,
    input_fasta_file TEXT,
    size INT,
    number_of_contigs INT,
    fraction_gc FLOAT,
    n50_contig_length INT,
    fraction_complete FLOAT,
    fraction_contaminated FLOAT,    
    number_of_features INT,
    number_of_proteins INT,
    number_of_ribosomal_rna INT,
    number_of_transfer_rna INT,
    number_of_transfer_messenger_rna INT,
    number_of_noncoding_rna INT,
    number_of_retrotransposons INT,
    number_of_crispr_repeats INT,
    number_of_other_repeats INT,
    fraction_coding FLOAT,
    fraction_repeats FLOAT,
    genetic_code INT,
    mean_protein_length FLOAT,
    top_taxon TEXT,
    fraction_classified FLOAT,
    fraction_classified_to_top_taxon FLOAT,
    codon_usage_bias FLOAT,
    doubling_time FLOAT,
    subsystems TEXT,
    subsystem_summary TEXT
)'''
SQLITE_UPDATE_GENOME_SYNTAX = '''UPDATE genomes SET
    name = ?,
    input_fasta_file = ?,
    size = ?,
    number_of_contigs = ?,
    fraction_gc = ?,
    n50_contig_length = ?,
    fraction_complete = ?,
    fraction_contaminated = ?
    number_of_features = ?,
    number_of_proteins = ?,
    number_of_ribosomal_rna = ?,
    number_of_transfer_rna = ?,
    number_of_transfer_messenger_rna = ?,
    number_of_noncoding_rna = ?,
    number_of_retrotransposons = ?,
    number_of_crispr_repeats = ?,
    number_of_other_repeats = ?,
    fraction_coding = ?,
    fraction_repeats = ?,
    genetic_code = ?,
    mean_protein_length = ?,
    top_taxon = ?,
    fraction_classified = ?,
    fraction_classified_to_top_taxon = ?,
    codon_usage_bias = ?,
    doubling_time = ?,
    subsystems = ?,
    subsystem_summary = ?
WHERE rowid = ?'''


GENOME_FORMATS = {'genome name': '<',
                  'input fasta file': '<',
                  '# contigs': ',',
                  'size': ',',
                  '% GC': '.2%',
                  'N50 contig length': ',',
                  '# proteins': ',',
                  '% coding': '.1%',
                  'mean protein length (aa)': '.0f',
                  '# ribosomal RNA': ',',
                  '# transfer RNA': ',',
                  '# transfer-messenger RNA'
                  '# non-coding RNA': ',',
                  '# retrotransposons': ',',
                  '# CRISPR repeats': ',',
                  '# other repeats': ',',
                  '# total features': ',',
                  '% repeats': '.1%',
                  'classification (top taxon)': '<',
                  '% of CDS classified to top taxon': '.1%',
                  '% of CDS that could be classified': '.1%',
                  'codon usage bias': '.3f',
                  'doubling_time (days)': '.1f'
                   }
class Genome:
    def __init__(self, rowid=0, name='', input_fasta_file='', size=0, number_of_contigs=0, fraction_gc=0.0,
                 n50_contig_length=0, fraction_complete=0.0, fraction_contaminated=0.0, number_of_features=0,
                 number_of_proteins=0, number_of_ribosomal_rna=0, number_of_transfer_rna=0,
                 number_of_transfer_messenger_rna=0, number_of_noncoding_rna=0, number_of_retrotransposons=0,
                 number_of_crispr_repeats=0, number_of_other_repeats=0, fraction_coding=0.0, fraction_repeats=0.0,
                 genetic_code=0, mean_protein_length=0, top_taxon='', fraction_classified=0.0,
                 fraction_classified_to_top_taxon=0.0, codon_usage_bias=0.0, doubling_time=0.0,
                 subsystems=None, subsystem_summary=None):
        self.rowid = rowid
        self.name = name
        self.input_fasta_file = input_fasta_file
        self.size = size
        self.number_of_contigs = number_of_contigs
        self.fraction_gc = fraction_gc
        self.n50_contig_length = n50_contig_length
        self.fraction_complete = fraction_complete
        self.fraction_contaminated = fraction_contaminated
        self.number_of_features = number_of_features
        self.number_of_proteins = number_of_proteins
        self.number_of_ribosomal_rna = number_of_ribosomal_rna
        self.number_of_transfer_rna = number_of_transfer_rna
        self.number_of_transfer_messenger_rna = number_of_transfer_messenger_rna
        self.number_of_noncoding_rna = number_of_noncoding_rna
        self.number_of_retrotransposons = number_of_retrotransposons
        self.number_of_crispr_repeats = number_of_crispr_repeats
        self.number_of_other_repeats = number_of_other_repeats
        self.fraction_coding = fraction_coding
        self.fraction_repeats = fraction_repeats
        self.genetic_code = genetic_code
        self.mean_protein_length = mean_protein_length
        self.top_taxon = top_taxon
        self.fraction_classified = fraction_classified
        self.fraction_classified_to_top_taxon = fraction_classified_to_top_taxon
        self.codon_usage_bias = codon_usage_bias
        self.doubling_time = doubling_time
        self.subsystems = eval(subsystems) if subsystems else {}
        self.subsystem_summary = eval(subsystem_summary) if subsystem_summary else {}


    def __iter__(self):
        for k, v in self.__dict__.items():
            if k == 'rowid':
                continue
            yield k, v

    def __repr__(self):
        return '{}({})'.format(type(self).__name__, ', '.join(f'{k}={v!r}' for k, v in self))

    def __str__(self):
        return '{}:\n{}'.format(type(self).__name__, '\n'.join(f'  {k:26}: {v}' for k, v in self))

    def to_dict_pretty(self) -> dict[str, str]:
        property_list = {}
        for k, v in self:
            if 'subsystem' in k:
                continue
            k = k.replace('_', ' ')
            if 'fraction' in k:
                property_list[k.replace('fraction', '%')] = f'{v:.1%}'
            elif isinstance(v, int):
                property_list[k] = f'{v:,}'
            elif isinstance(v, float):
                property_list[k] = f'{v:.2f}'
            else:
                property_list[k] = v
        return property_list

    def __len__(self):
        return len(self.__dict__.items())


def feature_factory(cursor, row) -> Feature:
    fields = [column[0] for column in cursor.description]
    if len(fields) < 10:
        return fields[0]
    return Feature(**{key: value for key, value in zip(fields, row)})

def genome_factory(cursor, row) -> Genome:
    fields = [column[0] for column in cursor.description]
    if len(fields) < 10:
        return fields[0]
    return Genome(**{key: value for key, value in zip(fields, row)})


# def connect_to_db(sql_db_file, target='Features'):
#     connection = sql.connect(sql_db_file)
#     if 'Features' == target:
#         connection.row_factory = feature_factory
#     elif 'Genomes' == target:
#         connection.row_factory = genome_factory
#     return connection

def create_db(target='Features'):
    #sql_db_file.unlink(missing_ok=True)
    connection = sql.connect(':memory:')  # sql_db_file
    cursor = connection.cursor()
    if 'Features' == target:
        cursor.execute(SQLITE_CREATE_FEATURE_TABLE_SYNTAX)
        connection.row_factory = feature_factory
    elif 'Genomes' == target:
        cursor.execute(SQLITE_CREATE_GENOME_TABLE_SYNTAX)
        connection.row_factory = genome_factory
    return connection


def write_db(sql_connection, db_file):
    db_file.unlink(missing_ok=True)
    db_write_connection = sql.connect(db_file)
    with db_write_connection:
        sql_connection.backup(db_write_connection)
    db_write_connection.close()


def add_new_feature_to_db(sql_connection, feature: Feature):
    feature_as_tuple = tuple(str(v) for k, v in feature)
    cursor = sql_connection.cursor()
    cursor.execute('INSERT INTO features VALUES(?,?,?,?, ?,?,?,?, ?,?,?,?, ?,?,?,?, ?,?,?,?, ?)', feature_as_tuple)
    sql_connection.commit()

def add_new_genome_to_db(sql_connection, genome: Genome):
    genome_as_tuple = tuple(str(v) for k, v in genome)
    cursor = sql_connection.cursor()
    cursor.execute('INSERT INTO genomes VALUES(?,?,?,?, ?,?,?,?, ?,?,?,?, ?,?,?,?, ?,?,?,?, ?,?,?,?, ?,?,?,?)', genome_as_tuple)
    sql_connection.commit()


def update_feature_in_db(sql_connection, feature: Feature):
    cursor = sql_connection.cursor()
    feature_as_list = list(str(v) for k, v in feature)
    feature_as_list.append(feature.rowid)
    cursor.execute(SQLITE_UPDATE_FEATURE_SYNTAX, tuple(feature_as_list))
    sql_connection.commit()


def update_genome_in_db(sql_connection, genome: Genome):
    cursor = sql_connection.cursor()
    genome_as_list = list(str(v) for k, v in genome)
    genome_as_list.append(genome.rowid)
    cursor.execute(SQLITE_UPDATE_FEATURE_SYNTAX, tuple(genome_as_list))
    sql_connection.commit()


def drop_feature(sql_connection, feature: Feature):
    cursor = sql_connection.cursor()
    cursor.execute('DELETE FROM features WHERE rowid = ?',
                   (feature.rowid,))
    sql_connection.commit()

def count_features(sql_connection):
    sql_connection.row_factory = None
    cursor = sql_connection.cursor()
    total = sum(1 for f in cursor.execute('SELECT start FROM features'))
    sql_connection.row_factory = feature_factory
    return total

def read_feature_by_id(sql_connection, feature_id) -> Feature:
    cursor = sql_connection.cursor()
    result = cursor.execute('SELECT rowid, * FROM features WHERE id = ?', (feature_id,))
    return result.fetchone()

def read_genome_by_id(sql_connection, genome_name) -> Genome:
    cursor = sql_connection.cursor()
    result = cursor.execute('SELECT rowid, * FROM genomes WHERE name = ?', (genome_name,))
    return result.fetchone()


def read_all_features(sql_connection, contig='', type=None, location=None, additional_sql=None):
    cursor = sql_connection.cursor()

    where_str = []
    fields = []
    if type:
        if isinstance(type, list) or isinstance(type, tuple)  or isinstance(type, set):
            where_str.append('type IN ({})'.format(','.join('?' for t in type)))
            fields.extend(type)
        else:
            where_str.append('type = ?')
            fields.append(type)
    if contig:
        where_str.append('contig = ?')
        fields.append(contig)
    if location:
        where_str.append('start < ? AND ? < end')
        fields.append(location[1])
        fields.append(location[0])
    if additional_sql:
        where_str.append(additional_sql)

    where_str = ' AND '.join(where_str)
    if where_str:
        #print('SELECT * FROM features WHERE ' + where_str + ' ORDER BY contig, start')
        #print(fields)
        return cursor.execute('SELECT rowid, * FROM features WHERE ' + where_str + ' ORDER BY contig, start', tuple(fields))
    else:
        return cursor.execute('SELECT rowid, * FROM features ORDER BY contig, start')
    #for columns in result.fetchall():
    #    yield Feature(*columns)

def read_all_genomes(sql_connection, sql_where=''):
    cursor = sql_connection.cursor()
    return cursor.execute('SELECT rowid, * FROM genomes ' + sql_where + ' ORDER BY name')
