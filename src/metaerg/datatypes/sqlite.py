import sqlite3 as sql

from metaerg.datatypes.blast import BlastResult, BlastHit, DBentry

FEATURE_FIELDS = tuple('id genome contig start end strand type inference subsystems descr taxon notes ' \
                    'aa_seq nt_seq antismash signal_peptide tmh tmh_topology blast cdd hmm'.split())
RNA_TARGETS = set("rRNA tRNA tmRNA ncRNA retrotransposon".split())
SQLITE_CREATE_TABLE_SYNTAX = '''CREATE TABLE features(
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
SQLITE_UPDATE_SYNTAX = '''UPDATE features SET
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
        self.subsystems = subsystems
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

def feature_factory(cursor, row) -> Feature:
    fields = [column[0] for column in cursor.description]
    if len(fields) < 10:
        return fields[0]
    return Feature(**{key: value for key, value in zip(fields, row)})

def connect_to_db(sql_db_file):
    connection = sql.connect(sql_db_file)
    connection.row_factory = feature_factory
    return connection

def create_db(sql_db_file):
    sql_db_file.unlink(missing_ok=True)
    connection = sql.connect(sql_db_file)
    cursor = connection.cursor()
    cursor.execute(SQLITE_CREATE_TABLE_SYNTAX)

def add_new_feature_to_db(sql_connection, feature: Feature):
    feature_as_tuple = tuple(str(v) for k, v in feature)
    cursor = sql_connection.cursor()
    cursor.execute('INSERT INTO features VALUES(?,?,?,?, ?,?,?,?, ?,?,?,?, ?,?,?,?, ?,?,?,?, ?)', feature_as_tuple)
    sql_connection.commit()

def update_feature_in_db(sql_connection, feature: Feature):
    cursor = sql_connection.cursor()
    feature_as_list = list(str(v) for k, v in feature)
    feature_as_list.append(feature.rowid)
    cursor.execute(SQLITE_UPDATE_SYNTAX, tuple(feature_as_list))
    sql_connection.commit()

#def update_feature_id_in_db(sql_connection, feature: Feature):
#    cursor = sql_connection.cursor()
#    cursor.execute('UPDATE features SET id = ? WHERE genome = ? AND contig = ? AND start = ? AND end = ? AND '
#                   'strand = ? AND type = ? AND inference = ?',
#                   (feature.id, feature.genome, feature.contig, feature.start, feature.end, feature.strand,
#                    feature.type, feature.inference))
#    sql_connection.commit()

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
