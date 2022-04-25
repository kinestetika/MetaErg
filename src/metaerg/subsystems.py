import re

from metaerg import utils
from metaerg import databases

SUBSYSTEM_LIST = []
SUBSYSTEM_HASH = {}

SUBSYSTEMS = '''>[secondary-metabolites]
# (populated by antismash)
>[ribosome]
ribosomal protein L1
ribosomal protein L2
ribosomal protein L3
ribosomal protein L4
ribosomal protein L5
ribosomal protein L6
ribosomal protein L9
ribosomal protein L11
ribosomal protein L11 methyltransferase
ribosomal protein L13
ribosomal protein L14
ribosomal protein L15
ribosomal protein L16
ribosomal protein L18
ribosomal protein L19
ribosomal protein L20
ribosomal protein L21
ribosomal protein L22
ribosomal protein L23
ribosomal protein L24
ribosomal protein L25
ribosomal protein L27
ribosomal protein L29
ribosomal protein L31
ribosomal protein L34
ribosomal protein L35
ribosomal protein S1
ribosomal protein S2
ribosomal protein S3
ribosomal protein S4
ribosomal protein S5
ribosomal protein S6
ribosomal protein S7
ribosomal protein S8
ribosomal protein S9
ribosomal protein S10
ribosomal protein S11
ribosomal protein S12
ribosomal protein S13
ribosomal protein S14
ribosomal protein S15
ribosomal protein S16
ribosomal protein S17
ribosomal protein S18
ribosomal protein S19
ribosomal protein S20
16S ribosomal RNA
23S ribosomal RNA
5S ribosomal RNA
ribosomal protein PSRP-3
small ribosomal subunit biogenesis GTPase RsgA
ribosomal protein S12 methylthiotransferase RimO
ribosomal protein S6--L-glutamate ligase
ribosome-binding factor RbfA
'''

def prep_subsystems():
    current_subsystem = None
    for line in SUBSYSTEMS.split('\n'):
        line = line.strip()
        if line.startswith("#") or not len(line):
            continue
        if line.startswith(">"):
            current_subsystem = line[1:]
            SUBSYSTEM_HASH[current_subsystem] = []
            continue
        SUBSYSTEM_HASH[current_subsystem].append(line)
        SUBSYSTEM_LIST.append((line, current_subsystem))
    SUBSYSTEM_LIST.sort(key=lambda x: -len(x[0]))
    utils.log(f'Prepared subsystem database with {len(SUBSYSTEM_HASH)} subsystems.')


def get_empty_subsystem_hash():
    h = {}
    for key in SUBSYSTEM_HASH:
        if len(SUBSYSTEM_HASH[key]):
            h[key] = dict()
            for phrase in SUBSYSTEM_HASH[key]:
                h[key][phrase] = list()
        else:
            h[key] = list()
    return h


def add_subsystem_to_feature(feature, subsystem, phrase, assignments):
    feature_id = utils.get_feature_qualifier(feature, 'id')
    prev = utils.get_feature_qualifier(feature, 'subsystem')
    if prev:
        utils.set_feature_qualifier(feature, 'subsystem', f'{prev}, {subsystem}')
    else:
        utils.set_feature_qualifier(feature, 'subsystem', subsystem)
    if phrase:
        assignments[subsystem][phrase].append(feature_id)
    else:
        assignments[subsystem].append(feature_id)


def match_subsystem(descr):
    for subsystem_entry in SUBSYSTEM_LIST:
        gene_phrase = subsystem_entry[0]
        if len(descr) > len(gene_phrase) + 20:
            return None, None
        match = re.search(r'\b' + gene_phrase + r'\b', descr)
        # match = re.search(gene_phrase, descr)
        if match and match.start() < 10:
            return subsystem_entry[1], gene_phrase
    return None, None


def match_feature_to_subsystems(feature, blast_results, subsystem_assignments):
    feature_id = utils.get_feature_qualifier(feature, 'id')
    for topic in ['cdd', 'diamond', 'blastn']:
        try:
            blast_result = blast_results[topic][feature_id]
            for hit in blast_result:
                # determine if blast alignment is good
                if topic == 'cdd':
                    cdd_id = int(hit["hit_id"][4:])
                    cdd_item = databases.CDD[cdd_id]
                    aligned = abs(hit["hit_start"] - hit["hit_end"]) / cdd_item[3]
                    hit_descr = cdd_item[2]
                else:
                    db_entry = databases.decipher_database_id(hit['hit_id'])
                    aligned = abs(hit["hit_start"] - hit["hit_end"]) / db_entry["length"]
                    hit_descr = db_entry["descr"]
                if aligned < 0.8:
                    continue
                subsystem_matched, phrase = match_subsystem(hit_descr)
                if subsystem_matched:
                    add_subsystem_to_feature(feature, subsystem_matched, phrase, subsystem_assignments)
                    return

        except KeyError:
            continue


def get_subsystem_stats(subsystem_hash):
    if isinstance(subsystem_hash, list):
        return 0, len(subsystem_hash), 1
    else:
        subsystem_gene_count = 0
        for phrase in subsystem_hash.keys():
            if len(subsystem_hash[phrase]):
                subsystem_gene_count += 1
        return len(subsystem_hash), subsystem_gene_count, subsystem_gene_count / len(subsystem_hash)

