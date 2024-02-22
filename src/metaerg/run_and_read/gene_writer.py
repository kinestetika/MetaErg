from metaerg import context
from metaerg.datatypes import fasta
from metaerg.datatypes import sqlite

def _run_programs(genome, contig_dict, db_connection, result_files):
    pass

def create_new_id(genome_name:str, contig_name:str, feature_number:int):
    if context.RENAME_CONTIGS:
        # contigs already contain genome name
        return context.DELIMITER.join((contig_name, f'{feature_number:05d}'))
    else:
        return context.DELIMITER.join((genome_name, contig_name, f'{feature_number:05d}'))


def _read_results(genome, contig_dict, db_connection, result_files) -> int:
    j = 0
    preliminary_id_mapping = {}
    prelim_id = ''
    for feature in sqlite.read_all_features(db_connection):
        #if feature.type == 'repeat_region':
        #    context.log(f"{feature.type}: {feature.id}")
        # some features already got an id assigned, to keep track of gene clusters
        if feature.parent:
            try:
                feature.parent = preliminary_id_mapping[feature.parent]
            except KeyError:
                context.log(f'WARNING: Unknown feature parent {feature.parent} reference for {feature.id}')
        if feature.id:
            prelim_id = feature.id
        feature.id = create_new_id(genome.name, feature.contig, j)
        preliminary_id_mapping[prelim_id] = feature.id
        j += 1
        sqlite.update_feature_in_db(db_connection, feature)
        # update genome
        genome.number_of_features += 1
        if feature.type == 'CDS':
            genome.number_of_proteins += 1
            genome.fraction_coding += feature.length_nt()
            genome.mean_protein_length += feature.length_aa()
        elif feature.type == 'rRNA':
            genome.number_of_ribosomal_rna += 1
        elif feature.type == 'tRNA':
            genome.number_of_transfer_rna += 1
        elif feature.type == 'tmRNA':
            genome.number_of_transfer_messenger_rna += 1
        elif feature.type == 'ncRNA':
            genome.number_of_noncoding_rna += 1
        elif feature.type == 'retrotransposon':
            genome.number_of_retrotransposons += 1
            genome.fraction_repeats += feature.length_nt()
        elif feature.type == 'CRISPR':
            genome.number_of_crispr_repeats += 1
            genome.fraction_repeats += feature.length_nt()
        elif feature.type == 'repeat_unit':
            genome.number_of_other_repeats += 1
            genome.fraction_repeats += feature.length_nt()

    genome.fraction_coding /= genome.size
    genome.fraction_repeats /= genome.size
    try:
        genome.mean_protein_length /= genome.number_of_proteins
    except ZeroDivisionError:
        pass
    total_rna = genome.number_of_ribosomal_rna + genome.number_of_transfer_rna + genome.number_of_noncoding_rna

    cds_file = context.spawn_file('cds.faa', genome.name)
    context.log(f'({genome.name}) Now writing {genome.number_of_proteins} proteins to fasta at {cds_file}...')
    fasta.write_features_to_fasta(db_connection, 'aa', cds_file, targets=('CDS',))
    rna_file = context.spawn_file('rna.fna', genome.name)
    context.log(f'({genome.name}) Now writing {total_rna} RNA genes and features to fasta at {rna_file}...')
    fasta.write_features_to_fasta(db_connection, 'nt', rna_file, targets=sqlite.RNA_TARGETS)
    return j


@context.register_annotator
def run_and_read_trf():
    return ({'pipeline_position': 66,
             'annotator_key': 'write_genes',
             'purpose': 'feature ID generation',
             'programs': (),
             'result_files': (),
             'run': _run_programs,
             'read': _read_results})
