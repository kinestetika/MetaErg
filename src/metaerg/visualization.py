from Bio import SeqIO
import pandas as pd


def open_genome(gbk_file, contig_dict, contig_list):
    features = pd.DataFrame(columns=['contig', 'id', 'start', 'end', 'strand', 'type', 'product', 'taxon', 'sequence' 
                                     'inference', 'note', 'cdd', '#TMH', 'TMH', 'signal peptide', 'rfam profile',
                                     'antismash_region', 'antismash_function', 'antismash_gene'])

    with open(gbk_file) as handle:
        for gb_record in SeqIO.parse(handle, "genbank"):
            contig_dict[gb_record.id] = gb_record
            contig_list.add(gb_record)
