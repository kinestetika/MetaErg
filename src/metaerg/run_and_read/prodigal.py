import re

from metaerg.data_model import Genome, SeqRecord, FeatureType, SeqFeature
from metaerg import context
from metaerg import bioparsers


def _run_programs(genome:Genome, result_files):
    fasta_file = context.spawn_file('masked', genome.id)
    bioparsers.write_genome_to_fasta_files(genome, fasta_file, mask=True)
    context.run_external(f'prodigal -g {genome.translation_table} -m -f gff -q -i {fasta_file} -a {result_files[0]}')


def _read_results(genome:Genome, result_files) -> int:
    cds_found = 0
    ORF_ID_PATTERN = re.compile(r'_(\d+?)$')
    with bioparsers.FastaParser(result_files[0], cleanup_seq=False) as fasta_reader:
        for seq_rec in fasta_reader:
            words = seq_rec.descr.split('#')
            try:
                m = ORF_ID_PATTERN.search(seq_rec.id)
                contig_id = seq_rec.id[0:m.start()]
                contig: SeqRecord = genome.contigs[contig_id]
            except KeyError:
                context.log(f'({genome.id}) Warning: Failed to find contig with "{seq_rec.id}"')
                continue
            start = int(words[1].strip()) - 1
            end =  int(words[2].strip())
            strand = int(words[3].strip())
            if seq_rec.seq.endswith('*'):
                seq_rec.seq = seq_rec.seq[:-1]
            feature = SeqFeature(start, end, strand, FeatureType.CDS, inference='prodigal', seq=seq_rec.seq)
            if 'partial=01' in seq_rec.descr or 'partial=01' in seq_rec.descr or 'partial=11' in seq_rec.descr:
                feature.notes.add('partial protein')
            contig.features.append(feature)
            cds_found += 1
    return cds_found


@context.register_annotator
def run_and_read_prodigal():
    return ({'pipeline_position': 61,
             'purpose': 'coding sequence prediction with prodigal',
             'programs': ('prodigal',),
             'result_files': ('prodigal',),
             'run': _run_programs,
             'read': _read_results})
