from pathlib import Path

from metaerg.data_model import FeatureType, SeqFeature, BlastResult, Genome
from metaerg import context


@context.register_html_writer
def write_html(genome: Genome, dir):
    """Writes a html file for each feature to dir <file>"""
    dir = Path(dir, genome.id, 'features')
    dir.mkdir(exist_ok=True, parents=True)
    for c in genome.contigs.values():
        for f in c.features:
            if f.type in (FeatureType.CDS, FeatureType.rRNA, FeatureType.ncRNA, FeatureType.retrotransposon):
                f_filename = Path(dir, f'{f.id}.html')
                with open(f_filename, 'w') as handle:
                    handle.write(make_feature_html(f, genome.properties['dominant taxon']))


def make_blast_table_html(blast_result: BlastResult, f_length, dominant_taxon) -> str:
    colors = 'id=cr id=cr id=co id=cb id=cg'.split()
    if blast_result:
        html = '<table>'
        html += '<thead><tr> '
        for column in ('percent id', 'query align' 'hit align', 'description', 'taxon'):
            if column == 'description':
                html += f'<th id=al>{column}</th>\n'
            else:
                html += f'<th>{column}</th>\n'
        html += '</tr><thead>\n<tbody>\n'
        for h in blast_result.hits:
            html += '<tr><td {}>{}</td {}><td {}>{}</td><td>{}</td><td>{}</td><td {}>{}</td></tr>'.format(
                colors[min(int(h.percent_id/20), len(colors)-1)],
                int(h.percent_id),
                colors[min(int(100 * h.aligned_length/f_length / 20), len(colors)-1)],
                int(100 * h.aligned_length/f_length),
                colors[min(int(100 * h.aligned_length/h.hit.length / 20), len(colors)-1)],
                int(100 * h.aligned_length/h.hit.length),
                h.hit.descr,
                colors[int(len(colors) * len(set(h.hit.taxon.split()) & set(dominant_taxon.split()))
                           / (len (h.hit.taxon.split()) + 1))],
                h.hit.taxon_at_genus()
            )
        html += '</tbody>\n'
        html += '</table>\n'
        return html
    else:
        return ''


def make_feature_html(f: SeqFeature, dominant_taxon) -> str:
    html = _make_html_template()
    html = html.replace('FEATURE_ID', f.id)
    if f.taxon:
        html = html.replace('HEADER', f'>{f.id} {f.descr} [{f.taxon_at_genus()}]')
    else:
        html = html.replace('HEADER', f'>{f.id} {f.descr}')
    html = html.replace('SEQUENCE', f.seq)
    html = html.replace('BLAST_TABLE', make_blast_table_html(f.blast, len(f), dominant_taxon))
    html = html.replace('CDD_TABLE', make_blast_table_html(f.cdd, len(f), dominant_taxon))
    attribute_html = '<table>\n'
    attribute_html += ''.join(f'<tr><td>{k}</td><td>{f.__dict__[k]}</td></tr>\n' for k in
                              SeqFeature.displayed_keys)
    attribute_html += '</table>\n'
    html = html.replace('ATTRIBUTE_TABLE', attribute_html)
    return html


def _make_html_template() -> str:
    return '''<!DOCTYPE html>
<html>
<head>
    <meta charset="utf-8">
    <title>FEATURE_ID</title>
</head>
<body>
    <style>
        div {font-family: Calibri, sans-serif; word-wrap: break-word;}
        h3 {font-family: Calibri, sans-serif; margin: 0px;}
        #al {text-align: left;}
        #cg {color: green;}
        #cr {color: red;}
        #cb {color: blue;}
        #co {color: orange;}
    </style>
    <h3>HEADER</h3>
    <div>SEQUENCE</div>
    BLAST_TABLE
    CDD_TABLE
    ATTRIBUTE_TABLE
</body>
</html>
'''