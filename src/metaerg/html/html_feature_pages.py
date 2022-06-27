from pathlib import Path
import pandas as pd
from metaerg.datatypes import blast
from metaerg import context


@context.register_html_writer
def write_html(genome_name, feature_data: pd.DataFrame, dir):
    """Writes a html file for each feature to dir <file>"""
    dir = Path(dir, genome_name, 'features')
    for f in feature_data.itertuples():
        if f.type in ('CDS', 'rRNA', 'ncRNA', 'retrotransposon'):
            f_filename = Path(dir, f'{f.id}.html')
            with open(f_filename, 'w') as handle:
                handle.write(make_feature_html(f, dominant_taxon))


def make_blast_table_html(blast_result: blast.BlastResult, f_length, dominant_taxon) -> str:
    colors = 'id=cr id=cr id=co id=cb id=cg'.split()
    if blast_result:
        html = '<table>'
        html += '<thead><tr> '
        for column in ('percent id', 'query align', 'hit align', 'description', 'taxon'):
            if column == 'description':
                html += f'<th id=al>{column}</th>\n'
            else:
                html += f'<th>{column}</th>\n'
        html += '</tr><thead>\n<tbody>\n'
        for h in blast_result.hits:
            html += '''<tr>
            <td {}>{}</td>
            <td {}>{}</td>
            <td {}>{}</td>
            <td id=al>{}</td>
            <td {}>{}</td>
            </tr>'''.format(
                colors[min(int(h.percent_id/20), len(colors)-1)],
                int(h.percent_id),

                colors[min(int(100 * h.aligned_length/f_length / 20), len(colors)-1)],
                int(100 * min(1.0, h.aligned_length/f_length)),

                colors[min(int(100 * h.aligned_length/h.hit.length / 20), len(colors)-1)],
                int(100 * min(1.0, h.aligned_length/h.hit.length)),

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


def make_feature_html(f, dominant_taxon) -> str:
    html = _make_html_template()
    html = html.replace('FEATURE_ID', f.id)
    if f.taxon:
        html = html.replace('HEADER', f'>{f.id} {f.descr} [{f.taxon_at_genus()}]')
    else:
        html = html.replace('HEADER', f'>{f.id} {f.descr}')
    html = html.replace('SEQUENCE', f.seq)
    if f.type is 'CDS':
        length = len(f) / 3
    else:
        length = len(f)
    html = html.replace('BLAST_TABLE', make_blast_table_html(f.blast, length, dominant_taxon))
    html = html.replace('CDD_TABLE', make_blast_table_html(f.cdd, length, dominant_taxon))
    attribute_html = '<table>\n'
    attribute_html += ''.join(f'<tr><td id=al>{k}</td><td id=al>{f.__dict__[k]}</td></tr>\n' for k in
                              ('start', 'end', 'strand', 'type', 'inference', 'subsystems', 'descr', 'taxon', 'notes',
                               'antismash ', 'signal_peptide', 'tmh', 'tmh_topology'))
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
        h3 {text-align: left; font-family: Calibri, sans-serif; margin: 10px 0px 0 px 0px;}
        #al {text-align: left;}
        #cg {color: green;}
        #cr {color: red;}
        #cb {color: blue;}
        #co {color: orange;}
        th {
          background-color: white; margin: 0px 10px 0 px 0px;
            }
        #al {
          text-align: left;
            }
        #f {
          font-family: Calibri, sans-serif;
          text-align: center;
          padding: 0px;
          margin: 0px;
           }
    </style>
    <h3>HEADER</h3>
    <div>SEQUENCE</div>
    <div id=f>
    <h3>Top Blast Hits</h3>
    BLAST_TABLE
    <h3>Top Conserved Domain Database Hits</h3>
    CDD_TABLE
    <h3>Other attributes</h3>
    ATTRIBUTE_TABLE
    </div>
</body>
</html>
'''