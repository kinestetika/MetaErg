from pathlib import Path
from math import log10

from metaerg.datatypes import sqlite
from metaerg.datatypes.blast import DBentry, BlastHit, BlastResult, taxon_at_genus
from metaerg import context

MAX_BLAST_HITS = 10
COLORS = 'id=cr id=cr id=co id=cb id=cg'.split()


@context.register_html_writer
def write_html(genome_name, db_connection, genome_properties:dict, dir):
    """Writes a html file for each feature to dir <file>"""
    dir = Path(dir, genome_name)
    dir.mkdir(exist_ok=True, parents=True)
    file = dir / 'feature-details.html'

    feature_html = ''
    for feature in sqlite.read_all_features(db_connection, type=('CDS', 'rRNA', 'ncRNA', 'retrotransposon')):
        if feature.type in ('CDS', 'rRNA', 'ncRNA', 'retrotransposon'):
            feature_html += make_feature_html(feature, genome_properties['dominant taxon'])

    html = _make_html_template().replace('FEATURE_HTML', feature_html)
    with open(file, 'w') as handle:
        handle.write(html)


def make_blast_table_html(blast_result: BlastResult, f_length, dominant_taxon, include_id, title: str,
                          headers: str ='percent id|query align|hit align|description|taxon') -> str:
    if blast_result:
        html = f'<h3>{title}</h3>\n'
        html += '<table><thead><tr>\n'
        columns = headers.split('|')
        for column in columns:
            if column == 'description':
                html += f'<th id=al>{column}</th>\n'
            else:
                html += f'<th>{column}</th>\n'
        html += '</tr><thead>\n<tbody>\n'
        for h in blast_result.hits[:MAX_BLAST_HITS]:
            if 'evalue' == columns[0]:
                percent_id = f'{h.evalue:.1e}'
                color_value = len(COLORS)-1 if not h.evalue else min(0, max(int(-log10(h.evalue)/20),len(COLORS)-1))
                percent_id_color = COLORS[color_value]
            else:
                percent_id = int(h.percent_id) if h.percent_id else ''
                percent_id_color = COLORS[min(int(h.percent_id/20), len(COLORS)-1)]
            if include_id:
                descr = h.hit.accession + ' ' + h.hit.descr
            else:
                descr = h.hit.descr
            html += '''<tr>
            <td {}>{}</td>
            <td {}>{}</td>
            <td {}>{}</td>
            <td id=al>{}</td>
            <td {}>{}</td>
            </tr>'''.format(
                percent_id_color,
                percent_id,

                COLORS[min(int(100 * h.aligned_length/f_length / 20), len(COLORS)-1)],
                int(100 * min(1.0, h.aligned_length/f_length)),

                COLORS[min(int(100 * h.aligned_length/h.hit.length / 20), len(COLORS)-1)],
                int(100 * min(1.0, h.aligned_length/h.hit.length)),

                descr,

                COLORS[int(len(COLORS) * len(set(h.hit.taxon.split()) & set(dominant_taxon.split()))
                           / (len (h.hit.taxon.split()) + 1))],
                taxon_at_genus(h.hit.taxon)
            )
        html += '</tbody>\n'
        html += '</table>\n'
        return html
    else:
        return ''


def make_feature_html(f, dominant_taxon) -> str:
    html = _make_feature_html_template()
    html = html.replace('FEATURE_ID', f.id)
    if f.taxon:
        html = html.replace('HEADER', f'>{f.id} {f.descr} [{taxon_at_genus(f.taxon)}]')
    else:
        html = html.replace('HEADER', f'>{f.id} {f.descr}')
    html = html.replace('SEQUENCE', f.aa_seq if 'CDS' == f.type else f.nt_seq)
    if f.type == 'CDS':
        length = (f.end - f.start) // 3
        length_unit = 'aa'
    else:
        length = f.end - f.start
        length_unit = 'nt'
    html = html.replace('BLAST_TABLE', make_blast_table_html(f.blast, length, dominant_taxon, include_id = False,
                                                             title='Top Blast Hits'))
    html = html.replace('CDD_TABLE', make_blast_table_html(f.cdd, length, dominant_taxon, include_id = True,
                                                           title='Top Conserved Domain Database Hits',
                                                           headers='percent id|query align|hit align|description| '))
    html = html.replace('HMM_TABLE', make_blast_table_html(f.hmm, length, dominant_taxon, include_id = True,
                                                           title='Top Functional Gene HMM Hits',
                                                           headers='evalue|query align|hit align|description| '))
    attribute_html = '<table>\n'
    f.tmh = '' if not f.tmh else f.tmh
    f_as_dict = {k:v for k,v in f}
    attribute_html += ''.join(f'<tr><td id=al>{k}</td><td id=al>{f_as_dict[k]}</td></tr>\n' for k in
                              ('start', 'end', 'strand', 'type', 'inference', 'subsystems', 'descr', 'taxon', 'notes',
                               'antismash', 'signal_peptide', 'tmh', 'tmh_topology'))
    attribute_html += f'<tr><td id=al>length</td><td id=al>{length} {length_unit}</td></tr>\n'
    attribute_html += '</table>\n'
    html = html.replace('ATTRIBUTE_TABLE', attribute_html)
    return html


def _make_html_template() -> str:
    return '''<!DOCTYPE html>
<html>
<head>
    <meta charset="utf-8">
    <title>Gene Details</title>
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
        #gn {
          background-color: gainsboro;
            }
    </style>
    FEATURE_HTML
</body>
</html>
'''


def _make_feature_html_template() -> str:
    return '''
    <div id=gn><h3 id="FEATURE_ID">HEADER</h3></div>
    <div>SEQUENCE</div>
    <div id=f>
    BLAST_TABLE
    CDD_TABLE
    HMM_TABLE
    <h3>Other attributes</h3>
    ATTRIBUTE_TABLE
    </div>
'''