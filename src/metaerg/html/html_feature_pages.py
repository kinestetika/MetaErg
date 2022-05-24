from pathlib import Path

from metaerg.run_and_read.data_model import FeatureType, MetaergSeqFeature, BlastResult, MetaergGenome
from metaerg.html.abc import HTMLwriter, register_html_writer
from metaerg.run_and_read.context import Executor


@register_html_writer
class HTMLFeaturePages(HTMLwriter):

    def __init__(self, genome, exec: Executor):
        super().__init__(genome, exec)
        self.genome: MetaergGenome = genome
        self.exec = exec


def make_blast_table_html(blast_result: BlastResult, f_length, dominant_taxon) -> str:
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
                self.get_color(h.percent_id),
                h.percent_id,
                self.get_color(100 * h.aligned_length/f_length),
                100 * h.aligned_length/f_length,
                self.get_color(100 * h.aligned_length/h.hit.length),
                100 * h.aligned_length/h.hit.length,
                h.hit.descr,
                self.get_color(20 * len(set(h.hit.taxon.split()) &
                                        set(dominant_taxon.split()))),
                h.hit.taxon_at_genus()
            )
        html += '</tbody>\n'
        html += '</table>\n'
        return html
    else:
        return ''


def make_feature_html(f: MetaergSeqFeature) -> str:
    html = _make_html_template()
    html = html.replace('FEATURE_ID', f.id)
    if f.taxon:
        html = html.replace('HEADER', f'>{f.id} {f.product} [{f.taxon_at_genus()}]')
    else:
        html = html.replace('HEADER', f'>{f.id} {f.product}')
    html = html.replace('SEQUENCE', f.sequence)
    html = html.replace('BLAST_TABLE', make_blast_table_html(f.blast, len(f)))
    html = html.replace('CDD_TABLE', make_blast_table_html(f.cdd, len(f)))
    attribute_html = '<table>\n'
    attribute_html += ''.join(f'<tr><td>{k}</td><td>{v}</td></tr>\n' for k, v in f self.genome.properties.items())
    attribute_html += '</table>\n'
    html = html.replace('ATTRIBUTE_TABLE', attribute_html)
    return html


@register_html_writer
def write_html(exec:Executor, genome: MetaergGenome, file=None):
    """Writes a html file for each feature to dir <file>"""
    if not file:
        file = Path(exec.html_dir, genome.id, "features")
    file.mkdir(exist_ok=True, parents=True)
    for c in genome.contigs.values():
        for f in c.features:
            if f.type in (FeatureType.CDS, FeatureType.rRNA, FeatureType.ncRNA, FeatureType.retrotransposon):
                f_filename = Path(file, f'{f.id}.html')
                with open(f_filename, 'w') as handle:
                    handle.write(make_feature_html(f))


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