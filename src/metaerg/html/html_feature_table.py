from pathlib import Path
from metaerg.run_and_read.data_model import FeatureType, MetaergGenome
from metaerg.html.abc import HTMLwriter, register_html_writer
from metaerg.run_and_read.abc import ExecutionEnvironment

@register_html_writer
class HTMLFeatureTable(HTMLwriter):

    def __init__(self, genome, exec: ExecutionEnvironment):
        super().__init__(genome, exec)
        self.genome: MetaergGenome = genome
        self.exec = exec

    def make_html(self) -> str:
        """Injects the content into the html base, returns the html."""
        html = self._make_html_template()
        html = html.replace('GENOME_NAME', self.genome.id)
        # table header
        table_headers = ''
        for column in 'id strand length type location subsystem CDD ident align recall description taxon'.split():
            if column in 'id description':
                table_headers += f'<th id=al>{column}</th>\n'
            else:
                table_headers += f'<th>{column}</th>\n'
        html = html.replace('TABLE_HEADERS', table_headers)
        # table body
        table_body = ''
        for c in self.genome.contigs.values():
            for f in c.features:
                format_hash = {'f_id': f.id,
                               'description': f.product,
                               'taxon': f.taxon_at_genus()}
                if f.type in (FeatureType.CDS, FeatureType.rRNA, FeatureType.ncRNA, FeatureType.retrotransposon):
                    format_hash['f_id'] = self.make_feature_link(f.id, f.id)
                if f.type in (FeatureType.CDS, FeatureType.tRNA, FeatureType.rRNA, FeatureType.ncRNA,
                              FeatureType.tmRNA, FeatureType.retrotransposon):
                    format_hash['strand'] = "+" if f.strand > 0 else "-"
                else:
                    format_hash['strand'] = ''
                format_hash['length'] = len(f) / 3 if f.type == FeatureType.CDS else len(f)
                match f.tmh_count(), f.signal_peptide, f.type:
                    case [_, 'LIPO', _]:
                        format_hash['destination'] = 'lipoprotein'
                    case [1, _, _]:
                        format_hash['destination'] = 'membrane anchor'
                    case [tmh, _, _] if tmh > 1:
                        format_hash['destination'] = 'membrane'
                    case [_, sp, _]:
                        format_hash['destination'] = 'envelope'
                    case [_, _, FeatureType.CDS]:
                        format_hash['destination'] = 'cytoplasm'
                    case [*_]:
                       format_hash['destination'] = ''
                format_hash['subsystem'] = ', '.join(f.subsystem)
                format_hash['has_cdd'] = 'Y' if len(f.cdd) else ''
                if len(f.blast):
                    format_hash['ident'] = f'{f.blast.hits[0].percent_id:.0f}'
                    format_hash['ci'] = self.get_color(f.blast.hits[0].percent_id)
                    format_hash['align'] = f'{f.blast.percent_aligned():0f}'
                    format_hash['ca'] = self.get_color(f.blast.percent_aligned())
                    format_hash['recall'] = f'{f.blast.percent_recall():0f}'
                    format_hash['cr'] = self.get_color(f.blast.percent_recall())
                else:
                    format_hash['ident'] = ''
                    format_hash['align'] = ''
                    format_hash['recall'] = ''
                    format_hash['ci'] = ''
                    format_hash['ca'] = ''
                    format_hash['cr'] = ''
                    format_hash['ct'] = self.get_color(20* len(set(f.taxon.split()) &
                                                               set(self.genome.properties['dominant taxon'].split())))
                table_body += ''''<tr>
                <td id=al>{f_id}</td> <td>{strand}</td> <td>{length}</td> <td>{destination}</td> <td>{subsystem}</td>
                <td>{has_cdd}</td> <td {ci}>{ident}</td> <td {ca}>{align}</td> <td {cr}>{recall}</td> 
                <td id=al>{descr}</td>
                <td {ct}>{taxon}</td>
                </tr>'''.format(**format_hash)
        html = html.replace('TABLE_BODY', table_body)
        return html

    def write_html(self, file=None):
        if not file:
            file = Path(self.exec.html_dir, self.genome.id, "index_of_features.html")
        file.parent.mkdir(exist_ok=True, parents=True)
        with open(Path(file), 'w') as handle:
            handle.write(self.make_html())

    def _make_html_template(self) -> str:
        """Creates and returns the html base for injecting the content in."""
        return '''<!DOCTYPE html>
<html>
<head>
    <meta charset="utf-8">
    <title>GENOME_NAME - all features</title>\n')
</head>
<body>

<script src="https://code.jquery.com/jquery-3.6.0.min.js" integrity="sha256-/xUj+3OJU5yExlq6GSYGSHk7tPXikynS7ogEvDej/m4=" crossorigin="anonymous"></script>
<link rel="stylesheet" type="text/css" href="https://cdn.datatables.net/1.11.5/css/jquery.dataTables.css">
<script type="text/javascript" charset="utf8" src="https://cdn.datatables.net/1.11.5/js/jquery.dataTables.js"></script>
<script type="text/javascript" charset="utf8">
$(document).ready( function () {
    $('#table_id').DataTable({
        "lengthMenu": [[20, 100, -1], [20, 100, "All"] ],
        "bSort" : false
        });
} );
</script>

<style>
  th {
    background-color: white;
      }
  #f {
    font-family: Calibri, sans-serif;
    text-align: center;
    padding: 0px;
    margin: 0px;
     }
  #cg {
    color: green;
     }
  #cr {
    color: red;
     }
  #cb {
    color: blue;
     }
  #co {
    color: orange;
     }
  #cw {
    color: white;
     }
  #al {
    text-align: left;
      }
</style>

<div id=f>
<table id="table_id" class="display">
    <thead>
        <tr>
TABLE_HEADERS
        </tr>
    </thead>
    <tbody>
TABLE_BODY
    </tbody>
</table> 
</div>
<div id=f>
<iframe src="" title="gene details" name="gene_details" style="border:none;width:100%;height:1000px;"></iframe>
</div>
</body>
</html>'''
