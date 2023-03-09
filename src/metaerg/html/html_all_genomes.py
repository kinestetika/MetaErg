from pathlib import Path

from metaerg.datatypes import sqlite

def write_html(db_connection, dir):
    dir.mkdir(exist_ok=True, parents=True)
    file = Path(dir, 'index.html')
    with open(file, 'w') as handle:
        handle.write(make_html(db_connection))


def make_html(db_connection) -> str:
    """Injects the content into the html base, returns the html."""
    html = _make_html_template()
    left_aligned = 'file name classification'
    headers = ''
    for column in 'file name (Mb) N50 code completeness contamination classification'.split():
        if column in left_aligned:
            headers += f'<th id=al>{column}</th>\n'
        else:
            headers += f'<th>{column}</th>\n'
    html = html.replace('HEADERS', headers)

    colors = 'id=cr id=cr id=co id=cb id=cg'.split()
    genomes_html = ''
    for genome in sqlite.read_all_genomes(db_connection):
        name_html = f'<a href="{Path(genome.name, "index.html")}">{genome.name}</a>'
        genomes_html += '<tr>\n<td>{}</td><td>{}</td><td>{:.1f}</td><td>{:,}</td><td>{}</td><td {}>{:.1f}</td>' \
                        '<td {}>{:.1f}</td><td>{}</td>\n</tr>'.format(genome.input_fasta_file,
                                                                       name_html,
                                                                       genome.size / 1000000,
                                                                       genome.n50_contig_length,
                                                                       genome.genetic_code,
                                                                       colors[int(max((genome.fraction_complete*100-76)/5,0))],
                                                                       genome.fraction_complete * 100,
                                                                       colors[int(max((24-genome.fraction_contaminated*100)/5, 0))],
                                                                       genome.fraction_contaminated * 100,
                                                                       genome.top_taxon)
    html = html.replace('GENOMES', genomes_html)
    return html


def _make_html_template() -> str:
    """Creates and returns the html base for injecting the content in."""
    return '''<!DOCTYPE html>
<html>
<head>
    <meta charset="utf-8">
    <title>all annotated genomes</title>
</head>
<body>

<script src="https://code.jquery.com/jquery-3.6.0.min.js" integrity="sha256-/xUj+3OJU5yExlq6GSYGSHk7tPXikynS7ogEvDej/m4=" crossorigin="anonymous"></script>

<link rel="stylesheet" type="text/css" href="https://cdn.datatables.net/1.11.5/css/jquery.dataTables.css">

<script type="text/javascript" charset="utf8" src="https://cdn.datatables.net/1.11.5/js/jquery.dataTables.js"></script>

<script type="text/javascript" charset="utf8">
$(document).ready( function () {
    $('#table_id').DataTable({
        "lengthMenu": [[20, 100, -1], [20, 100, "All"] ],
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
HEADERS
      </tr><
</thead>
<tbody>
GENOMES
</tbody>
</table> 
</div>
<div id=f>
<iframe src="" title="gene details" name="gene_details" style="border:none;width:100%;height:1000px;"></iframe>
</div>
</body>
</html>'''
