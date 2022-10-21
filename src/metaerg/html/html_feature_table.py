from pathlib import Path
import pandas as pd
from metaerg import context
from metaerg.datatypes.blast import taxon_at_genus, DBentry, BlastHit, BlastResult


@context.register_html_writer
def write_html(genome_name, feature_data: pd.DataFrame, genome_properties:dict, dir):
    dir.mkdir(exist_ok=True, parents=True)
    file = Path(dir, genome_name, "index_of_features.html")
    file.parent.mkdir(exist_ok=True, parents=True)
    with open(Path(file), 'w') as handle:
        handle.write(make_html(genome_name, feature_data, genome_properties))


def get_empty_format_dict():
    return {'f_id': '',
            'strand': '',
            'length': 0,
            'type': '',
            'destination': '',
            'subsystem': '',
            'has_cdd': '',
            'ident': '',
            'align': '',
            'recall': '',
            'description': '',
            'taxon': '',
            'ci': '', 'ca': '', 'cr': '', 'ct': ''}


def format_feature(f, format_hash, dominant_taxon, colors, path_to_feature_html):
    format_hash['f_id'] = f.id
    format_hash['taxon'] = taxon_at_genus(f.taxon)
    format_hash['type'] = f.type
    if f.type in ('CDS', 'rRNA', 'ncRNA', 'retrotransposon'):
        format_hash['description'] = '<a target="gene details" href="{}">{}</a>'.format(
            Path(path_to_feature_html, 'features', f'{f.id}.html'), f.descr)
    else:
        format_hash['description'] = f.descr
    if f.type in ('CDS', 'tRNA', 'rRNA', 'ncRNA', 'tmRNA', 'retrotransposon'):
        format_hash['strand'] = "+" if f.strand > 0 else "-"
    else:
        format_hash['strand'] = ''
    format_hash['length'] = len(f.seq)
    match f.tmh, f.signal_peptide, f.type:
        case [_, 'LIPO', _]:
            format_hash['destination'] = 'lipoprotein'
        case [1, _, _]:
            format_hash['destination'] = 'membrane anchor'
        case [tmh, _, _] if tmh > 1:
            format_hash['destination'] = 'membrane'
        case [_, sp, _] if len(sp):
            format_hash['destination'] = 'envelope'
        case [_, _, 'CDS']:
            format_hash['destination'] = 'cytoplasm'
        case [*_]:
            format_hash['destination'] = ''
    format_hash['subsystem'] = f.subsystems
    format_hash['has_cdd'] = 'Y' if f.cdd is not None else ''
    try:
        if len(f.blast):
            blast_result = eval(f.blast)
            format_hash['ident'] = f'{blast_result.hits[0].percent_id:.1f}'
            format_hash['ci'] = colors[min(int(blast_result.hits[0].percent_id / 20), len(colors) - 1)]
            format_hash['align'] = f'{blast_result.percent_aligned():.1f}'
            format_hash['ca'] = colors[min(int(blast_result.percent_aligned() / 20), len(colors) - 1)]
            format_hash['recall'] = f'{blast_result.percent_recall():.1f}'
            format_hash['cr'] = colors[min(int(blast_result.percent_recall() / 20), len(colors) - 1)]
    except SyntaxError:
        print(f.blast)
    dominant_taxon = dominant_taxon.split()
    taxon = f.taxon.split()
    format_hash['ct'] = colors[int(len(colors) * len(set(taxon) & set(dominant_taxon)) / (len(taxon) + 1))]


def format_hash_to_html(format_hash):
    return '''<tr>
    <td id=al>{f_id}</td> <td>{strand}</td> <td>{length:,}</td> <td>{type}</td> <td>{destination}</td> 
    <td>{subsystem}</td> <td>{has_cdd}</td> <td {ci}>{ident}</td> <td {ca}>{align}</td> <td {cr}>{recall}</td> 
    <td id=al>{description}</td>
    <td {ct}>{taxon}</td>
    </tr>'''.format(**format_hash)


def make_html(genome_name, feature_data: pd.DataFrame, genome_properties:dict, path_to_feature_html='') -> str:
    """Injects the content into the html base, returns the html."""
    html = _make_html_template()
    html = html.replace('GENOME_NAME', genome_name)
    colors = 'id=cr id=cr id=co id=cb id=cg'.split()

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
    previous_repeats = []
    prev_f = None
    for f in feature_data.itertuples():
        if f.type in ('crispr_repeat', 'repeat') and (not len(previous_repeats) or (f and f.type is prev_f.type)):
            previous_repeats.append(f)
        elif len(previous_repeats):
            format_hash = get_empty_format_dict()
            format_hash['type'] = f'[{len(previous_repeats)} {prev_f.type}s]' if len(previous_repeats) > 1 \
                                  else prev_f.type
            format_hash['length'] = previous_repeats[-1].end - previous_repeats[0].start
            previous_repeats.clear()
            table_body += format_hash_to_html(format_hash)
            if f.type in ('crispr_repeat', 'repeat'):
                previous_repeats.append(f)
            else:
                format_hash = get_empty_format_dict()
                format_feature(f, format_hash, genome_properties['dominant taxon'], colors, path_to_feature_html)
                table_body += format_hash_to_html(format_hash)
        else:
            format_hash = get_empty_format_dict()
            if type(f.seq) == float:
                print(f)
            format_feature(f, format_hash, genome_properties['dominant taxon'], colors, path_to_feature_html)
            table_body += format_hash_to_html(format_hash)
        prev_f = f
    html = html.replace('TABLE_BODY', table_body)
    return html


def _make_html_template() -> str:
    """Creates and returns the html base for injecting the content in."""
    return '''<!DOCTYPE html>
<html>
<head>
    <meta charset="utf-8">
    <title>GENOME_NAME - all features</title>
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
