from pathlib import Path

from metaerg import context
from metaerg.datatypes.functional_genes import format_list_of_subsystem_genes
from metaerg.datatypes import sqlite
from metaerg.datatypes.blast import taxon_at_genus

CURRENT_REGION_ID = ''
CURRENT_REGION_COLOR_DONGLE = True

@context.register_html_writer
def write_html(genome, db_connection, dir):
    dir.mkdir(exist_ok=True, parents=True)
    file = Path(dir, genome.name, "feature_table.html")
    file.parent.mkdir(exist_ok=True, parents=True)
    with open(Path(file), 'w') as handle:
        handle.write(make_html(genome, db_connection))


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
            'is_part_of_cog': '',
            'ci': '', 'ca': '', 'cr': '', 'ct': '', 'cc': ''}


def format_feature(f: sqlite.Feature, format_hash, top_taxon, colors):
    global CURRENT_REGION_ID, CURRENT_REGION_COLOR_DONGLE
    if f.homologous_group_id > 0:  # these ids start with 1
        format_hash['is_part_of_cog'] = 'al_bold'
        paralogue_warning = '!' if f.homologous_group_feature_is_paralogue else ''
        expression_level = f'id =gl{str(f.homologous_group_codon_usage_bias_quantile)}'
        representation_level = f'id =gl{str(int(f.homologous_group_taxon_representation*10-1))}'
        if -1 < f.homologous_group_selective_pressure_quantile <= 1:
            selection_type = ' id=sp>\u25CF'
        elif f.homologous_group_selective_pressure_quantile >= 8:
            selection_type = ' id=sd>\u25CF'
        else:
            selection_type = '>'
        format_hash['f_id'] = f'{f.id}' \
                              f'<b id=symbols>' \
                                f'<b id=cg><b {representation_level}>\u25A0</b></b>' \
                                f'<b id=co>{paralogue_warning}</b>' \
                                f'<b{selection_type}</b>' \
                              f'</b>' \
                              f'<b id=cblack><b {expression_level}>\u25BC</b></b>'
    else:
        format_hash['is_part_of_cog'] = 'al'
        format_hash['f_id'] = f.id
    format_hash['taxon'] = taxon_at_genus(f.taxon)
    format_hash['type'] = f.type
    if f.type in ('CDS', 'rRNA', 'ncRNA', 'retrotransposon'):
        format_hash['description'] = '<a target="Gene Details" href="feature-details.html#{}">{}</a>'.format(
            f'{f.id}', f.descr)
    else:
        format_hash['description'] = f.descr
    if f.type in ('CDS', 'tRNA', 'rRNA', 'ncRNA', 'tmRNA', 'retrotransposon'):
        format_hash['strand'] = "+" if f.strand > 0 else "-"
    else:
        format_hash['strand'] = ''
    if 'CDS' == f.type:
        format_hash['length'] = f'{len(f.aa_seq):,} aa'
    else:
        format_hash['length'] = f'{f.end-f.start:,} nt'
    match f.tmh, f.signal_peptide, f.type:
        case [1, _, _]:
            format_hash['destination'] = 'membrane anchored'
        case [tmh, _, _] if tmh > 1:
            format_hash['destination'] = 'membrane'
        case [_, sp, _] if len(sp):
            format_hash['destination'] = f'envelope ({sp})'
        case [_, _, 'CDS']:
            format_hash['destination'] = 'cytoplasm'
        case [*_]:
            format_hash['destination'] = ''
    format_hash['subsystem'] = format_list_of_subsystem_genes(f.subsystems)
    format_hash['has_cdd'] = 'Y' if f.cdd is not None else ''
    if f.blast:
        format_hash['ident'] = f'{f.blast.hits[0].percent_id:.1f}'
        format_hash['ci'] = colors[min(int(f.blast.hits[0].percent_id / 20), len(colors) - 1)]
        format_hash['align'] = f'{f.blast.percent_aligned():.1f}'
        format_hash['ca'] = colors[min(int(f.blast.percent_aligned() / 20), len(colors) - 1)]
        format_hash['recall'] = f'{f.blast.percent_recall():.1f}'
        format_hash['cr'] = colors[min(int(f.blast.percent_recall() / 20), len(colors) - 1)]
    top_taxon = top_taxon.split()
    taxon = f.taxon.split()
    format_hash['ct'] = colors[int(len(colors) * len(set(taxon) & set(top_taxon)) / (len(taxon) + 1))]
    if f.parent:
        if f.parent != CURRENT_REGION_ID:
            CURRENT_REGION_COLOR_DONGLE = not CURRENT_REGION_COLOR_DONGLE
        CURRENT_REGION_ID = f.parent
        format_hash['cc'] = ' id=cc1' if CURRENT_REGION_COLOR_DONGLE else ' id=cc2'


def format_hash_to_html(format_hash):
    return '''<tr>
    <td id={is_part_of_cog}>{f_id}</td> <td{cc}>{strand}</td> <td>{length}</td> <td>{type}</td> <td>{destination}</td> 
    <td>{subsystem}</td> <td>{has_cdd}</td> <td{ci}>{ident}</td> <td{ca}>{align}</td> <td{cr}>{recall}</td> 
    <td id=al>{description}</td>
    <td{ct}>{taxon}</td>
    </tr>'''.format(**format_hash)


def make_html(genome, db_connection) -> str:
    """Injects the content into the html base, returns the html."""
    html = _make_html_template()
    html = html.replace('GENOME_NAME', genome.name)
    colors = [' id=cr', ' id=cr', ' id=co', ' id=cb', ' id=cg']

    if context.DO_CLUSTER_GENOMES:
        html = html.replace('LEGEND', '''<p id=f>Legend:
<b id=symbols><b id=cg>\u25A0</b></b> Fraction of taxa with this gene. 
<b id=symbols><b id=sp>\u25CF</b></b> Purifying selection.
<b id=symbols><b id=sd>\u25CF</b></b> Diversifying selection.
<b id=symbols><b id=co>!</b></b> This is a paralogue.
<b id=cblack>\u25BC</b> Estimated expression level.
</p>
''')
    else:
        html = html.replace('LEGEND','')

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
    repeat_types = set()
    for f in sqlite.read_all_features(db_connection):
        if f.type == 'region':
            continue
        if f.type == 'CRISPR':
            previous_repeats.append(f)
            repeat_types.add(f.type)
        elif f.type == 'repeat_unit':
            previous_repeats.append(f)
            repeat_types.add('repeat')
        elif f.type == 'binding_site':
            previous_repeats.append(f)
            repeat_types.add('spacer')
        else:
            if previous_repeats:
                format_hash = get_empty_format_dict()
                format_hash['description'] = previous_repeats[0].id + ' ... ' + previous_repeats[-1].id
                format_hash['type'] = f'{len(previous_repeats)} {", ".join(repeat_types)}(s)'
                format_hash['length'] = f'{previous_repeats[-1].end - previous_repeats[0].start:,} nt'
                if previous_repeats[-1].parent:
                    format_hash['cc'] = ' id=cc1' if CURRENT_REGION_COLOR_DONGLE else ' id=cc2'
                previous_repeats.clear()
                repeat_types.clear()
                table_body += format_hash_to_html(format_hash)
        if not f.type in ('repeat_region', 'spacer', 'binding_site', 'repeat_unit', 'CRISPR'):
            format_hash = get_empty_format_dict()
            format_feature(f, format_hash, genome.top_taxon, colors)
            table_body += format_hash_to_html(format_hash)
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

LEGEND

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
  #cblack {
    color: black;
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
  #al_bold {
    text-align:  left;
    font-weight: bold;
      }
  #cc1 {
    background-color:#CCDDEE;
      }
  #cc2 {
    background-color:#CCEEDD;
      }
  #symbols {
      font-size: 1.5em;
      }
  #gl-1 {
    opacity: 0.0;
      }
  #gl0 {
    opacity: 0.1;
      }
  #gl1 {
    opacity: 0.2;
      }
  #gl2 {
    opacity: 0.3;
      }
  #gl3 {
    opacity: 0.4;
      }
  #gl4 {
    opacity: 0.5;
      }
  #gl5 {
    opacity: 0.6;
      }
  #gl6 {
    opacity: 0.7;
      }
  #gl7 {
    opacity: 0.8;
      }
  #gl8 {
    opacity: 0.9;
      }
  #gl9 {
    opacity: 1.0;
      }
  #gl9 {
    opacity: 1.0;
      }
  #sp {
    color: darkgoldenrod
  }
  #sd {
    color: darkslateblue
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
