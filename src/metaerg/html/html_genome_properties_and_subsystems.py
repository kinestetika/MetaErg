from pathlib import Path

from metaerg import context
from metaerg.datatypes.genome_properties import GENOME_PROPERTY_FORMATS

@context.register_html_writer
def write_html(genome_name, db_connection, genome_properties:dict, dir):
    dir.mkdir(exist_ok=True, parents=True)
    file = Path(dir, genome_name, 'index.html')
    with open(file, 'w') as handle:
        handle.write(make_html(genome_name, genome_properties))


def make_html(genome_name, genome_properties:dict) -> str:
    """injects the content into the html base, returns the html"""
    html = _make_html_template()
    html = html.replace('GENOME_NAME', genome_name)
    # genome properties
    html = html.replace('CONTENT_PROPERTIES', ''.join((f'<tr><td>{k}</td><td>{genome_properties[k]:{GENOME_PROPERTY_FORMATS[k]}}</td></tr>\n'
                                                       for k in GENOME_PROPERTY_FORMATS.keys())))
    # subsystem_summary
    sssummary = genome_properties['subsystem_summary']
    html = html.replace('CONTENT_SUBSYSTEM_SUMMARY', ''.join((f'<tr><td>{k}</td><td><a href="#{k}">{v:{"," if isinstance(v, int) else ".1%"}}</a></td></tr>\n'
                                                             for k, v in sssummary.items() if (isinstance(v, float) and v>=0.5) or isinstance(v, int) and v>0)))
    # subsystem_data
    subsystem_html = ''
    subsystem_data = genome_properties['subsystems']
    for subsystem, subsystem_genes in subsystem_data.items():
        subsystem_html += f'<a id="{subsystem}"><b id=f>{subsystem}</b></a>'
        if 'Secondary-metabolites' == subsystem:
            if len(subsystem_genes):
                subsystem_html += '<p id=f><a href="antismash/index.html" target="">View antismash results.</a></p>\n'
        else:
            subsystem_html += '<table id=f>\n'
            for gene_name, feature_ids in subsystem_genes.items():
                subsystem_html += f'<tr><td>{gene_name}</td><td>\n'
                subsystem_html += ', '.join('<a target="Gene Details" href="feature-details.html#{}">{}</a>'
                                            .format(f_id, f_id) for f_id in feature_ids)
                subsystem_html += f'</td></tr>\n'
            subsystem_html += '</table>\n'
    html = html.replace('CONTENT_SUBSYSTEMS', subsystem_html)
    return html


def _make_html_template() -> str:
    """should return the html base for injecting the content in. Returns the html"""
    return '''<!DOCTYPE html>
<html>
<head>
    <meta charset="utf-8">
        <title>GENOME_NAME - properties and subsystems</title>
</head>
<body>

<style>
  th {
    background-color: white;
     }
  #f {
    font-family: Calibri, sans-serif;
    text-align: left;
    padding: 0px;
    margin: 0px;
     }
  .accordion {
    background-color: #eee;
    color: #444;
    cursor: pointer;
    padding: 10px;
    width: 100%;
    border: none;
    text-align: left;
    outline: none;
    font-size: 15px;
    transition: 0.4s;
  }
  .active, .accordion:hover {
    background-color: #ccc; 
  }
  .panel {
    font-family: Calibri, sans-serif;
    padding: 0 18px;
    display: none;
    background-color: white;
    overflow: hidden;
  }
</style>

<h4 id=f><a href="feature_table.html">(View table with all genes here)</a></h4>
<br>
<h3 id=f>Genome properties</h4>
<table id=f>
    <thead>
        <tr>
          <th>property</th>
          <th>value</th>
        </tr>
    </thead>
    <tbody>
CONTENT_PROPERTIES
CONTENT_SUBSYSTEM_SUMMARY
    </tbody></table>
    <p></p>
    <p></p>
    <h3 id=f>Subsystems in Detail</h3>
CONTENT_SUBSYSTEMS
</body></html>'''
