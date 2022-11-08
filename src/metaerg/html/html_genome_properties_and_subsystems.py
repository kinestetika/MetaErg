from pathlib import Path
import pandas as pd
from metaerg import context


GENOME_PROPERTY_FORMATS = {'size': ',',
                           '% GC': '.1%',
                           'N50': ',',
                           '# proteins': ',',
                           '% coding': '.1%',
                           'mean protein length (aa)': '.1f',
                           '# ribosomal RNA': ',',
                           '# transfer RNA': ',',
                           '# non-coding RNA': ',',
                           '# retrotransposons': ',',
                           '# CRISPR repeats': ',',
                           '# other repeats': ',',
                           '# total features': ',',
                           '% repeats': '.1%',
                           '% CDS classified to taxon': '.1%',
                           '% of CDS classified to dominant taxon': '.1%',
                           'dominant taxon': '<',
                           }

@context.register_html_writer
def write_html(genome_name, feature_data: pd.DataFrame, genome_properties:dict, dir):
    dir.mkdir(exist_ok=True, parents=True)
    file = Path(dir, genome_name, 'index.html')
    with open(file, 'w') as handle:
        handle.write(make_html(genome_name, feature_data, genome_properties))


def make_html(genome_name, feature_data: pd.DataFrame, genome_properties:dict, path_to_feature_html='') -> str:
    """injects the content into the html base, returns the html"""
    html = _make_html_template()
    html = html.replace('GENOME_NAME', genome_name)
    # genome properties
    html = html.replace('CONTENT_PROPERTIES', ''.join((f'<tr><td>{k}</td><td>{genome_properties[k]:{GENOME_PROPERTY_FORMATS[k]}}</td></tr>\n'
                                                       for k in GENOME_PROPERTY_FORMATS.keys())))
    subsystem_html = ''
    subsystem_data = genome_properties['subsystems']
    # subsystem_data
    for subsystem in subsystem_data.index.unique('subsystem'):
        sub_data = subsystem_data.loc[subsystem, :]
        subsystem_html += f'<b id=f>{subsystem}</b>'
        if 'secondary-metabolites' == subsystem:
            if len(sub_data.index):
                subsystem_html += '<p id=f><a href="antismash/index.html" target="">View antismash results.</a></p>\n'
        else:
            subsystem_html += '<table id=f>\n'
            for function in sub_data.itertuples():
                if function.Index:
                    subsystem_html += f'<tr><td>{function.Index}</td><td>\n'
                    subsystem_html += ', '.join('<a target="_blank" href="{}">{}</a>'
                                                .format(Path(path_to_feature_html, 'features', f'{f_id.split("@")[0]}.html'),
                                                        f_id) for f_id in function.genes.split())
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
    </tbody></table>
    <p></p>
    <h4 id=f><a href="index_of_features.html">View table with all genes</a></h4>
    <p></p>
    <h3 id=f>Subsystems overview</h3>
CONTENT_SUBSYSTEMS
</body></html>'''
