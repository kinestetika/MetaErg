import ast
import re
from pathlib import Path

from metaerg import context
from metaerg.datatypes import fasta

def write_html(dir):
    dir.mkdir(exist_ok=True, parents=True)
    file = Path(dir, 'index.html')
    with open(file, 'w') as handle:
        handle.write(make_html())


def make_html() -> str:
    """Injects the content into the html base, returns the html."""
    genome_list = []
    with open(context.GENOME_NAME_MAPPING_FILE) as handle:
        for line in handle:
            genome_list.append(line.split())

    checkm_results = {}
    checkm_result_file = Path(context.CHECKM_DIR, 'storage', 'bin_stats_ext.tsv')
    if checkm_result_file.exists():
        with open(checkm_result_file) as handle:
            for line in handle:
                words = line.split('\t')
                checkm_results[words[0]] = ast.literal_eval(words[1])

    gtdbtk_results = {}
    for file in (Path(context.GTDBTK_DIR, 'gtdbtk.ar53.summary.tsv'),
                 Path(context.GTDBTK_DIR, 'gtdbtk.bac120.summary.tsv')):
        if file.exists():
            with open(file) as handle:
                for line in handle:
                    if line.startswith("user_genome\t"):
                        continue
                    words = line.split("\t")
                    words[1] = re.sub('[a-z]__', ' ', words[1])
                    gtdbtk_results[words[0]] = words[1]

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
    for new_name, old_name, contig_file in genome_list:
        with fasta.FastaParser(contig_file) as fasta_reader:
            contigs = [c for c in fasta_reader]
        contigs.sort(key=lambda c: len(c['seq']), reverse=True)
        genome_size = sum(len(c['seq']) for c in contigs) / 1e6
        N50 = 0
        for c in contigs:
            N50 += len(c)
            if N50 >= genome_size / 2:
                N50 = len(c)
                break
        completeness = 0
        contamination = 0
        completeness_as_str = ''
        contamination_as_str = ''
        code = ''
        gtdbtk_classification = ''
        try:
            checkm_result = checkm_results[old_name]
            if not checkm_result:
                checkm_result = checkm_results[new_name]
            if checkm_result:
                completeness = float(checkm_result["Completeness"])
                contamination = float(checkm_result["Contamination"])
                completeness_as_str = f'{completeness:.1f}'
                contamination_as_str = f'{contamination:.1f}'
                code = int(checkm_result["Translation table"])
        except KeyError:
            pass
        try:
            gtdbtk_classification = gtdbtk_results[old_name]
            if not gtdbtk_classification:
                gtdbtk_classification = gtdbtk_results[new_name]
        except KeyError:
            pass
        link_target = Path(new_name, 'index.html')
        new_name = f'<a href="{link_target}">{new_name}</a>'
        genomes_html += '<tr>\n<td>{}</td><td>{}</td><td>{:.1f}</td><td>{}</td><td>{}</td><td {}>{}</td>' \
                        '<td {}>{}</td><td>{}</td>\n</tr>'.format(old_name,
                                                                       new_name,
                                                                       genome_size,
                                                                       N50,
                                                                       code,
                                                                       colors[int(max((completeness-76)/5,0))],
                                                                       completeness_as_str,
                                                                       colors[int(max((24-contamination)/5, 0))],
                                                                       contamination_as_str,
                                                                       gtdbtk_classification)
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
