import os
import re
import pandas as pd
import ast

from pathlib import Path

import predict
import data_model
from metaerg import databases
from metaerg import utils
from metaerg import subsystems

from Bio import SeqIO


PRODUCT_RE = re.compile(r'\[(\d+)/(\d+)\w\w@([\d,.]+)%\] \[(\d+)/(\d+)\] (.+)')
CENTER_DOT = u'\u25CF'


def make_feature_product(feature):
    product = utils.get_feature_qualifier(feature, "product")
    if not product:
        product = utils.get_feature_qualifier(feature, "profile")
    return product


def make_feature_taxon(feature):
    taxon = utils.get_feature_qualifier(feature, 'taxonomy')
    if "~" in taxon:
        for t in reversed(taxon.split("~ ")):
            if " " not in t:
                taxon = f'[{t}]'
                break
    return taxon


def make_feature_short_description(feature):
    return f'{make_feature_product(feature)} {make_feature_taxon(feature)}'


def html_make_link(feature_id, description):
    return '<a target="_blank" href="features/{}.html">{}</a>'.format(feature_id, description)


def html_write_genome_stats_and_subsystems(writer, mag_name, genome_stats, subsystems_hash):
    writer.write_html('''<!doctype html>
    <html>
    <head>
        <meta charset="utf-8">''')
    writer.write_html(f'<title>{mag_name} - properties and subsystems</title>\n')
    writer.write_html('''</head>
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

    <h4 id=f>Genome properties</h4>
    <table id=f>
        <thead>
            <tr>
              <th>key</th>
              <th>value</th>
            </tr>
        </thead>
        <tbody>''')
    for (key, value) in genome_stats.items():
        if key == 'total # features':
            value = '<a href="index_of_features.html">{}</a>'.format(genome_stats['total # features'])
        writer.write_html(f'          <tr><td>{key}</td><td>{value}</td></tr>\n')
    writer.write_html('        </tbody></table>\n')
    writer.write_html('        <h4 id=f>Subsystems overview</h4>\n')

    for subsystem in subsystems_hash.keys():
        s = subsystems.get_subsystem_stats(subsystems_hash[subsystem])
        subsystem_txt = ')'
        if s[0]:
            subsystem_txt = f'/{s[0]}): {s[2]*100:.0f}%'
        writer.write_html(f'          <button class="accordion">{subsystem} ({s[1]}{subsystem_txt}</button>\n')
        writer.write_html('          <div class="panel">\n')
        if '[secondary-metabolites]' == subsystem:
            if len(subsystems_hash[subsystem]):
                writer.write_html('          <p><a href="antismash/index.html" target="">View antismash results.</a></p>\n')
        elif isinstance(subsystems_hash[subsystem], list):
            writer.write_html('          <p>\n')
            for feature_id in subsystems_hash[subsystem]:
                writer.write_html(f'{html_make_link(feature_id, feature_id)} ')
            writer.write_html('\n          </p>\n')
        else:
            writer.write_html('          <table>\n')
            for phrase in subsystems_hash[subsystem].keys():
                writer.write_html(f'<tr><td>{phrase}</td><td>{" ".join(html_make_link(g, g) for g in subsystems_hash[subsystem][phrase])}</td></tr>\n')
            writer.write_html('          </table>\n')
        writer.write_html('          </div>\n')
    writer.write_html('''<script>
    var acc = document.getElementsByClassName("accordion");
    var i;
    
    for (i = 0; i < acc.length; i++) {
      acc[i].addEventListener("click", function() {
        this.classList.toggle("active");
        var panel = this.nextElementSibling;
        if (panel.style.display === "block") {
          panel.style.display = "none";
        } else {
          panel.style.display = "block";
        }
      });
    }
    </script>
    </div></body></html>''')


def html_create_blast_table_for_feature_page(feature, blast_results, is_cdd, dom_taxon=None, max_hits=0):
    query_length = int(len(feature.location) / 3)
    feature_id = utils.get_feature_qualifier(feature, 'id')
    if is_cdd:
        try:
            blast_result = blast_results['cdd'][feature_id]
        except KeyError:
            return None
        columns = ['% id', 'query pos', 'hit pos', 'description']
    else:
        if 'CDS' == feature.type:
            try:
                blast_result = blast_results['diamond'][feature_id]
            except KeyError:
                return None
        else:
            try:
                blast_result = blast_results['blastn'][feature_id]
            except KeyError:
                return None
            query_length = len(feature.location)
        columns = ['% id', 'query pos', 'hit pos', 'description', 'taxon']

    blast_result_copy = []
    styles = []
    for hit in blast_result:
        style = {}
        blast_query_aligned = abs(hit["query_start"] - hit["query_end"]) / query_length
        if is_cdd:
            cdd_id = int(hit["hit_id"][4:])
            cdd_item = databases.CDD[cdd_id]
            blast_hit_aligned = abs(hit["hit_start"] - hit["hit_end"]) / cdd_item[3]
            table_row = {'% id': hit['percent_id'],
                         'query pos': f'{hit["query_start"]}-{hit["query_end"]}/{query_length}',
                         'hit pos': f'{hit["hit_start"]}-{hit["hit_end"]}/{cdd_item[3]}',
                         'description': f'{cdd_item[0]}|{cdd_item[1]} {cdd_item[2]}'}
        else:
            db_entry = databases.decipher_database_id(hit['hit_id'])
            blast_hit_aligned = abs(hit["hit_start"] - hit["hit_end"]) / db_entry["length"]
            table_row = {'% id': hit['percent_id'],
                         'query pos': f'{hit["query_start"]}-{hit["query_end"]}/{query_length}',
                         'hit pos': f'{hit["hit_start"]}-{hit["hit_end"]}/{db_entry["length"]}',
                         'description': db_entry["descr"],
                         'taxon': f'[{db_entry["taxon"].replace("~", ";")}]'}
        if table_row['% id'] < 30:
            style['% id'] = 'color: red; text-align: center;'
        elif table_row['% id'] < 50:
            style['% id'] = 'color: orange; text-align: center;'
        else:
            style['% id'] = 'color: green; text-align: center;'
        if 'taxon' in table_row.keys():
            taxon_list = table_row['taxon'].split("; ")
            score = 0
            for j in range(min(len(taxon_list), len(dom_taxon))):
                if taxon_list[j] == dom_taxon[j]:
                    score += 1
            if score < 2:
                style['taxon'] = 'color: red;'
            elif score < 4:
                style['taxon'] = 'color: orange;'
            elif score < 5:
                style['taxon'] = 'color: blue;'
            else:
                style['taxon'] = 'color: green;'
        if blast_hit_aligned < 0.5:
            style['hit pos'] = 'color: red; text-align: center;'
        elif blast_hit_aligned < 0.8:
            style['hit pos'] = 'color: orange; text-align: center;'
        else:
            style['hit pos'] = 'color: green; text-align: center;'
        if blast_query_aligned < 0.5:
            style['query pos'] = 'color: red; text-align: center;'
        elif blast_query_aligned < 0.8:
            style['query pos'] = 'color: orange; text-align: center;'
        else:
            style['query pos'] = 'color: green; text-align: center;'

        blast_result_copy.append(table_row)
        styles.append(style)
        if max_hits and len(blast_result_copy) >= max_hits:
            break
    df = pd.DataFrame(blast_result_copy, columns=columns)
    table_style = pd.DataFrame(styles, columns=columns)
    s = df.style.format(precision=1)
    s.set_table_styles({'description': [{'selector': 'th.col_heading', 'props': 'text-align: left;'}],
                        'taxon': [{'selector': 'th.col_heading', 'props': 'text-align: left;'}]})
    s.set_table_styles([{'selector': 'th.col_heading', 'props': 'text-align: center;'}], overwrite=False)
    s.set_table_styles([{'selector': 'td', 'props': 'font-family: Calibri, sans-serif;'}], overwrite=False)
    s.set_table_styles([{'selector': 'th.col_heading', 'props': 'font-family: Calibri, sans-serif;'}], overwrite=False)
    s.set_table_styles([{'selector': 'th.col_heading', 'props': 'border-bottom: 1px solid black;'}], overwrite=False)
    s.set_table_styles([{'selector': 'td', 'props': 'padding-left: 10px;'}], overwrite=False)
    s.apply(lambda df: table_style, axis=None)
    s.set_sticky(axis=1)
    s.hide(axis="index")
    return s


def html_write_page_for_feature(feature, contig, blast_results, genome_stats):
    feature_id = utils.get_feature_qualifier(feature, 'id')
    header = f'>{feature_id} {make_feature_short_description(feature)}'
    if 'CDS' == feature.type:
        seq = data_model.pad_seq(feature.extract(contig)).translate(table=utils.TRANSLATION_TABLE)[:-1].seq
    else:
        seq = feature.extract(contig).seq

    with open(f'{feature_id}.html', 'w') as writer:
        writer.write('<!DOCTYPE html>\n<html>\n<body>\n<meta charset="utf-8">\n\n')
        writer.write(f'<h2>{header}</h2>')
        writer.write('<style>\ndiv {\n  font-family: Calibri, sans-serif;\n  word-wrap: break-word;\n}\n</style>\n')
        writer.write(f'<div>{seq}</div>')

        writer.write('<style>\nh3 {\n  font-family: Calibri, sans-serif;\n  margin: 0px;\n}\n</style>\n')

        database_blast_table = html_create_blast_table_for_feature_page(feature, blast_results, is_cdd=False, max_hits=10,
                                                                        dom_taxon=genome_stats['dominant taxon'].split("~ "))
        if database_blast_table:
            writer.write(f'<hr><h3>database hits</h3>')
            writer.write(database_blast_table.to_html(None))

        cdd_blast_table = html_create_blast_table_for_feature_page(feature, blast_results, is_cdd=True, max_hits=10,
                                                                   dom_taxon=genome_stats['dominant taxon'].split("~ "))
        if cdd_blast_table:
            writer.write(f'<hr><h3>conserved domain database hits</h3>')
            writer.write(cdd_blast_table.to_html(None))

        writer.write('<table>\n')
        for key in feature.qualifiers:
            writer.write('  <tr>\n')
            writer.write(f'    <td>{key}</td><td>{utils.get_feature_qualifier(feature, key)}</td>\n')
            writer.write('  </tr>\n')
        writer.write('</table></body>\n</html>\n')


# To make the header sticky, add the following lines to th style:
# position: sticky;
# position: -webkit - sticky;
# top: 0
# px;
# z - index: 2;


def html_write_feature_overview(writer, mag_name, contig_dict, genome_stats, blast_results):
    writer.write_html('''<!doctype html>
<html>
<head>
    <meta charset="utf-8">''')
    writer.write_html(f'<title>{mag_name} - all features</title>\n')
    writer.write_html('''</head>
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
''')
    left_aligned = 'id description'
    for column in 'id strand length type location subsystem CDD ident align recall description taxon'.split():
        if column in left_aligned:
            writer.write_html(f'          <th id=al>{column}</th>\n')
        else:
            writer.write_html(f'          <th>{column}</th>\n')
    writer.write_html('''        </tr>
    </thead>
    <tbody>''')
    for contig in contig_dict.values():
        for feature in contig.features:
            product = make_feature_product(feature)
            feature_id = utils.get_feature_qualifier(feature, "id")
            description = ''
            signal_peptide = utils.get_feature_qualifier(feature, 'signal_peptide')
            writer.write_html('      <tr>\n')
            # id
            writer.write_html(f'          <td id=al>{feature_id}</td>\n')
            # strand
            if feature.type in ['CDS', 'tRNA', 'rRNA', 'ncRNA']:
                if feature.location.strand > 0:
                    writer.write_html('          <td>+</td>\n')
                else:
                    writer.write_html('          <td>-</td>\n')
            else:
                writer.write_html('          <td></td>\n')
            # length
            if feature.type == 'CDS':
                writer.write_html(f'          <td>{int(len(feature.location) / 3 - 1)}</td>\n')
            else:
                writer.write_html(f'          <td>{len(feature.location)}</td>\n')
            # type:
            writer.write_html(f'          <td>{feature.type}</td>\n')
            # location
            try:
                number_of_tmh = int(utils.get_feature_qualifier(feature, 'transmembrane_helixes'))
            except ValueError:
                number_of_tmh = 0
            if number_of_tmh > 1:
                writer.write_html('          <td>membrane</td>\n')
            elif number_of_tmh == 1:
                writer.write_html('          <td>membrane anchor</td>\n')
            elif 'LIPO' in signal_peptide:
                writer.write_html('          <td>(lipoprotein)</td>\n')
            elif signal_peptide:
                writer.write_html('          <td>envelope</td>\n')
            elif 'CDS' == feature.type and product:
                writer.write_html('          <td>cytoplasm</td>\n')
            else:
                writer.write_html('          <td></td>\n')
            # subsystem
            subsystem = utils.get_feature_qualifier(feature, 'subsystem')
            writer.write_html(f'          <td>{subsystem}</td>\n')
            # homology blast hits (cdd)
            if utils.get_feature_qualifier(feature, 'cdd'):
                writer.write_html(f'          <td>Y</td>\n')
                # construct a description in case of no blast hit below
                cdd_result = blast_results['cdd'][feature_id]
                if len(cdd_result):
                    cdd_id = int(cdd_result[0]["hit_id"][4:])
                    cdd_item = databases.CDD[cdd_id]
                    txt = cdd_item[2]
                    if len(txt) > 30:
                        txt = txt[:30] + '...'
                    description = f'{cdd_item[0]}|{cdd_item[1]} {txt}'
            else:
                writer.write_html(f'          <td></td>\n')
            # homology blast hits (metaerg database)
            match = re.match(PRODUCT_RE, product)
            if match:
                blast_percent_id = float(match.group(3))
                color = 'cg'
                if blast_percent_id < 30:
                    color = 'cr'
                elif blast_percent_id < 50:
                    color = 'co'
                writer.write_html(f'<td id={color}>{blast_percent_id:.0f}</td>\n')
                blast_aligned = (int(match.group(1)) / int(match.group(2))) * 100
                color = 'cg'
                if blast_aligned < 60:
                    color = 'cr'
                elif blast_aligned < 80:
                    color = 'co'
                writer.write_html(f'<td id={color}>{blast_aligned:.0f} </td>\n')
                blast_hit_count = (int(match.group(4)) / int(match.group(5))) * 100
                color = 'cg'
                if blast_hit_count < 50:
                    color = 'cr'
                elif blast_hit_count < 80:
                    color = 'co'
                writer.write_html(f'<td id={color}>{blast_hit_count:.0f} </td>\n')
                description = match.group(6)
            else:
                writer.write_html('          <td></td><td></td><td></td>\n')

            if description:
                writer.write_html('<td id=al><a target="gene_details" href="features/{}.html">{}</a></td>\n'.format(feature_id, description))
            else:
                writer.write_html(f'          <td id=al>{product}</td>\n')
            # taxon
            taxonomy = utils.get_feature_qualifier(feature, 'taxonomy')
            if taxonomy:
                taxon = make_feature_taxon(feature)
                score = 0
                taxon_list = taxonomy.split("~ ")
                dominant_taxon = genome_stats['dominant taxon'].split("~ ")
                for j in range(min(len(taxon_list), len(dominant_taxon))):
                    if taxon_list[j] == dominant_taxon[j]:
                        score += 1
                taxon_color = 'cg'
                if score < 2:
                    taxon_color = 'cr'
                elif score < 4:
                    taxon_color = 'co'
                elif score < 5:
                    taxon_color = 'cb'
                writer.write_html(f'          <td id={taxon_color}>{taxon}</td>\n')
            else:
                writer.write_html('          <td></td>\n')
            writer.write_html('      </tr>\n')

    writer.write_html('''    </tbody>
</table> 
</div>
<div id=f>
<iframe src="" title="gene details" name="gene_details" style="border:none;width:100%;height:1000px;"></iframe>
</div>
</body>
</html>''')


def html_write_all(mag_name, genome_stats, contig_dict, blast_results, subsystem_hash):
    run_and_read.execution.log(f'({mag_name}) Writing html index and overview...')
    with open('index.html', 'w') as html_writer:
        html_write_genome_stats_and_subsystems(html_writer, mag_name, genome_stats, subsystem_hash)
    with open('index_of_features.html', 'w') as html_writer:
        html_write_feature_overview(html_writer, mag_name, contig_dict, genome_stats, blast_results)

    feature_dir = Path('features')
    feature_dir.mkdir()
    os.chdir(feature_dir)
    run_and_read.execution.log(f'({mag_name}) Writing html page for each gene...')
    for contig in contig_dict.values():
        for feature in contig.features:
            html_write_page_for_feature(feature, contig, blast_results, genome_stats)
    os.chdir('..')


def html_write_mag_table(writer, target_dir, checkm_dir='checkm', gtdbtk_dir='gtdbtk'):
    mag_list = []
    with open(Path(target_dir, 'mag.name.mapping.txt')) as handle:
        for line in handle:
            mag_list.append(line.split())

    checkm_results = {}
    checkm_result_file = Path(checkm_dir, 'storage', 'bin_stats_ext.tsv')
    #print(checkm_result_file.absolute())
    if checkm_result_file.exists():
        with open(checkm_result_file) as handle:
            for line in handle:
                words = line.split('\t')
                checkm_results[words[0]] = ast.literal_eval(words[1])
    #print(checkm_results)
    gtdbtk_results = {}
    for file in (Path(gtdbtk_dir, 'gtdbtk.ar53.summary.tsv'), Path(gtdbtk_dir, 'gtdbtk.bac120.summary.tsv')):
        if file.exists():
            with open(file) as handle:
                for line in handle:
                    if line.startswith("user_genome\t"):
                        continue
                    words = line.split("\t")
                    words[1] = re.sub('[a-z]__', ' ', words[1])
                    gtdbtk_results[words[0]] = words[1]

    writer.write_html('''<!doctype html>
    <html>
    <head>
        <meta charset="utf-8">''')
    writer.write_html(f'<title>all annotated genomes</title>\n')
    writer.write_html('''</head>
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
          <tr>''')
    left_aligned = 'file name classification'
    for column in 'file name (Mb) N50 code completeness contamination classification'.split():
        if column in left_aligned:
            writer.write_html(f'          <th id=al>{column}</th>\n')
        else:
            writer.write_html(f'          <th>{column}</th>\n')
    writer.write_html('''        </tr>
    </thead>
    <tbody>''')

    for new_name, old_name in mag_list:
        contig_dict = {}
        gbk_file = Path(target_dir, 'gbk', new_name)
        with open(gbk_file) as handle:
            for gb_record in SeqIO.parse(handle, "genbank"):
                contig_dict[gb_record.id] = gb_record
        genome_stats = predict.compile_genome_stats(new_name, contig_dict)
        writer.write_html('        <tr>')
        writer.write_html(f'            <td id=al>{old_name}</td>')
        writer.write_html(f'            <td id=al><a href="{new_name}/index.html">{new_name}</a></td>')
        writer.write_html(f'            <td>{genome_stats["size"] / 1e6:.2f}</td>')
        writer.write_html(f'            <td>{genome_stats["N50"]}</td>')
        completeness = ''
        contamination = ''
        code = ''
        gtdbtk_classification = ''
        try:
            checkm_result = checkm_results[Path(old_name).stem]
            if not checkm_result:
                checkm_result = checkm_results[new_name]
            if checkm_result:
                completeness = f'{float(checkm_result["Completeness"]):.1f}'
                contamination = f'{float(checkm_result["Contamination"]):.1f}'
                code = int(checkm_result["Translation table"])
        except KeyError:
            pass
        try:
            gtdbtk_classification = gtdbtk_results[Path(old_name).stem]
            if not gtdbtk_classification:
                gtdbtk_classification = gtdbtk_results[new_name]
        except KeyError:
            pass
        writer.write_html(f'            <td>{code}</td>')
        writer.write_html(f'            <td>{completeness}</td>')
        writer.write_html(f'            <td>{contamination}</td>')
        writer.write_html(f'            <td id=al>{gtdbtk_classification}</td>')
        writer.write_html('        </tr>')

    writer.write_html('''    </tbody>
</table> 
</div>
<div id=f>
<iframe src="" title="gene details" name="gene_details" style="border:none;width:100%;height:1000px;"></iframe>
</div>
</body>
</html>''')
