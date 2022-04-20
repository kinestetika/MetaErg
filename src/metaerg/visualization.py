import os
import re
import pandas as pd
from pathlib import Path
from Bio import SeqIO

from metaerg import databases
from metaerg import predict
from metaerg import utils

PRODUCT_RE = re.compile('\[(\d+)/(\d+)\w\w@([\d,.]+)%\] \[(\d+)/(\d+)\] (.+)')
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
            if not " " in t:
                taxon = f'[{t}]'
                break
    return taxon


def make_feature_short_description(feature):
    return (f'{make_feature_product(feature)} {make_feature_taxon(feature)}')


def html_make_link(feature_id, description):
    return '<a target="_blank" href="{}.html">{}</a>'.format(feature_id, description)


def html_write_genome_stats(contig_dict):
    genome_stats = {}
    genome_stats["#contigs"] = len(contig_dict)
    total_size = 0
    for contig in contig_dict.values():
        total_size += len(contig)
    genome_stats["size"] = total_size
    cum_size = 0
    for contig in contig_dict.values():
        cum_size += len(contig)
        if cum_size > total_size/2:
            genome_stats["N50"] = len(contig)
            break
    genome_stats['#CDS'] = 0
    genome_stats['% coding'] = 0
    genome_stats['#rRNA'] = 0
    genome_stats['#tRNA'] = 0
    genome_stats['#ncRNA'] = 0
    genome_stats['#repeats'] = 0
    genome_stats['#retrotransposons'] = 0
    genome_stats['#CRISPR repeats'] = 0
    genome_stats['% repeats'] = 0
    genome_stats['total # features'] = 0
    taxon_dict = {}
    for contig in contig_dict.values():
        for feature in contig.features:
            if 'CDS' == feature.type:
                genome_stats['#CDS'] += 1
                genome_stats['% coding'] += feature.location.end - feature.location.start
            elif 'rRNA' == feature.type:
                genome_stats['#rRNA'] += 1
            elif 'tRNA' == feature.type:
                genome_stats['#tRNA'] += 1
            elif 'ncRNA' == feature.type:
                genome_stats['#ncRNA'] += 1
            elif 'repeat_region' == feature.type:
                genome_stats['#repeats'] += 1
                genome_stats['% repeats'] += feature.location.end - feature.location.start
            elif 'retrotransposon' == feature.type:
                genome_stats['#retrotransposons'] += 1
                genome_stats['% repeats'] += feature.location.end - feature.location.start
            elif 'crispr_repeat' == feature.type:
                genome_stats['#CRISPR repeats'] += 1
                genome_stats['% repeats'] += feature.location.end - feature.location.start
            t = utils.get_feature_qualifier(feature, 'taxonomy')
            if len(t):
                if t in taxon_dict.keys():
                    taxon_dict[t] += 1
                else:
                    taxon_dict[t] = 1
            genome_stats['total # features'] += 1

    genome_stats['% coding'] = f'{genome_stats["% coding"]/total_size*100:.2f}%'
    genome_stats['% repeats'] = f'{genome_stats["% repeats"]/total_size*100:.2f}%'
    genome_stats['dominant taxon'] = max(taxon_dict, key=taxon_dict.get)
    genome_stats['total # features'] = '<a target="_blank" href="index_of_features.html">{}</a>'.format(genome_stats['total # features'])

    # create html
    genome_stats_for_viz = []
    for (key, value) in genome_stats.items():
        genome_stats_for_viz.append({'property': key, 'value': value})
        print(f'{key:20}: {value}')
    df = pd.DataFrame(genome_stats_for_viz, columns=['property', 'value'])
    s = df.style.format(precision=1)
    s.set_table_styles([{'selector': 'td', 'props': 'font-family: Calibri, sans-serif;'}], overwrite=False)
    s.set_table_styles([{'selector': 'th.col_heading', 'props': 'font-family: Calibri, sans-serif;'}], overwrite=False)
    s.set_table_styles([{'selector': 'td', 'props': 'padding-left: 10px;'}], overwrite=False)
    s.set_sticky(axis=1)
    s.hide(axis="index")
    s.hide(axis="columns")
    s.to_html(f'index.html', doctype_html=True)
    return genome_stats


def html_create_blast_table_for_feature_page(feature, blast_results, is_cdd, dom_taxon=None, max_hits=0):
    blast_result = None
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
    seq = ''
    if 'CDS' == feature.type:
        seq = utils.pad_seq(feature.extract(contig)).translate(table=predict.TRANSLATION_TABLE)[:-1].seq
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

        writer.write('</body>\n</html>\n')


def html_write_feature_overview(writer, mag_name, contig_dict, genome_stats):
    writer.write('''<!doctype html>
<html>
<head>
    <meta charset="utf-8">''')
    writer.write(f'<title>{mag_name} - all features</title>\n')
    writer.write('''</head>
<body>

<script src="https://code.jquery.com/jquery-3.6.0.min.js" integrity="sha256-/xUj+3OJU5yExlq6GSYGSHk7tPXikynS7ogEvDej/m4=" crossorigin="anonymous"></script>

<link rel="stylesheet" type="text/css" href="https://cdn.datatables.net/1.11.5/css/jquery.dataTables.css">
  
<script type="text/javascript" charset="utf8" src="https://cdn.datatables.net/1.11.5/js/jquery.dataTables.js"></script>


<script type="text/javascript" charset="utf8">
$(document).ready( function () {
    $('#table_id').DataTable({
        "paging": false,
        "bSort" : false
        });
} );
</script>

<style>
  th {
	position: sticky;
	position: -webkit-sticky;
	top: 0px;
	z-index: 2;
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
    for column in 'id strand length type location CDD|ident|align|recall description taxon'.split():
        if column in left_aligned:
            writer.write(f'          <th id=al>{column}</th>\n')
        else:
            writer.write(f'          <th>{column}</th>\n')
    writer.write('''        </tr>
    </thead>
    <tbody>''')
    for contig in contig_dict.values():
        for feature in contig.features:
            product = make_feature_product(feature)
            feature_id = utils.get_feature_qualifier(feature, "id")
            description = ''
            signal_peptide = utils.get_feature_qualifier(feature, 'signal_peptide')
            writer.write('      <tr>\n')
            # id
            writer.write(f'          <td id=al>{feature_id}</td>\n')
            # strand
            if feature.type in ['CDS', 'tRNA', 'rRNA', 'ncRNA']:
                if feature.location.strand > 0:
                    writer.write('          <td>+</td>\n')
                else:
                    writer.write('          <td>-</td>\n')
            else:
                writer.write('          <td></td>\n')
            # length
            if feature.type == 'CDS':
                writer.write(f'          <td>{int(len(feature.location)/3-1)}</td>\n')
            else:
                writer.write(f'          <td>{len(feature.location)}</td>\n')
            # type:
            writer.write(f'          <td>{feature.type}</td>\n')
            # location
            try:
                number_of_tmh = int(utils.get_feature_qualifier(feature, 'transmembrane_helixes'))
            except ValueError:
                number_of_tmh = 0
            if number_of_tmh > 1:
                writer.write('          <td>membrane</td>\n')
            elif number_of_tmh == 1:
                writer.write('          <td>membrane anchor</td>\n')
            elif 'LIPO' in signal_peptide:
                writer.write('          <td>(lipoprotein)</td>\n')
            elif signal_peptide:
                writer.write('          <td>envelope</td>\n')
            elif 'CDS' == feature.type and product:
                writer.write('          <td>cytoplasm</td>\n')
            else:
                writer.write('          <td></td>\n')
            # homology blast hits (cdd)
            writer.write('          <td>')
            if utils.get_feature_qualifier(feature, 'cdd'):
                writer.write(f'          <b>{CENTER_DOT} </b>')
            else:
                writer.write(f'          <b id=cw>{CENTER_DOT} </b>')
            # homology blast hits (metaerg database)
            match = re.match(PRODUCT_RE, product)
            if (match):
                blast_percent_id = float(match.group(3))
                percent_id_color = 'cg'
                if blast_percent_id < 30:
                    percent_id_color = 'cr'
                elif blast_percent_id < 50:
                    percent_id_color = 'co'
                writer.write(f'<b id={percent_id_color}>{CENTER_DOT} </b>')
                blast_aligned = int(match.group(1)) / int(match.group(2))
                alignment_color = 'cg'
                if blast_aligned < 0.6:
                    alignment_color = 'cr'
                elif blast_aligned < 0.8:
                    alignment_color = 'co'
                writer.write(f'<b id={alignment_color}>{CENTER_DOT} </b>')
                blast_hit_count = int(match.group(4)) / int(match.group(5))
                hit_count_color = 'cg'
                if blast_hit_count < 0.5:
                    hit_count_color = 'cr'
                elif blast_hit_count < 0.8:
                    hit_count_color = 'co'
                writer.write(f'<b id={hit_count_color}>{CENTER_DOT} </b>')
                description = '<a target="_blank" href="{}.html">{}</a>'.format(feature_id, match.group(6))
            writer.write('</td>\n')
            if description:
                writer.write(f'          <td id=al>{description}</td>\n')
            else:
                writer.write(f'          <td id=al>{product}</td>\n')
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
                writer.write(f'          <td id={taxon_color}>{taxon}</td>\n')
            else:
                writer.write('          <td></td>\n')
            writer.write('      </tr>\n')

    writer.write('''    </tbody>
</table> 
</div>
</body>

</html>''')


def html_save_all(mag_name, contig_dict, blast_results):
    genome_stats = html_write_genome_stats(contig_dict)
    with open('index_of_features.html', 'w') as html_writer:
        html_write_feature_overview(html_writer, mag_name, contig_dict, genome_stats)

    #mag_dir = Path(mag_name)
    #os.chdir(mag_dir)
    # write html
    for contig in contig_dict.values():
        for feature in contig.features:
            html_write_page_for_feature(feature, contig, blast_results, genome_stats)