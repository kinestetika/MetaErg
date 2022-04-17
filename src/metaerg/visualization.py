import os
import re
import pandas as pd
from pathlib import Path
from Bio import SeqIO

from metaerg import databases
from metaerg import predict
from metaerg import utils

PRODUCT_RE = re.compile('\[(\d+)/(\d+)\w\w@([\d,.]+)%\] \[(\d+)/(\d+)\] (.+)')


def html_make_link(feature_id, description):
    return '<a target="_blank" href="{}.html">{}</a>'.format(feature_id, description)


def html_write_genome_stats(mag_name, contig_dict):
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
    taxon = utils.get_feature_qualifier(feature, 'taxonomy')
    feature_id = utils.get_feature_qualifier(feature, 'id')
    if "~" in taxon:
        for t in reversed(taxon.split("~ ")):
            if not " " in t:
                taxon = f'[{t}]'
                break
    header = (f'>{feature_id} {utils.get_feature_qualifier(feature, "product")} {taxon}')
    seq = ''
    if 'CDS' == feature.type:
        seq = utils.get_feature_qualifier(feature, 'translation')
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


def html_write_feature_overview(mag_name, contig_dict, genome_stats):
    columns='id strand length type location D % A R description taxon'.split()
    dominant_taxon = genome_stats['dominant taxon'].split("~ ")
    feature_list = []
    style_list = []
    for contig in contig_dict.values():
        for feature in contig.features:
            f = {'id': utils.get_feature_qualifier(feature, 'id'),
                 'strand': '',
                 'type': feature.type,
                 'location': '',
                 'D': '',
                 '%': '',
                 'A': '',
                 'R': '',
                 'description': '',
                 'taxon': ''}
            style = {'id': 'text-align: left;',
                     'strand': 'text-align: center;',
                     'length': 'text-align: center;',
                     'type': '' 'text-align: center;',
                     'location': 'text-align: center;',
                     'D':'',
                     '%': '',
                     'A': '',
                     'R': '',
                     'description': 'text-align: left;',
                     'taxon': 'text-align: left;'}

            # strand
            if feature.type in ['CDS', 'tRNA', 'rRNA', 'ncRNA']:
                if feature.location.strand > 0:
                    f['strand'] = '+'
                else:
                    f['strand'] = '-'

            # length
            if feature.type == 'CDS':
                f['length'] = len(utils.get_feature_qualifier(feature, 'translation'))
            else:
                f['length'] = len(feature.location)

            # taxonomy
            taxonomy = utils.get_feature_qualifier(feature, 'taxonomy')
            if taxonomy:
                taxon_list = taxonomy.split("~ ")
                for t in reversed(taxon_list):
                    if not " " in t:
                        f['taxon'] = f'[{t}]'
                        break
                score = 0
                for j in range(min(len(taxon_list), len(dominant_taxon))):
                    if taxon_list[j] == dominant_taxon[j]:
                        score += 1
                if score < 2:
                    style['taxon'] = 'color: red;'
                elif score < 4:
                    style['taxon'] = 'color: orange;'
                elif score < 5:
                    style['taxon'] = 'color: blue;'
                else:
                    style['taxon'] = 'color: green;'

            # blast hit, product
            f['description'] = utils.get_feature_qualifier(feature, 'product')
            match = re.match(PRODUCT_RE, utils.get_feature_qualifier(feature, 'product'))
            if (match):
                blast_aligned = int(match.group(1)) / int(match.group(2))
                blast_percent_id = float(match.group(3))
                blast_hit_count = int(match.group(4)) / int(match.group(5))
                f['description'] = match.group(6)
                f['%'] = u'\u25CF'
                f['A'] = u'\u25CF'
                f['R'] = u'\u25CF'

                if blast_percent_id < 30:
                    style['%'] = 'color: red; text-align: center;'
                elif blast_percent_id < 50:
                    style['%'] = 'color: orange; text-align: center;'
                else:
                    style['%'] = 'color: green; text-align: center;'
                if blast_aligned < 0.5:
                    style['A'] = 'color: red; text-align: center;'
                elif blast_aligned < 0.8:
                    style['A'] = 'color: orange; text-align: center;'
                else:
                    style['A'] = 'color: green; text-align: center;'
                if blast_hit_count < 0.5:
                    style['R'] = 'color: red; text-align: center;'
                elif blast_hit_count < 0.8:
                    style['R'] = 'color: orange; text-align: center;'
                else:
                    style['R'] = 'color: green; text-align: center;'
            if utils.get_feature_qualifier(feature, 'cdd'):
                f['D'] = u'\u25CF'
                style['D'] = 'text-align: center;'
            if f['description']:
                f['description'] = html_make_link(f['id'], f['description'])
            #location
            try:
                number_of_tmh = int(utils.get_feature_qualifier(feature, 'transmembrane_helixes'))
            except ValueError:
                number_of_tmh = 0
            if number_of_tmh > 1:
                f['location'] = 'membrane'
            elif number_of_tmh == 1:
                f['location'] = 'membrane anchor'
            elif utils.get_feature_qualifier(feature, 'signal_peptide'):
                f['location'] = 'envelope'
            elif 'CDS' == feature.type and f['description']:
                f['location'] = 'cytoplasm'
            feature_list.append(f)
            style_list.append(style)
    dataframe_with_data = pd.DataFrame(feature_list, columns=columns)
    dataframe_with_styles = pd.DataFrame(style_list, columns=columns)
    styled_table = dataframe_with_data.style
    styled_table.set_table_styles([{'selector': 'th.col_heading', 'props': 'text-align: center; border-left: 1px solid black;'}])
    styled_table.set_table_styles({'description': [{'selector': 'th.col_heading', 'props': 'text-align: left;'}],
                        'taxon': [{'selector': 'th.col_heading', 'props': 'text-align: left;'}],
                        'id': [{'selector': 'th.col_heading', 'props': 'text-align: left;'}]})
    styled_table.set_table_styles([{'selector': 'td', 'props': 'font-family: Calibri, sans-serif;'}], overwrite=False)
    styled_table.set_table_styles([{'selector': 'th.col_heading', 'props': 'border-bottom: 1px solid black;'}], overwrite=False)
    styled_table.set_sticky(axis=1)
    styled_table.hide(axis="index")
    styled_table.apply(lambda table: dataframe_with_styles, axis=None)
    styled_table.to_html('index_of_features.html', doctype_html=True)


def html_save_all(mag_name, contig_dict):
    metaerg_dir = os.getcwd()
    # load blast results for visualization
    blast_results = {'diamond': predict.spawn_file('diamond', mag_name),
                     'blastn': predict.spawn_file('blastn', mag_name),
                     'cdd': predict.spawn_file('cdd', mag_name)}
    for key in blast_results.keys():
        with utils.TabularBlastParser(blast_results[key]) as handle:
            blast_result_hash = {}
            for br in handle:
                blast_result_hash[br[0]] = br[1]
        blast_results[key] = blast_result_hash
    # make dirs
    html_dir = Path('html')
    html_dir.mkdir(exist_ok=True)
    os.chdir(html_dir)
    mag_dir = Path(mag_name)
    os.chdir(mag_dir)
    # write html
    genome_stats = html_write_genome_stats(mag_name, contig_dict)
    html_write_feature_overview(mag_name, contig_dict, genome_stats)
    for contig in contig_dict.values():
        for feature in contig.features:
            html_write_page_for_feature(feature, contig, blast_results, genome_stats)
    os.chdir(metaerg_dir)