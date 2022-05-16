from metaerg.html import abc


class HTMLOneGenoeOverviewPage(abc.AbstractBaseClass):
    def __init__(self, genome):
        super().__init__(genome)

    def make_html(self) -> str:
        '''injects the content into the html base, returns the html'''
        html = self._make_html_template()
        html = html.replace('GENOME_NAME', self.genome.name)
        # genome properties
        html = html.replace('CONTENT_PROPERTIES', ''.join((f'<tr><td>{key}</td><td>{value}</td></tr>\n'
                                                           for key, value in self.genome.properties.items())))
        subsystem_html = ''
        for subsystem in self.genome.subsystems.values():
            subsystem_html += f'<button class="accordion">{str(subsystem)}</button>\n<div class="panel">\n'
            if '[secondary-metabolites]' == subsystem.name:
                if len(subsystem.hits):
                    subsystem_html += '<p><a href="antismash/index.html" target="">View antismash results.</a></p>\n'
                elif not len(subsystem.targets):
                    subsystem_html += '<p>\n'
                    for feature_id in subsystem.hits.values():
                        subsystem_html += '<a target="_blank" href="features/{}.html">{}</a>\n'.format(feature_id, feature_id)
                    subsystem_html += '</p>\n'
            else:
                subsystem_html += '<table>\n'
                for gene, hits in self.genome.subsystems.hits.items():
                    subsystem_html += f'<tr><td>{gene}</td><td>\n'
                    for h in hits:
                        subsystem_html += '<a target="_blank" href="features/{}.html">{}</a>'.format(h, h)
                    subsystem_html += f'</td></tr>\n'
                subsystem_html += '</table>\n'
            subsystem_html += '</div>\n'
        html = html.replace('CONTENT_SUBSYSTEMS', subsystem_html)
        return html


    def _make_html_template(self) -> str:
        '''should return the html base for injecting the content in. Returns the html'''
        return '''<!doctype html>
    <html>
    <head>
        <meta charset="utf-8">
            <title>GENOME_NAME - properties and subsystems</title>\n')
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

    <h4 id=f>Genome properties</h4>
    <table id=f>
        <thead>
            <tr>
              <th>key</th>
              <th>value</th>
            </tr>
        </thead>
        <tbody>
CONTENT_PROPERTIES
        </tbody></table>
        <p><a href="index_of_features.html">View table with all genes</a></p>
        <h4 id=f>Subsystems overview</h4>
CONTENT_SUBSYSTEMS
    <script>
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
    </div></body></html>'''

