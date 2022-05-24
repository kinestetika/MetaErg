from metaerg.run_and_read.data_model import MetaergGenome
from metaerg.run_and_read.context import Executor

html_registry = []

class HTMLwriter:
    def __init__(self, genome, exec: Executor):
        self.genome = genome
        self.exec = exec

    def _make_html_template(self) -> str:
        """Creates and returns the html base for injecting the content in."""
        pass

    def make_html(self) -> str:
        """Injects the content into the html base, returns the html"""
        pass

    def write_html(self, file=None):
        """writes the html to a file"""
        pass

    def make_feature_link(self, feature_id, description)-> str:
        return '<a target="gene details" href="features/{}.html">{}</a>'.format(feature_id, description)

    def get_color(self, value, thresholds = (80, 60, 40, 20)):
        if value > thresholds[0]:
            return 'id=cg'
        elif value > thresholds[1]:
            return 'id=cb'
        elif value > thresholds[2]:
            return 'id=co'
        elif value > thresholds[3]:
            return 'id=cr'

def register_html_writer(html_writer):
    html_registry.append(html_writer)

    return html_writer
