

from metaerg import utils
from metaerg.run_and_read.data_model import MetaergGenome


class AbstractBaseClass:
    def __init__(self, genome: MetaergGenome):
        self.genome = genome

    def _make_html_template(self) -> str:
        """should return the html base for injecting the content in. Returns the html"""
        pass

    def make_html(self) -> str:
        """injects the content into the html base, returns the html"""
        pass

    def write(self):
        """writes the html to a file"""
        pass

    def make_feature_link(self, feature_id, description)-> str:
        return '<a target="_blank" href="features/{}.html">{}</a>'.format(feature_id, description)
