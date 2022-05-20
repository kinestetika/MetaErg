from metaerg.run_and_read.data_model import FeatureType
from metaerg.html.abc import HTMLwriter, register

@register
class HTMLFeaturePages(HTMLwriter):
    def __init__(self, genome):
        super().__init__(genome)

    def make_html(self) -> str:
        pass

    def _make_html_template(self) -> str:
        pass
