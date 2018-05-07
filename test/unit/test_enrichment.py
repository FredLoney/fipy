import os
import pandas as pd
from nose.tools import assert_equal
from .. import ROOT
from reactome.fipy import (reactome, enrichment)

# The test fixtures directory.
FIXTURES = os.path.join(ROOT, 'fixtures')


class TestEnrichment(object):
    """Overlap analysis tests."""

    def setUp(self):
        index = [0, 3, 5]
        genes = [
            ['AKNA', 'ANKRD1', 'BMPR2', 'CNOT10', 'EIF4G1'],
            ['FANCB', 'FOXA1', 'FRYL', 'GATA3', 'GLI2'],
            ['GYG2', 'HERC2', 'HIVEP1', 'JAG1', 'KLHL8']
        ]
        self.clusters = pd.Series(data=genes, index=index)
        reactome.load_pathway_hierarchy()

    def test_enrich(self):
        enrichments = enrichment.enrich(self.clusters)
        expected = [38, 13, 26]
        for i, enriched in enumerate(enrichments):
            assert_equal(expected[i], enriched.index.size,
                     "%s enrichment size incorrect: %d" %
                     (enriched.values, enriched.index.size))


if __name__ == "__main__":
    import nose
    nose.main(defaultTest=__name__)
