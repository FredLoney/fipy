import os
import numpy as np
from nose.tools import assert_equal
from .. import ROOT
from reactome.fipy import (network, overlap)

# The test fixtures directory.
FIXTURES = os.path.join(ROOT, 'fixtures')


class TestOverlap(object):
    """Overlap analysis tests."""

    def setUp(self):
        self.clusters = network.prepare('A', 'B', in_dir=FIXTURES)

    def test_analyse(self):
        overlaps = overlap.analyse(self.clusters)
        assert_equal(20, overlaps.index.size,
                     "Overlap size incorrect: %d" % overlaps.index.size)
        # Restrict by FDR.
        shared = overlaps[overlaps.FDR.ge(.01)]
        assert_equal(11, shared.index.size,
                     "Shared size incorrect: %d" % shared.index.size)


if __name__ == "__main__":
    import nose
    nose.main(defaultTest=__name__)
