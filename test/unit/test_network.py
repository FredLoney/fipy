import os
from nose.tools import (assert_true, assert_equal)
from .. import ROOT
from reactome.fipy import network

# The test fixtures directory.
FIXTURES = os.path.join(ROOT, 'fixtures')


class TestNetwork(object):
    """Network build and cluster tests."""

    def test_prepare(self):
        expected = dict(A=27, B=24)
        cohorts = expected.keys()
        clusters = network.prepare(*cohorts, in_dir=FIXTURES)
        for cluster in clusters:
            assert_true(cluster.name in expected,
                        "Cluster name incorrect: %s" % cluster.name)
            actual = cluster.index.size
            assert_equal(expected[cluster.name], actual,
                         "%s cluster size incorrect: %d" %
                         (cluster.name, actual))

    def test_sample_cutoff(self):
        maf_file = os.path.join(FIXTURES, "A.maf")
        network.build_network(maf_file, min_sample_proportion=0.005)
        cluster = network.cluster('A')
        expected = 6
        actual = cluster.index.size
        assert_equal(expected, actual,
                     "Cluster size with cut-off incorrect: %d"
                     % actual)


if __name__ == "__main__":
    import nose
    nose.main(defaultTest=__name__)
