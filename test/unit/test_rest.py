import os
import shutil
from nose.tools import assert_equal
from .. import ROOT
from reactome.fipy import rest


class TestRest(object):
    """REST tests."""

    def test_get_url(self):
        expected = 'http://localhost:1234/v1/foo'
        actual = rest.get_url('foo')
        assert_equal(expected, actual,
                     "Incorrect URL: %s" % actual)

    def test_get_fi_url(self):
        expected = 'http://localhost:1234/reactomefiviz/v1/foo'
        actual = rest.get_fi_url('foo')
        assert_equal(expected, actual,
                     "Incorrect URL: %s" % actual)


if __name__ == "__main__":
    import nose
    nose.main(defaultTest=__name__)
