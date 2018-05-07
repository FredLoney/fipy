import os
import shutil
from nose.tools import (assert_true, assert_equal)
from .. import ROOT
from reactome.fipy import maf
import csv
import logging

# The test results directory.
RESULTS = os.path.join(ROOT, 'results', 'maf')


class TestMaf(object):
    """MAF download tests."""

    def setUp(self):
        #shutil.rmtree(RESULTS, True)
        #os.makedirs(RESULTS)
        logging.disable(logging.DEBUG)

    def test_download(self):
        cohort = 'OV'
        out_file = os.path.join(RESULTS, "%s.maf" % cohort)
        try:
            maf_file = maf.download(cohort, out_file)
        except Exception as e:
            self.fail("%s download unsuccessful: %s" % (out_file, e))
        assert_equal(out_file, maf_file,
                     "Incorrect output file: %s" % maf_file)
        # Confirm that the file can be parsed.
        with open(maf_file) as f:
            for row in csv.DictReader(f, delimiter='\t'):
                barcode = row['Tumor_Sample_Barcode']
                assert_true("Tumor sample barcode incorrect: %s" % barcode,
                            barcode.startswith('TCGA'))


if __name__ == "__main__":
    import nose
    nose.main(defaultTest=__name__)
