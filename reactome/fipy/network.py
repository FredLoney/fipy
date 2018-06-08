import os
import logging
import requests
import pandas as pd
from .constants import DEF_MIN_SAMPLE_COUNT
from . import (cytoscape, rest)

# The cluster content value parsers.
PARSERS = {'Module': int, 'Node List': lambda s: s.split(',')}


def prepare(*cohorts, **kwargs):
    """
    This convenience function builds the network and clusters
    modules.

    If there is only one input cohort, then the result is the
    cluster data frame. Otherwise, the result is a list
    of cluster data frames in the same order as the input
    cohorts.

    :param cohorts: the cancer type cohort names
    :param kwargs: the following options:
    :option in_dir: the directory with the MAF files
        (default current directory)
    :option min_sample_count: the absolute minimum number
        of samples
    :option min_sample_proportion: the minimum proportion of samples
    :return: the cluster or list of clusters
    """
    in_dir = kwargs.pop('in_dir', os.getcwd())
    clusters = [_prepare_cohort(cohort, in_dir, **kwargs)
                for cohort in cohorts]
    return clusters[0] if len(clusters) == 1 else clusters


def limit_module_size(cluster, cutoff):
    """
    :param cluster: the input cluster series
    :param cutoff: the minimum cluster size
    :return: the cluster series whose modules are at least as
        large as the cut-off
    """
    return cluster.loc[cluster.apply(len).ge(cutoff)]


def _prepare_cohort(cohort, in_dir, **kwargs):
    """
    Utility function to build the network and cluster modules
    for the given cohort.

    :param cohort: the cancer type cohort name
    :param kwargs: the following options:
    :option min_sample_count: the absolute minimum number
        of samples
    :option min_sample_proportion: the minimum proportion of samples
    """
    maf = os.path.join(in_dir, "%s.maf" % cohort)
    build_network(maf, **kwargs)
    return cluster(name=cohort)


def build_network(maf_file, min_sample_count=None, min_sample_proportion=None):
    """
    Builds the network from the given MAF file. Only gene modules
    whose sample size exceeds a threshold are included. The
    threshold is the greater of the minimum module size and
    sample_count * cutoff, where *sample_count* is the total
    number of MAF samples.

    :param maf_file: the downloaded MAF file
    :option min_sample_count: the absolute minimum number
        of samples
    :option min_sample_proportion: the minimum proportion of samples
    """
    if not min_sample_count:
        min_sample_count = DEF_MIN_SAMPLE_COUNT
    if not min_sample_proportion:
        min_sample_proportion = 0

    # Clear the current Cytoscape session, if any.
    cytoscape.clear()
    # Count the samples.
    sample_cnt = _get_sample_count(maf_file)
    # The proportional sample cut-off.
    min_prop_sample_count = int(min_sample_proportion * sample_cnt)
    # The sample cut-off is the larger of the absolute minimum
    # and the proportional minimum.
    sample_cutoff = max(min_sample_count, min_prop_sample_count)
    # Build the network.
    _load_network(maf_file, sample_cutoff)


def cluster(name):
    """
    Cluster the network currently displayed in Cytoscape into
    gene modules.

    :param name: the data series name
    :return: the cluster data series with index *Module*
    """
    # Cluster the currently displayed network.
    logging.info("Clustering the %s FI network..." % name)
    resp = requests.get(rest.get_fi_url('cluster'))
    # Parse the response JSON into a data frame.
    parsed = rest.parse_fi_table_response(resp, parsers=PARSERS,
                                          index='Module')
    logging.info("The %s cluster table is available in Cytoscape." % name)
    genesets = parsed['Node List']
    genesets.name = name

    # Sort by module size.
    sizes = genesets.apply(len).sort_values(ascending=False)
    return genesets.reindex(sizes.index)


def _get_sample_count(maf_file):
    """
    :param: the input MAF file
    :return: the number of samples
    """
    maf_df = pd.read_csv(maf_file, sep='\t',
                         usecols=['Tumor_Sample_Barcode'])
    # The number of samples.
    sample_cnt = len({tsb for tsb in maf_df.Tumor_Sample_Barcode})
    logging.debug("Sample Count: %d" % sample_cnt)

    return sample_cnt


def _load_network(maf_file, sample_cutoff):
        """Loads the network into Cytoscape."""
        logging.info("Building the FI network with sample cut-off %d..." %
              sample_cutoff)
        # Build the FI network with 1% sample cut-off.
        body = dict(fiVersion='2016', format='MAF', file=maf_file,
                    enteredGenes='', chooseHomoGenes=False,
                    userLinkers=False, showUnLinked=False,
                    fetchFIAnnotations=True,
                    sampleCutoffValue=sample_cutoff)
        requests.post(rest.get_fi_url('buildFISubNetwork'), json=body)
        logging.info("The FI network is loaded to Cytoscape.")
