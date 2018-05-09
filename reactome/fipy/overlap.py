from functools import reduce
import pandas as pd
from scipy.stats import hypergeom
from statsmodels.sandbox.stats.multicomp import multipletests
from IPython.display import HTML
from .constants import BACKGROUND_GENE_CNT


def analyse(clusters, n=None):
    """
    Performs pair-wise overlap analysis on the given cluster
    data frames. The resulting data frame has a multi-index
    with levels taken from the inputs and columns *Genes*,
    *p-value* and *FDR*. Rows without overlap are omitted.

    :param clusters: the cluster series collection
    :param n: the background gene count
        (default :const:`reactome.fipy.constants.BACKGROUND_GENE_CNT`)
    :return the overlap data frame
    """
    if not n:
        n = BACKGROUND_GENE_CNT
    # Make the overlap multi-index.
    index_names = [cluster.name for cluster in clusters]
    index_values = [cluster.index.values for cluster in clusters]
    multi_index = pd.MultiIndex.from_product(index_values,
                                             names=index_names)
    # Take the pair-wise intersections.
    genesets = [[set(gs) for gs in cluster.values]
                for cluster in clusters]
    data = [gs1.intersection(gs2)
            for gs1 in genesets[0]
            for gs2 in genesets[1]]
    intersections = pd.Series(data=data, index=multi_index)
    # Omit the empty intersections.
    filtered = intersections[intersections.apply(len).gt(0)]

    # Combine the geneset series.
    shared = filtered.to_frame(name='Shared')
    renamed = [ds.rename_axis(ds.name) for ds in clusters]
    join = lambda accum, ds: accum.join(ds)
    combined = reduce(join, renamed, shared)

    # Calculate the p-values.
    pvalue = lambda row: _overlap_pvalue(n, *row.values)
    pvalues = combined.apply(pvalue, axis=1)

    # Correct the p-values for multiple comparison hypothesis
    # testing by applying the Benjamini–Hochberg FDR procedure.
    _, fdrs, _, _ = multipletests(pvalues.values, method='fdr_bh')

    # Assemble the data frame.
    cols = {'p-value': pvalues, 'FDR': fdrs}
    return shared.assign(**cols)


def limit_overlap(overlaps, column, cutoff):
    """
    Restricts the given overlaps data frame to those module
    pairs whose FDR is no greater than the given cut-off.
    
    :param overlaps: the overlaps data frame
    :param column: the overlaps column name (`FDR` or `p-value`)
    :param cutoff: the threshold value
    """
    return overlaps[overlaps[column].le(cutoff)]


def unshared(clusters, shared):
    """
    Finds each cluster modules subset which is not
    in the given shared overlap data frame. The input
    clusters must be in the same order as the analyse
    *clusters* argument.

    :param clusters: the list of cluster data series
    :param shared: the shared overlap data frame
    :return: the filtered clusters
    """
    return [_difference(cluster, shared.index, i)
            for i, cluster in enumerate(clusters)]


def _difference(cluster, multi_index, level):
    """
    Finds the subset of the given cluster whose index is
    not in the given index level.

    :param cluster: the input cluster
    :param multi_index: the multi-index to select from
    :param level: the level number
    :return: the filtered cluster
    """
    shared = multi_index.get_level_values(level).unique()
    return cluster[cluster.index.difference(shared)]


def _index_level_values(multi_index, level):
    """
    :param multi_index: the multi-index to select from
    :param level: the level number
    :return: the level values
    """
    return multi_index.get_level_values(level).unique()


def print_overlap(overlaps, qualifier=None, format='html'):
    """
    Represents the given overlap data frame as a string
    suitable for printing.

    :param overlaps: the overlap data frame
    :option qualifier: the optional qualifier to print
        in the HTML-formatted title
    :option format: 'html' or 'text' (default is html)
    :return: the printable string
    """
    flat = overlaps.reset_index()

    columns = [col for col in flat.columns if not col.endswith("Genes")]
    cohort_col_grps = [[col for col in columns if col.startswith(name)]
                       for name in overlaps.index.names]
    cohort_cols = reduce(lambda x,y: x + y, cohort_col_grps)
    non_cohort_cols = [col for col in columns if col not in cohort_cols]
    ordered_cols = cohort_cols + non_cohort_cols
    printable = flat.loc[:, ordered_cols]
    if format == 'text':
        return printable.to_string(index=False)
    elif format == 'html':
        start = '<h4>Module Overlap'
        if qualifier:
            style = 'font-size:normal;font-weight:normal;'
            middle = "<span style=%s> (%s)</span>" % (style, qualifier)
        else:
            middle = ''
        end = '</h4>'
        heading = start + middle + end
        return HTML(heading + printable.to_html(index=False))
    else:
        raise ValueError("Unrecognized overlap print format: %s" % format)

def _overlap_pvalue(N, overlap, gs1, gs2):
    """
    :param N: the number of background genes
    :param overlap: the genes in common
    :param gs1: the first gene set
    :param gs2: the second gene set
    :return: the probability of gene set overlap from
      a background population of *N* genes
    """
    # The hypergeometric distribution parameters.
    k = len(overlap) - 1
    n1 = len(gs1)
    n2 = len(gs2)
    # Return the overlap probability.
    return hypergeom.sf(k, N, n1, n2)
