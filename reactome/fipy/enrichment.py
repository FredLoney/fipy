# 2.x compatibility to print to stderr with a function.
from __future__ import print_function
from functools import reduce
import sys
import os
import re
import logging
import requests
from . import rest


PARSERS = {
    'ReactomePathway': lambda p: p,
    'RatioOfProteinInPathway': float,
    'NumberOfProteinInPathway': int,
    'ProteinFromGeneSet': int,
    'P-value': float,
    'FDR': float,
    'HitGenes': lambda s: s.split(',')
}

def enrich(*clusters):
    """
    Perform pathway enrichment analysis on the given clusters.
    The resulting data frame(s) have index *Pathway* and columns
    *p-value* and *FDR*.

    If there is only one input cluster, then the result is the
    enrichments data frame. Otherwise, the result is a list
    of enrichments data frames in the same order as the input
    clusters. If an input cluster could not be enriched, then
    the result for that cluster is `None`.

    :param clusters: the cluster series to enrich
    :return: the result data frame(s)
    """
    enrichments = [cluster.apply(_enrich_genes) for cluster in clusters]
    return enrichments[0] if len(enrichments) == 1 else enrichments


def shared_pathways(enrichments):
    """
    Finds the pathways shared among the given list of
    {module: enrichment} series.

    :param enrichments: the enriched series list
    :return: the common pathways
    """
    intersect = lambda s1, s2: s1.intersection(s2)
    pathways = [_collect_pathways(enriched) for enriched in enrichments]
    return reduce(intersect, pathways)


def transpose_pathways(*enrichments, exclude=None):
    """
    Transposes one or more {module: {pathway: results}} series
    into {pathway: {module: results}} series.

    If there is only one input enrichment, then this method
    returns the transposed series. Otherwise, the return
    value is a list of transposed series in the same order
    as the input enrichments.

    :param enrichments: the enrichments to transpose
    :param exclude: pathways to exclude
    :return: the data frame(s) with pathway index
    """
    transposed = [_transpose_enrichment_pathways(enrichment)
                  for enrichment in enrichments]
    return transposed[0] if len(transposed) == 1 else transposed


def flatten_transposed(column, *transposed):
    """
    Converts the given {pathway: {module: results}} series
    into the {pathway: module} series by selecting the
    nested module with lowest value of the given column.

    If there is only one input series, then this method
    returns the flattened series. Otherwise, the return
    value is a list of flattened series in the same order
    as the inputs.

    :param column: the enrichment result column to filter on
    :param transposed: the transposed series
    :return: the flattend series or list of series
    """
    flattened = [_flatten_transposed_series(column, ds)
                 for ds in transposed]
    return flattened[0] if len(flattened) == 1 else flattened


def distinct(series_list):
    """
    Restricts each of the given series to those index
    values that do not occur in the other series in the
    list. Each series must be indexed solely by pathways.

    :param series_list: the list of pathway-indexed series
    :return: the distinct restrictions
    """
    indexes = [s.index for s in series_list]
    criteria = [_distinct_index(i, indexes)
                for i in range(len(indexes))]
    return [s[criteria[i]] for i, s in enumerate(series_list)]


def exportable(pathways, hierarchy):
    """
    Selects the pathways which have a diagram.

    :param pathways: the pathways to choose from
    :param hierarchy: the Reactome pathway hierarchy
    :yield: the (pathway, database id) tuple
    """
    for pathway in pathways:
        node = _get_hierarchy_node(pathway, hierarchy)
        if node and node['hasDiagram']:
            db_id = node['dbId']
            yield (pathway, db_id)


def export_diagram(pathway, db_id, genes, out_dir=None):
    """
    Exports a diagram PDF for the given pathway. The PDF
    is placed in the target output directory. The file name
    capitalizes the pathway name and removes spaces and
    punctuation, e.g. `DEKBindsTFAP2Homodimers.pdf`.

    :param pathway: the pathway name to export
    :param db_id: the Reactome pathway db id
    :param genes: the genes for the pathway
    :option out_dir: the target directory (default current directory)
    :return: the exported PDF file name
    """
    # Re-enrich the genes in order to get the proper diagram
    # highlighting.
    _enrich_genes(genes)
    if not out_dir:
        out_dir = os.getcwd()
    separator = re.compile('[^\w]+')
    capitalized = [word[0].upper() + word[1:]
                   for word in separator.split(pathway) if word]
    base_name = ''.join(capitalized)
    file_name = os.path.join(out_dir, "%s.pdf" % base_name)
    body = dict(dbId=db_id, pathwayName=pathway, fileName=file_name)
    requests.post(rest.get_fi_url('exportPathwayDiagram'), json=body)
    logging.info("Exported pathway '%s' to %s." % (pathway, file_name))
    return file_name


def _enrich_genes(genes):
    """
    Perform pathway enrichment analysis on the given gene
    list or set. The resulting data frame has index *Pathway*
    and columns *p-value* and *FDR*

    :param genes: the gene list or set to enrich
    :return: the result data frame, or an empty data frame
      if the genes could not be enriched
    """
    data = ','.join(genes)
    # Perform the Reactome enrichment analysis.
    resp = requests.post(rest.get_fi_url('ReactomePathwayEnrichment'),
                         data=data)
    if not resp.ok:
        print("Enrichment unsuccessful: was the Reactome hierarchy loaded?",
              file=sys.stderr)
        resp.raise_for_status()
    parsed = rest.parse_fi_table_response(resp, PARSERS,
                                        index='ReactomePathway')
    # Rename the index and P-value column for consistency.
    parsed.index.rename('Pathway', inplace=True)
    parsed.rename({'P-value': 'p-value'}, inplace=True, axis=1)
    # We only want the p-value and FDR.
    return parsed.loc[:, ['p-value', 'FDR']]


def _get_hierarchy_node(pathway, hierarchy):
    """
    Returns the Reactome hierarchy node for the given pathway.

    :param pathway: the pathway to check
    :param hierarchy: the Reactome pathway hierarchy to check
    :return: the hierarchy node
    """
    if hierarchy['name'] == pathway:
        return hierarchy
    children = hierarchy.get('children')
    if children:
        for subtree in children:
            target = _get_hierarchy_node(pathway, subtree)
            if target:
                return target


def _distinct_index(i, indexes):
    copied = indexes.copy()
    target = copied.pop(i)
    diff = lambda idx1, idx2: idx1.difference(idx2)
    return reduce(diff, copied, target)


def _flatten_transposed_series(column, transposed):
    """
    Converts the given {pathway: {module: results}} series
    into the {pathway: module} series by selecting the
    nested module with lowest value of the given column.

    :param column: the enrichment result column to filter on
    :param transposed: the transposed series
    :return: the flattend series
    """
    grouper = transposed[column].groupby(level='Pathway')
    # Pull the module number from the (pathway, module)
    # maximum index value.
    return grouper.idxmax().apply(lambda idx: idx[1])


def _transpose_enrichment_pathways(enrichment, exclude=None):
    """
    Transposes the given {module: {pathway: results}} series
    into {pathway: {module: results}} series.

    :param enriched: a {module: data frame} series
    :param exclude: pathways to exclude
    :return: the data frame(s) with pathway index
    """
    # The ith module number.
    module = lambda i: enrichment.index.values[i]
    # Inject the module number into each enhanced series.
    dfs = [_add_module_index(df, module(i))
           for i, df in enumerate(enrichment)]
    # Accumulate the multi-indexed enrichments.
    combined = dfs[0].append(dfs[1:])
    if exclude:
        return _exclude_transposed_pathways(combined, exclude)
    else:
        return combined


def _exclude_transposed_pathways(transposed, exclude):
    pathways = transposed.index.droplevel('Module')
    return transposed.loc[~pathways.isin(exclude), :]


def _add_module_index(enriched, module):
    """
    Adds the given module number as the enriched data frame index
    second level.

    :param enriched: a {pathway: result} data frame
    :param module: the module number
    :return: the data frame with (pathway, module) multi-index
    """
    return enriched.assign(Module=module).set_index([enriched.index, 'Module'])


def _collect_pathways(enrichment):
    """
    :param enrichment: the enriched series
    :return: the union of pathways in the series
    """
    pathways = set()
    for df in enrichment.values:
        pathways.update(df.index.get_values())
    return pathways
