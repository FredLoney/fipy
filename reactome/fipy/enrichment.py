# 2.x compatibility to print to stderr with a function.
from __future__ import print_function
from functools import reduce
import sys
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

def enrich(modules):
    """
    Perform pathway enrichment analysis on the given gene
    modules. The resulting data frame has index *Pathway*
    and columns *p-value* and *FDR*

    :param modules: the series to enrich
    :return: the result data frame, or an empty data frame
      if the genes could not be enriched
    """
    return modules.apply(_enrich_genes)


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


def transpose_pathways(enriched, exclude=None):
    """
    Transposes the given {module: {pathway: results}} series
    into a {pathway: {module: results}} series.

    :param enriched: a {module: data frame} series
    :param exclude: pathways to exclude
    :return: the data frame with pathway index
    """
    # The ith module number.
    module = lambda i: enriched.index.values[i]
    # Inject the module number into each enhanced series.
    dfs = [_add_module_index(df, module(i))
           for i, df in enumerate(enriched)]
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
    Adds the given module number as the given enriched data
    frame index second level.

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
