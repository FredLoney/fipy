import requests
from . import rest

def load_pathway_hierarchy():
    resp = requests.get(rest.get_fi_url('pathwayTree'))
    return resp.json()['data']

def get_hierarchy_node(pathway, hierarchy):
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
            target = get_hierarchy_node(pathway, subtree)
            if target:
                return target
