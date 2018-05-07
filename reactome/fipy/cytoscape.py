import requests
from . import rest

def clear():
    """
    Clears the current Cytoscape session, if any.
    """
    try:
        requests.delete(rest.get_url('session')).ok
    except requests.exceptions.ConnectionError:
        print("Error: Connection refused: is Cytoscape started?")
        raise
