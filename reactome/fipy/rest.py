import pandas as pd
from .constants import SERVICE_URL


def get_url(*params, **kwargs):
    """
    :param params: the URL path component strings
    :option app: the CyREST application name
    :return: the REST app URL
    """
    path = [SERVICE_URL]
    app = kwargs.get('app')
    if app is not None:
        path.append(app)
    path.append('v1')
    path.extend(params)
    return '/'.join(path)


def get_fi_url(*params):
    """
    A convenience function to format a URL path for
    the `reactomefiviz` application.

    :param params: the URL path component strings
    :return: the CyREST Reactome FI url
    """
    return get_url(*params, app='reactomefiviz')


def parse_fi_table_response(resp, parsers, index=None):
    """
    Returns the data frame for the given CyREST Reactome
    FI response. The response `data` property object must
    be a JSON object with properties *tableHeaders* and
    *tableContent*. If the response `data` is empty, then
    this function returns an empty data frame with columns
    given by the *parsers* keys.

    The required *parsers* dictionary argument associates
    a parser with each column.

    :param resp: the CyREST response
    :param parsers: the column parser dictionary
    :option index: the index column name
    :return: the parsed data frame
    """
    # The response JSON data object.
    data = resp.json()['data']
    # The data columns.
    columns = data.get('tableHeaders') if data else parsers.keys()
    # The default "parser".
    identity = lambda value: value
    parsers_list = [parsers.get(col, identity) for col in columns]
    # Parses a content row.
    parse_row = lambda row: tuple(parsers_list[i](value)
                                  for i, value in enumerate(row))
    # The parsed content list.
    content = map(parse_row, data['tableContent']) if data else []
    # Return the data frames.
    return pd.DataFrame.from_records(content, index=index, columns=columns)
