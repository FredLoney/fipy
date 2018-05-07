import os
import requests


def download(cohort, out_file=None):
    """
    Downloads the MAF file for the given cancer type cohort.
    The default output file is `<cohort>.maf` in the
    sandbox directory.

    *Note*: As of April 2018, the Firebrowse server is unstable.
    Calling this method is not recommended until the server
    stabilizes.

    :param cohort: the cancer type cohort name
    :option out_file: the optional output file path
    :return: the output file path
    """
    if out_file:
        maf_file = out_file
    else:
        maf_file = os.path.join(os.getcwd(), "%s.maf" % cohort)
    # Download the MAF.
    url = 'http://firebrowse.org/api/v1/Analyses/Mutation/MAF'
    print("Downloading the %s MAF file to %s..." % (cohort, maf_file))
    params = dict(format='tsv', cohort=cohort, page_size=200)
    eof = False
    page = 1
    with open(maf_file, 'w') as f:
        while not eof:
            params['page'] = page
            resp = requests.get(url, params=params)
            if not resp.ok:
                print("Error encountered downloading the %s MAF file. Please retry." %
                      cohort)
                resp.raise_for_status()
            text = resp.text
            if text:
                f.write(text)
                page = params['page'] = page + 1
            else:
                eof = True
            print('+', end='')
    print('')
    print("MAF file downloaded.")

    return maf_file
