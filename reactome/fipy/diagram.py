def export_diagram(db_id, pathway, genes, out_dir=None):
    """
    Exports a diagram PDF for the given pathway. The PDF
    is placed in the target output directory. The file name
    capitalizes the pathway name and removes spaces and
    punctuation, e.g. `DEKBindsTFAP2Homodimers`.

    :param db_id: the Reactome pathway db id
    :param pathway: the pathway name to export
    :param genes: the hit genes for the pathway
    :option out_dir: the target directory (default current directory)
    :return: the exported PDF file name
    """
    # Re-enrich the genes in order to get the proper diagram
    # highlighting.
    enrich_genes(genes)
    if not out_dir:
        out_dir = os.getcwd()
    separator = re.compile('[^\w]+')
    capitalized = [word[0].upper() + word[1:]
                   for word in separator.split(pathway) if word]
    base_name = ''.join(capitalized)
    file_name = os.path.join(out_dir, "%s.pdf" % base_name)
    body = dict(dbId=db_id, pathwayName=pathway, fileName=file_name)
    requests.post(get_fi_url('exportPathwayDiagram'), json=body)
    print("Exported pathway '%s' to %s." % (pathway, file_name))
    return file_name
