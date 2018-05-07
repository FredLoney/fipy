#! /usr/bin/env python
import sys
import os
import argparse
import pandas as pd
from reactome.fipy import (
    cytoscape, network, overlap, reactome, enrichment, diagram
)

"""The workflow default option values."""
DEFAULTS = dict(
    min_sample_proportion=.01,
    max_module_count=20, min_module_size=3,
    max_overlap_fdr=.01, max_enrichment_fdr=.001
)


def run(*cohorts, **kwargs):
    """
    Executes the workflow example shown in the corresponding
    (notebook)[https://github.com/reactome-fi/workflows/blob/master/CytoscapeReactomeFI.ipynb].
    This method generalizes the example for any two cohort
    inputs and customized threshold options.

    :param cohorts: the cancer type names
    :param kwargs: the following keyword arguments:
    :option sandbox: the directory containing the MAF files
      and generated diagrams
    :param kwargs: the following options:
    :option min_sample_count: the minimum module size
    :option min_sample_proportion: the minimum proportion
         of samples
    :option max_module_count: the maximum number of modules
    :option min_module_size: the minimum number of genes
         per module
    :option max_overlap_fdr: the overlap FDR cutoff
    :option max_enrichment_fdr: the enrichment FDR cutoff
    """

    # Add the option defaults.
    kwargs.update(DEFAULTS)

    # The sandbox directory contains the MAF files and generated
    # diagrams. The default is the working directory.
    sandbox = kwargs.get('sandbox', os.getcwd())

    # Prepare the inputs for analysis.
    prep_opts = ('min_sample_count', 'min_sample_proportion',
                 'max_module_count', 'min_module_size')
    prep_kwargs = {k: kwargs.get(k) for k in prep_opt}
    clusters = [_prepare(k, sandbox, **prep_opts) for k in cohorts]
    for cluster in clusters:
        print("%s module count: %d" % (cluster.name, cluster.index.size))

    # Get the overlap.
    overlaps = overlap.analyse(clusters)
    # Restrict by FDR cutoff.
    overlaps = overlaps[overlaps.FDR.ge(kwargs['max_overlap_fdr'])]
    # Print the overlap.
    print(overlap.print_overlap(overlaps,
                                cutoff=opts.get('overlap_threshold'),
                                format='text'))



    # Partition into shared and unshared.
    partition = overlap.partition_shared(overlap_df, inputs)

    # Load the Reactome hierarchy into Cytoscape.
    reactome.load_pathway_hierarchy()

    # Perform pathway enrichment on the gene modules.
    rename_opts = {'Shared Genes': 'Genes'}
    columns = ['Genes', 'p-value', 'FDR']
    shared_df = partition['shared'].rename(columns=rename_opts).loc[:, columns]
    groups = dict(Shared=shared_df)
    groups.update(partition['unshared'])
    enriched = {name: enrichment.enrich_modules(name, df.Genes)
                for name, df in groups.items()}
    enriched_unshared = {name: enriched[name] for name in inputs.keys()}
    distinct = enrichment.distinct_pathways(enriched_unshared)

    # Display a sample module.
    cancer = inputs.keys()[0]
    # The enriched module numbers.
    module_numbers = enriched[cancer].index.get_values()
    print("Enriched %s module numbers:" % cancer)
    print(module_numbers)
    # Show the enriched pathways for the sample module.
    module_number = module_numbers[0]
    df = enriched[cancer].loc[module_number]
    print("Module %d Genes:" % module_number)
    print(','.join(df.Genes))
    print()
    print("Module %d Pathways:" % module_number)
    print(df.Pathways.loc[:, ['p-value', 'FDR']]
            .reset_index()
            .to_string(index=False))

    # Export one diagram per cancer type.
    for name, pathways in distinct.items():
        print("Exporting a %s pathway..." % name)
        df = enriched[name]
        diagram.export_first_pathway_diagram(pathways, df, reactome_tree)


def main(argv):
    # There are no positional arguments.
    cohorts, opts = _parse_arguments()
    run(*cohorts, **opts)


def _parse_arguments():
    """
    Parses the command line arguments.

    :return: the (positional, options) tuple, where *positional*
        is the non-option positional arguments and *options* is
        an {option: value} dictionary
    """
    parser = argparse.ArgumentParser()

    # The log options.
    parser.add_argument('-l', '--log', help='the log file (default stdout)',
                        metavar='FILE', dest='log_file')
    verbosity_grp = parser.add_mutually_exclusive_group()
    verbosity_grp.add_argument(
        '-q', '--quiet', help="only log error messages", dest='log_level',
        action='store_const', const='ERROR')
    verbosity_grp.add_argument(
        '-d', '--debug', help='log debug messages', dest='log_level',
        action='store_const', const='DEBUG')

    # The output option.
    parser.add_argument('cohorts', metavar='name', type=str, nargs=2)

    # The sandbox directory option.
    parser.add_argument('-s', '--sandbox', metavar='DIR', dest='sandbox',
                        help='the directory holding MAF and diagram files'
                             ' (default is the current working directory)')
    # Tuning parameters.
    parser.add_argument('--cluster-threshold', type=int, metavar='FLOAT',
                        help='the cluster p-value threshold (default .01)')
    parser.add_argument('--max-module-count', type=int, metavar='INT',
                        help='the maximum number of cluster modules'
                             ' (default 20)')
    parser.add_argument('--min-module-size', type=int, metavar='INT',
                        help='the minimum cluster module size'
                             ' (default 3)')

    args = vars(parser.parse_args())
    nonempty_args = dict((k, v) for k, v in args.items() if v != None)

    return nonempty_args.pop('cohorts'), nonempty_args


def _prepare(cohort, maf_dir, min_sample_count, min_sample_proportion,
             max_module_count, min_module_size):
    """
    Prepares the given cancer type cohort for analysis.

    :param cohort: the cancer type cohort name
    :param maf_dir: the MAF directory
    :return: the clustered data frame
    """
    maf_file = os.path.join(maf_dir, "%s.maf" % cohort)
    # Download the MAF, if necessary.
    if (not os.path.exists(maf_file)):
        # Download disabled due to Firebrowse server instability.
        #download_maf(cohort, out=maf_file)
        raise IOError("The %s MAF file was not found: %s" %
                      (cohort, maf_file))

    # Build the network.
    network.build_network(maf_file, min_sample_count, min_sample_proportion)
    # Cluster into gene modules.
    cluster = network.cluster()

    # Restrict by module count.
    if max_module_count < cluster.index.size:
        cluster = cluster[:max_module_count]
    # Restrict by module size.
    return cluster[cluster.apply(len).ge(min_module_size)]


if __name__ == "__main__":
    main(sys.argv[1:])
