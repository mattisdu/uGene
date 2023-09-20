"""
uGeneCore.py
Version: 0.5
Beschreibung: This is the core file of uGene and performs the uMap clustering. All inputs being made by command line input.
              Example: $> python uGeneCore.py -v -f exampleData.phyloprofile -t [{'y_axis':'geneID','x_axis':'ncbiID','jobs':'gene','values':['FAS_F','FAS_B']}]

Autor: Mattis Kaumann

The MIT License (MIT)

Copyright (c) 2023 Mattis Kaumann, Goethe-Universität Frankfurt am Main
Read more on LICENSE.txt
"""
UGENE_CORE_VERSION = "0.5"

import pandas as pd
import json
import logging
import argparse
import warnings

# Depress numba warning before import umap.
from numba.core.errors import NumbaDeprecationWarning

# Since change form python 3.10 to 3.11 the current numba version throws the following warning:
# NumbaDeprecationWarning: The 'nopython' keyword argument was not supplied to the 'numba.jit' decorator.
# The implicit default value for this argument is currently False, but it will be changed to True in Numba 0.59.0.
warnings.simplefilter('ignore', category=NumbaDeprecationWarning)
import umap


def printHelp(arg_parser):
    print("\n------------ Help ------------")
    arg_parser.print_help()
    print("Full Example:")
    print("\tWindows-> python uGeneCore.py -f yourFile.csv -t [{'y_axis':'geneID','x_axis':'ncbiID','jobs':'gene'," +
          "'values':['FAS_F','FAS_B']}]")
    print("\tLinux-> python3 uGeneCore.py -f yourFile.csv -t \"[{'y_axis':'geneID','x_axis':'ncbiID','jobs':'gene'," +
          "'values':['FAS_F','FAS_B']}]\"")
    print("About one task:")
    print("\tOne task is shall be one dictionary with the following keys:")
    print("\t y_axis : Strict Necessary! Column name which column will become y_axis of data matrix.")
    print("\t x_axis : Strict Necessary! Column name which column will become x_axis of data matrix.")
    print("\t value : Strict Necessary! Column name or list of column names within used raw data.")
    print("\t jobs : Strict Necessary! List of jobs or just a job name for standard use. Standard: 'gene' | 'taxa' .")
    print("\t dev_report : Enable report about dataframe key conflicts. Case of key y_axis and x_axis are duplicated.")
    print("\t fill : Float to fill nan into the data matrix.")
    print("\t pattern : Enable pattern mode. Scores will be translated into one zero pattern. Overrides fill value.")
    print("\t drop : Enable drop mode for dataframe key conflicts.")
    print("About jobs:")
    print("\tYou be abele to hand over a list of Jobs or one job name.")
    print("\tOne job name will converted to a job list of three jobs with your job name and the n_components 1, 2, 3.")
    print("\tInto a job you be able to set sklearn umap cluster parameter. All parameter of clusterDF() are supported.")
    print("------------------------------")


def jsonLoads(json_str, error_return):
    """
    Save use of json.load() to convert a json string to an object.
    :param str json_str: A simpy json string
    :param object error_return: A object which one should get returned in a case of an error.
    :return: Give back the parsed json sting or the error_return object.
    """
    try:
        json_str = json_str.replace("'", '"')
        return json.loads(json_str)
    except:
        logging.error("Can´t interpret" + str(json_str) + " arguments.")
    # Error case gives origin string back.
    return error_return


def clusterDF(df_mat, job_name='unnamedJob', angular_rp_forest=False, b=None,
              force_approximation_algorithm=False, init='spectral', learning_rate=1.0,
              local_connectivity=1.0, low_memory=False, metric='euclidean',
              metric_kwds=None, min_dist=0.1, n_components=2, n_epochs=None,
              n_neighbors=15, negative_sample_rate=5, output_metric='euclidean',
              output_metric_kwds=None, random_state=42, repulsion_strength=1.0,
              set_op_mix_ratio=1.0, spread=1.0, target_metric='categorical',
              target_metric_kwds=None, target_n_neighbors=-1, target_weight=0.5,
              transform_queue_size=4.0, transform_seed=42, unique=False, verbose=False):
    """
    Perform the main cluster job with the uMap cluster algorithm form umap-sklearn.
    Most of these descriptions are more or less taken from the Sklearn-Umap-Wiki. Read the origin wiki for more details.
    :param pandas df_mat: Takes a pandas dataframe except only filled with raw data and organised like a matrix.
    :param str job_name: Just a name of this job. It is needed to create the result column name. Use one name just once.
    :param bool angular_rp_forest: Enable an angular random forest projection to initialize nearest neighbor search.
    :param float b: Controls the embedding. If not given b will be scaled by min_dist.
    :param force_approximation_algorithm:
    :param string init: String to set up low dimensional embedding.
    :param float learning_rate: Float for initial learning rate.
    :param int local_connectivity: Adjust local manifold connectivity.
    :param bool low_memory: Avoids excessive memory use for umap as well for multiprocessing jobs.
    :param str metric: Means the metric which one calculates the distance of data from each other. Ex.'euclidean'
    :param dict metric_kwds: Additional required for some metrics.
    :param float min_dist: Minimal distance of to data cluster points.
    :param int n_components: Defines how many dimensions this uMap should have.
    :param int n_epochs: Number of training epochs.
    :param int n_neighbors: Heuristic value, which effect a global vs local vision of the algorithm.
    :param negative_sample_rate: Amount of negative samples per positive selected sample.
    :param function output_metric: Post process metric function.
    :param dict output_metric_kwds: Arguments for the metric function.
    :param int|array random_state: Numpy random state seed.
    :param repulsion_strength:
    :param float set_op_mix_ratio: Parameter between 0.0 and 1.0 to interpolate between Unions.
    :param float spread: Scale spread of points. Effects how the cluster clumps like min_dist.
    :param string target_metric: Metric to measure distance for a target arrays.
    :param dict target_metric_kwds: Arguments for the target metric.
    :param int target_n_neighbors: Count of neighbors for target simplicity. Set -1 to use the main n_neighbors.
    :param float target_weight: Weighting factor between data topology and target topology.
    :param float transform_queue_size: Control how aggressively perform nearest neighbors search.
    :param int transform_seed: Seed for transform operations.
    :param bool unique: Enable support for high rate of duplicates.
    :param bool verbose: Enables print out of logging information.
    :return list: List within tuple of len = 2. First item represent column name and the second keeps cluster data.
    """
    try:
        umap_result = umap.UMAP(
            angular_rp_forest=angular_rp_forest, b=b,
            force_approximation_algorithm=force_approximation_algorithm,
            init=init, learning_rate=learning_rate,
            local_connectivity=local_connectivity,
            low_memory=low_memory, metric=metric,
            metric_kwds=metric_kwds, min_dist=min_dist,
            n_components=n_components, n_epochs=n_epochs,
            n_neighbors=n_neighbors, negative_sample_rate=negative_sample_rate,
            output_metric=output_metric,
            output_metric_kwds=output_metric_kwds, random_state=random_state,
            repulsion_strength=repulsion_strength,
            set_op_mix_ratio=set_op_mix_ratio, spread=spread,
            target_metric=target_metric,
            target_metric_kwds=target_metric_kwds, target_n_neighbors=target_n_neighbors,
            target_weight=target_weight,
            transform_queue_size=transform_queue_size, transform_seed=transform_seed,
            unique=unique, verbose=verbose,
        ).fit_transform(df_mat)

    except Exception as error:
        logging.error(error)
        # An empty list will have no effect because the list of results get iterated.
        return []

    # Save the UMAP output data
    cord = ['x', 'y', 'z', 'a', 'b', 'c']

    # Iterate over umap results. Give back a list within tuple. First column name, second data vector.
    res = [(job_name + str(n_components) + 'd_' + it[0], it[1]) for it in zip(cord, umap_result.transpose())]
    return res


def pivotDF(df, index, column, values, fill=float("nan"), report=False, drop=False, multilevel=False, pattern=False):
    """
    Manage the pandas.pivot() calculation. Allows column combinations and do conflict management. Index and column are
    somthing like key in a dataframe. Conflicts appear whenever these keys are not unique.

    :param df: Standard molten dataframe.
    :param list index: List of column names which one should affect the y-axis. Cluster subjects.
    :param list column: List of column names which one should affect the x-axis.
    :param values: List with columns which containing raw date for the matrix. All these columns should contain values.
    :param float|int fill: Value to fill not defined values into the data matrix. Depressed by pattern.
    :param bool report: Simple flac. If True, a simple report about conflict into the dataframe get printed.
    :param bool drop: If True duplicated will be cut. This could in some cases improve the runtime.
    :return: Returns a pandas dataframe only filled with raw data and organised like a matrix.
    :param bool multilevel: Enables origin multi level column structure. Disable avoids multilevel columns.
    :param bool|float pattern: Enable convert data to presence absence patter.
    """
    #  Cut unimportant data columns
    df = df.drop([col for col in df.columns if col not in index + column + values], axis=1)

    # Report conflict deviation. Conflicts are given by the chosen key combination.
    # To prevent conflicts, the index and column must be genuine keys.
    if report:
        groups_df = df.groupby(index + column)
        # Report max deviation, mean deviation and count of conflicts for the given genuine key.
        std_df = groups_df.std()
        size_series = groups_df.size()
        logging.critical("\n->Deviation report for database conflicts:")
        logging.critical("\tChosen keys : " + str(index) + "\t" + str(column))
        logging.critical("\tTotal count conflicts : " + str(size_series[size_series > 1].count()))
        logging.critical("\tMean of standard deviation : " + std_df.mean().to_string())
        logging.critical("\tTotal max deviation : " + std_df.max().to_string())
        if not drop:
            df = groups_df.mean().reset_index()
        else:
            df = df.drop_duplicates(subset=index + column).reset_index(drop=True)
    else:
        # Reduce Conflicts without doing a conflict report.
        if not drop:
            df = df.groupby(index + column).mean().reset_index()
        else:
            df = df.drop_duplicates(subset=index + column).reset_index(drop=True)

    df = df.pivot(index=index, columns=column, values=values)

    if not multilevel:
        # Clean up with multilevel columns
        while df.columns.nlevels > 1:
            df = df.droplevel(0, axis=1)

        # Reset all column names to a distinct number
        df.columns = range(len(df.columns))

    if pattern:
        if type(pattern) != float:
            df = df.applymap(lambda x: True if x > 0 else False)
        # Floats are allowed to define a cutoff border.
        else:
            df = df.applymap(lambda x: True if x > pattern else False)

    elif fill == fill and bool(fill):
        # Fill all nan cells with fill value. If pattern is used it will handle nan to avoid dataframe operations.
        df = df.fillna(fill)

    return df


def mainAnalytics(df, x_axis=[], y_axis=[], values=[], jobs=[], dev_report=False, fill=-1, pattern=False, drop=False):
    """
    This is a kind of main function which one is called for every task. All over it effect program and processing
    structure.

    :param df: Standard molten pandas dataframe.
    :param str|list x_axis: One column name or list of column names which one get used to pivot.
    :param str|list y_axis: One column name or list of column names which one get used to pivot.
    :param str|list values: One column name or list of column names which one get used to pivot.
    :param str|list jobs: One sting which effect the job_name or a job list within a dictionary for each job.
    :param bool dev_report: Simple flac. If True, a simple report about conflict into the dataframe get printed.
    :param float|int fill: Value to fill not defined values into the data matrix. Depressed by pattern.
    :param bool|float pattern: Enable convert data to presence absence patter.
    :param bool drop: If True duplicated will be cut. This could in some cases improve the runtime.
    :return: Returns a standard molten pandas dataframe with all results.
    """

    # Check inputs
    if len(x_axis) == 0 or len(y_axis) == 0 or len(values) == 0 or len(jobs) == 0:
        raise Exception("Error mainAnalytics()! Unable arguments x_axis, y_axis, values or jobs with len zero. ")

    # Process x_axis, y_axis and values given by a simple string. After these have to
    if type(x_axis) == str:
        x_axis = [str(x_axis)]
    if type(y_axis) == str:
        y_axis = [str(y_axis)]
    if type(values) == str:
        values = [str(values)]
    if type(jobs) == str:
        # Enable standard cluster jobs by one string
        jobs = [{'job_name': str(jobs), 'n_components': 1},
                {'job_name': str(jobs), 'n_components': 2},
                {'job_name': str(jobs), 'n_components': 3}]

    # Verify that all axes are lists and contain the same data type.
    if not type(y_axis) == type(x_axis) == type(values) == type(jobs) == list:
        raise Exception("Error mainAnalytics()! Wrong data type of given arguments.")

    # Check if all used column names are present into the dataframe.
    if not all([col_name in df.columns for col_name in x_axis + y_axis + values]):
        raise Exception("Error mainAnalytics()! Hand over column names are not present into the dataframe.")

    mat_df = pivotDF(df, column=x_axis, index=y_axis, values=values, fill=fill, report=dev_report, pattern=pattern,
                     drop=drop)

    # Do all the cluster jobs. All new data columns will be stored into res_df
    res_df = pd.DataFrame(index=mat_df.index)

    # Calculate all jobs
    for it_job in jobs:
        data_job = clusterDF(mat_df, **it_job)

        # Keep in mind. One column is a tuple with (column name, data list).
        for data_col in data_job:
            res_df[data_col[0]] = data_col[1]

    # Drop all column level. Old version bug fix and a case which should never appear now.
    while res_df.columns.nlevels > 1:
        logging.error("Drop res_df multi level columns.")
        res_df = res_df.droplevel(1, axis=1)

    if res_df.empty:
        logging.error("No cluster data is produced in total.")
        return df

    # Override duplicated column names.
    df = df.drop([col for col in res_df.columns if col in df.columns], axis=1)

    # Merge results into the origin database.
    df = df.merge(res_df, left_on=y_axis, right_index=True, how="left")

    return df


def main():
    # Process given arguments.
    arg_parser = argparse.ArgumentParser(description="uGeneCore.py", conflict_handler="resolve")
    arg_parser.add_argument('-f', "--file", type=str, default=None,
                            help="Filename with path to a .csv file of data. Ex: ./exampleData.csv")
    arg_parser.add_argument('-l', '--logfile', type=str, default=None, help="Path with filename to a given logfile.")
    arg_parser.add_argument('-t', '--tasks', type=str, default="[]",
                        help="List of Tasks. Ex:[{'y_axis':'geneID','x_axis':'ncbiID','values':'FAS_F','jobs':'gene'}]")
    arg_parser.add_argument('-v', '--verbose', action='store_true', help="Enable verbose logging. ")
    arg_parser.add_argument('-h', '--help', action='store_true', help="Show help thread.")
    arg_parser.add_argument('-i', '--info', '--version', action='store_true', help="Show current uGeneCore version.")

    args = arg_parser.parse_args()

    # Setup logging system. Depending on logfile and verbose option.
    if args.verbose:
        logging.basicConfig(
            filename=args.logfile,
            filemode='w',
            level=logging.DEBUG
        )
    else:
        logging.basicConfig(
            filename=args.logfile
        )

    logging.info("\n--- Start uGeneCore.py --- ")

    if args.info:
        print("uGeneCore Version: ", UGENE_CORE_VERSION)

    if args.help:
        printHelp(arg_parser)

    # Argument post processing
    file_name = args.file
    tasks = jsonLoads(args.tasks, [])

    logging.debug("Filename: " + str(file_name))
    logging.debug("Tasks: " + str(tasks))

    # ----------------------- Debug ------------------------------------------
    # Bypass the command line input with the following lines.
    # job_name='unnamedJob', n_components=1, n_neighbors=15, min_dist=0.1, metric='euclidean'
    # file_name = "PipelineData/F297.csv"
    # tasks = [{'dev_report': 1, 'x_axis': 'ncbiID', 'y_axis': 'geneID', 'values': 'FAS_F', 'jobs': 'gene'}]

    # Finish by no given tasks.
    if not tasks:
        logging.warning("\n--- Exit uGeneCore.py by no tasks --- ")
        return

    # Analytics processing
    try:
        logging.debug("Load Data " + str(file_name))

        if file_name.endswith(".phyloprofile"):
            df = pd.read_csv(file_name, sep="\t")
            file_name = file_name[:-len(".phyloprofile")] + ".csv"
        else:
            df = pd.read_csv(file_name)

        logging.debug(df.head())

    except Exception as error:
        logging.error("Cant´s open '" + str(file_name) + "' file.\n")
        logging.error(error)
        logging.warning("\n--- Exit uGeneCore.py --- ")
        return;

    for cur_task in tasks:
        try:
            logging.debug("Process task : " + str(cur_task))
            df = mainAnalytics(df, **cur_task)

        except Exception as error:
            logging.error("Task " + str(cur_task) + " fail! \t")
            logging.error(error)

    df.to_csv(file_name.replace(".csv", ".cluster.csv"), index=False)
    logging.info("\n--- Finish uGeneCore.py --- ")


if __name__ == "__main__":
    main()
