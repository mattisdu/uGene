import pandas as pd
import uGeneCore as uGene
import numpy as np
import random
import time
import plotly.express as px
import os

# Not working on Windows.
if os.name != "nt":
    import resource


def figuresToHtml(figs, html_add="", filename="dashboard.html"):
    """ Save plotly figures to html file.
    :param list figs: List with in plotly figures
    :param str html_add: Any kind of html text which should be added to the html body.
    :param str filename: Filename and path of the html-file.
    :return void:
    """
    # https://stackoverflow.com/questions/45577255/plotly-plot-multiple-figures-as-subplots
    with open(filename, 'a', encoding='utf-8') as dashboard:
        dashboard.write("<html><head></head><body>" + "\n")
        dashboard.write(html_add)
        for fig in figs:
            inner_html = fig.to_html().split('<body>')[1].split('</body>')[0]
            dashboard.write(inner_html)
        dashboard.write("</body></html>" + "\n")


def cutValue(value):
    """Function which adjust random values into the most common uGene data values. FAS scores are numbers between zero
    and one. No ortholog will calculate as -1.0. This effects no range between zero and minus one.
    :param float value: Any float.
    :return: Float between zero and one or minus one."""
    return value if 0 <= value < 1.0 else float("nan") if value < 0 else 1.0


def newSample(y, x):
    """ Simulate random data of a phylogenetic profile.
    :param int y: Amount of generated genes.
    :param int x: Amount of generated taxa.
    :return: Dataframe with the following columns. ["geneID", "ncbiID", "FAS_F", "FAS_B"]
    """
    # Create a y times x matrix within random numbers.
    df = pd.DataFrame(
        np.random.randn(y, x),
        columns=["taxa_" + str(it_x) for it_x in range(x)],
        index=["gene_" + str(it_y) for it_y in range(y)]
    )
    # Correct values standard occurring values.
    df = df.applymap(cutValue)
    df = df.reset_index(names=['geneID'])
    df = df.melt(id_vars=['geneID'], value_vars=df.columns[1:]).dropna()
    df = df.rename(columns={'variable': 'ncbiID', 'value': 'FAS_F'})
    # FAS_F and FAS_B scores are different sometimes but just by chance. Not the other way around, similar by chance.
    df['FAS_B'] = df['FAS_F'].apply(lambda it: cutValue(it + (random.random() - 0.5) * 0.9))
    return df


def subSample(df, y=None, x=None):
    """ Simulate random sub samples of phylogenetic profile. Columns ["geneID", "ncbiID"] are required.
    :param obj df: Pandas dataframe of the full phylogenetic profile
    :param int y: Amount of generated genes.
    :param int x: Amount of generated taxa.
    :return: Dataframe with random sub sampled number of geneID and ncbiID.
    """
    if y:
        gene_ids = df['geneID'].drop_duplicates()
        gene_ids = gene_ids.sample(y) if len(gene_ids) > y else gene_ids
        df = df[df['geneID'].isin(gene_ids)].reset_index(drop=True)
    if x:
        taxa_ids = df['ncbiID'].drop_duplicates()
        taxa_ids = taxa_ids.sample(x) if len(taxa_ids) > x else taxa_ids
        df = df[df['ncbiID'].isin(taxa_ids)].reset_index(drop=True)
    return df


def runtimeTest(func, **args):
    """ Function to calculate the runtime of a given function.
    :param function func: Function to perform the workload of the test.
    :param dict args: Takes keyword named arguments for the function func.
    :return: Tuple with (real_runtime, user_runtime, system_runtime)
    """

    # Save start timings
    start_time = time.perf_counter()

    # The resource feature is not available on Windows system
    win = os.name == "nt"
    if not win:
        start_sys_time = resource.getrusage(resource.RUSAGE_SELF)

    # Execute the given function with all given arguments.
    res = func(**args)

    # Save finish timings
    end_time = time.perf_counter()
    if not win:
        end_sys_time = resource.getrusage(resource.RUSAGE_SELF)

    # print(res) # Debug to check funtion results.

    # Give back real_runtime, user_runtime, sys_runtime to show efficiency.
    real_runtime = end_time - start_time

    if not win:
        user_runtime = end_sys_time.ru_utime - start_sys_time.ru_utime
        sys_runtime = end_sys_time.ru_stime - start_sys_time.ru_stime
        return real_runtime, user_runtime, sys_runtime
    else:
        return real_runtime, float("nan"), float("nan")


def getTasks():
    """ Create standard gene and taxa task. These are used by uGeneGUI.py as well.
    :return: Tuple of standard gene and taxa task and as well with low memory option.
    """
    gene_task = {'y_axis': 'geneID',
                 'x_axis': 'ncbiID',
                 'values': ['FAS_F', 'FAS_B'],
                 'jobs': [{'job_name': 'gene', 'n_components': 1},
                          {'job_name': 'gene', 'n_components': 2},
                          {'job_name': 'gene', 'n_components': 3}]}

    gene_task_low_mem = {'y_axis': 'geneID',
                         'x_axis': 'ncbiID',
                         'values': ['FAS_F', 'FAS_B'],
                         'jobs': [{'job_name': 'gene', 'n_components': 1, 'low_memory': True},
                                  {'job_name': 'gene', 'n_components': 2, 'low_memory': True},
                                  {'job_name': 'gene', 'n_components': 3, 'low_memory': True}]}

    taxa_task = {
        'x_axis': 'geneID',
        'y_axis': 'ncbiID',
        'values': ['FAS_F', 'FAS_B'],
        'jobs': [
            {'job_name': 'taxa', 'n_components': 1},
            {'job_name': 'taxa', 'n_components': 2},
            {'job_name': 'taxa', 'n_components': 3}
        ]}

    taxa_task_low_mem = {
        'x_axis': 'geneID',
        'y_axis': 'ncbiID',
        'values': ['FAS_F', 'FAS_B'],
        'jobs': [
            {'job_name': 'taxa', 'n_components': 1, 'low_memory': True},
            {'job_name': 'taxa', 'n_components': 2, 'low_memory': True},
            {'job_name': 'taxa', 'n_components': 3, 'low_memory': True}
        ]}

    return gene_task, gene_task_low_mem, taxa_task, taxa_task_low_mem


def runtimeConstGenes():
    # Create standard gene and taxa task. These are used by uGeneGUI.py as well.
    gene_task, gene_low_mem, taxa_task, taxa_low_mem = getTasks()

    # Dataframe to store the runtime results.
    try:
        df_stat = pd.read_csv("runtimeGene25000.csv")
    except Exception as error:
        df_stat = pd.DataFrame(columns=['setting', 'size', 'runtime'])
        print("Create new dataframe ! ", error)
    # Change the test setup by iterated workload list.
    for it_setting in ["vanilla", "low_memory"]:
        for it_xy in [10, 2500, 5000, 7500, 10000, 12500, 15000, 17500, 20000, 22500, 25000, 27500, 30000, 32500, 35000,
                      37500, 40000]:

            try:
                df = newSample(25000, it_xy)
                if it_setting == "low_memory":
                    real, user, system = runtimeTest(uGene.mainAnalytics, df=df, **gene_low_mem)
                else:
                    real, user, system = runtimeTest(uGene.mainAnalytics, df=df, **gene_task)

                print("Sample Done with ", it_xy, " in r:", real // 60, " u: ", user // 60, " s: ", system // 60)
                df_stat = pd.concat([df_stat, pd.DataFrame({'setting': [it_setting], 'size': [it_xy],
                                                            'runtime': [real // 60]})], ignore_index=True)

                time.sleep(2)
                df_stat.to_csv("runtimeGene25000.csv", index=False)
            except Exception as error:
                print("Runtime test Fail! ", error)

    # Create a line plot with the given runtimes
    fig = px.line(df_stat, x='size', y="runtime", color="setting", symbol="setting")

    # Manual update traces.
    for it_scatter in fig['data']:
        # Increase marker size.
        it_scatter['marker'].update({'size': 10})

    # fig.show()

    # Save the runtime plot to reopen and share the results.
    figuresToHtml([fig], filename="runtimeGene25000.html")


def runtimeConstTaxa():
    # Create standard gene and taxa task. These are used by uGeneGUI.py as well.
    gene_task, gene_low_mem, taxa_task, taxa_low_mem = getTasks()

    # Dataframe to store the runtime results.
    try:
        df_stat = pd.read_csv("runtimeTaxa25000.csv")
    except Exception as error:
        df_stat = pd.DataFrame(columns=['setting', 'size', 'runtime'])
        print("Create new dataframe ! ", error)
    # Change the test setup by iterated workload list.
    for it_setting in ["vanilla", "low_memory"]:
        for it_xy in [10, 2500, 5000, 7500, 10000, 12500, 15000, 17500, 20000, 22500, 25000, 27500, 30000, 32500, 35000,
                      37500, 40000]:

            try:
                df = newSample(it_xy, 25000)
                if it_setting == "low_memory":
                    real, user, system = runtimeTest(uGene.mainAnalytics, df=df, **gene_low_mem)
                else:
                    real, user, system = runtimeTest(uGene.mainAnalytics, df=df, **gene_task)

                print("Sample Done with ", it_xy, " in r:", real // 60, " u: ", user // 60, " s: ", system // 60)
                df_stat = pd.concat([df_stat, pd.DataFrame({'setting': [it_setting], 'size': [it_xy],
                                                            'runtime': [real // 60]})], ignore_index=True)

                time.sleep(2)
                df_stat.to_csv("runtimeTaxa25000.csv", index=False)
            except Exception as error:
                print("Runtime test Fail! ", error)

    # Create a line plot with the given runtimes
    fig = px.line(df_stat, x='size', y="runtime", color="setting", symbol="setting")

    # Manual update traces.
    for it_scatter in fig['data']:
        # Increase marker size.
        it_scatter['marker'].update({'size': 10})

    # fig.show()

    # Save the runtime plot to reopen and share the results.
    figuresToHtml([fig], filename="runtimeTaxa25000.html")


def runtimeConstTaxa():
    # Create standard gene and taxa task. These are used by uGeneGUI.py as well.
    gene_task, gene_low_mem, taxa_task, taxa_low_mem = getTasks()

    # Read real empiric data.
    df = pd.read_csv("PipelineData/F297.csv")

    # Dataframe to store the runtime results.
    try:
        df_stat = pd.read_csv("runtimeRealTaxa.csv")
    except Exception as error:
        df_stat = pd.DataFrame(columns=['setting', 'size', 'runtime'])
        print("Create new dataframe ! ", error)
    # Change the test setup by iterated workload list.
    for it_xy in [1312, 2624, 3936, 5248, 6560, 7872, 9184, 10496, 11808, 13120, 14432, 15744, 17056]:
        for repeats in range(2):
            df_sub = subSample(df, it_xy)
            for it_setting in ["vanilla", "low_memory"]:
                try:
                    if it_setting == "low_memory":
                        real, user, system = runtimeTest(uGene.mainAnalytics, df=df_sub, **gene_low_mem)
                    else:
                        real, user, system = runtimeTest(uGene.mainAnalytics, df=df_sub, **gene_task)

                    print("Sample Done with ", it_xy, " in r:", real, " u: ", user, " s: ",
                          system, " setting:", it_setting)

                    # Save results to df_stat.
                    df_stat = pd.concat([df_stat, pd.DataFrame({'setting': [it_setting], 'size': [it_xy],
                                                                'runtime': [real]})], ignore_index=True)

                    time.sleep(2)
                    df_stat.to_csv("runtimeRealTaxa.csv", index=False)
                except Exception as error:
                    print("Runtime test Fail! ", error)

    #print(df_stat)

    # Calculate mean an standard deviation.
    df_group = df_stat.groupby(['setting', 'size'])
    df_mean = df_group.mean().reset_index()
    df_std = df_group.std().reset_index()

    # Create a line plot with the given runtimes
    fig = px.line(df_mean, x='size', y="runtime", color="setting", symbol="setting")

    # Manual update traces.
    for it_scatter in fig['data']:

        # Keep in mind 'runtime' means the standard deviation of the runtime.
        it_scatter["error_y"] = dict(
            type='data',
            array=df_std[df_std['setting'] == it_scatter['legendgroup']]['runtime'],
            visible=True)

        # Increase marker size.
        it_scatter['marker'].update({'size': 10})

        # Update line color and shape.
        if it_scatter['legendgroup'] == "low_memory":
            it_scatter["line"] = dict(color='firebrick', width=3, dash='dot')
        else:
            it_scatter["line"] = dict(color='royalblue', width=3, dash='dash')

    # Save the runtime plot to reopen and share the results.
    figuresToHtml([fig], filename="runtimeRealTaxa.html")


def runtimeConstGene():
    # Create standard gene and taxa task. These are used by uGeneGUI.py as well.
    gene_task, gene_low_mem, taxa_task, taxa_low_mem = getTasks()

    # Read real empiric data.
    df = pd.read_csv("PipelineData/F297.csv")

    # Dataframe to store the runtime results.
    try:
        df_stat = pd.read_csv("runtimeRealGene.csv")
    except Exception as error:
        df_stat = pd.DataFrame(columns=['setting', 'size', 'runtime'])
        print("Create new dataframe ! ", error)
    # Change the test setup by iterated workload list.
    for it_xy in [13, 26, 39, 52, 65, 78, 91, 104, 117, 130, 143, 156, 169]:
        for repeats in range(2):
            df_sub = subSample(df, None, it_xy)
            for it_setting in ["vanilla", "low_memory"]:
                try:
                    if it_setting == "low_memory":
                        real, user, system = runtimeTest(uGene.mainAnalytics, df=df_sub, **gene_low_mem)
                    else:
                        real, user, system = runtimeTest(uGene.mainAnalytics, df=df_sub, **gene_task)

                    print("Sample Done with ", it_xy, " in r:", real, " u: ", user, " s: ",
                          system, " setting:", it_setting)

                    # Save results to df_stat.
                    df_stat = pd.concat([df_stat, pd.DataFrame({'setting': [it_setting], 'size': [it_xy],
                                                                'runtime': [real]})], ignore_index=True)

                    time.sleep(2)
                    df_stat.to_csv("runtimeRealGene.csv", index=False)
                except Exception as error:
                    print("Runtime test Fail! ", error)

    #print(df_stat)

    # Calculate mean an standard deviation.
    df_group = df_stat.groupby(['setting', 'size'])
    df_mean = df_group.mean().reset_index()
    df_std = df_group.std().reset_index()

    # Create a line plot with the given runtimes
    fig = px.line(df_mean, x='size', y="runtime", color="setting", symbol="setting")

    # Manual update traces.
    for it_scatter in fig['data']:

        # Keep in mind 'runtime' means the standard deviation of the runtime.
        it_scatter["error_y"] = dict(
            type='data',
            array=df_std[df_std['setting'] == it_scatter['legendgroup']]['runtime'],
            visible=True)

        # Increase marker size.
        it_scatter['marker'].update({'size': 10})

        # Update line color and shape.
        if it_scatter['legendgroup'] == "low_memory":
            it_scatter["line"] = dict(color='firebrick', width=3, dash='dot')
        else:
            it_scatter["line"] = dict(color='royalblue', width=3, dash='dash')

    # Save the runtime plot to reopen and share the results.
    figuresToHtml([fig], filename="runtimeRealGene.html")

if __name__ == "__main__":
    runtimeConstTaxa()
    runtimeConstGenes()
    runtimeConstTaxa()
    runtimeConstGenes()