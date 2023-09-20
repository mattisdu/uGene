"""
uGeneGUI.py
Version: 0.5
Beschreibung: This is the uGene graphical interface. It's built on Python-Dash and requires a browser for use..
              Example: $> python uGeneGUI.py

Autor: Mattis Kaumann

The MIT License (MIT)

Copyright (c) 2023 Mattis Kaumann, Goethe-Universität Frankfurt am Main
Read more on LICENSE.txt
"""

import numpy as np
import math
import os
import socket
import json

import dash
import dash_bootstrap_components as dbc
import pandas as pd
import tkinter as tk
import tkinter.filedialog as fd

import plotly.graph_objs as go
import plotly.express as px
import subprocess
import time
import concurrent.futures
import scipy.stats as sst

import uGeneCore as uGene


# ------------------------------------------ Function declarations -----------------------------------------------------
def rgbStrToVec(color):
    """ Converts a hex color code string into a numpy 3 vector.
    :param color: Color code string. An example would be "#1A05FF".
    :return: Returns a numpy 3 vector with red green and blue value.
    """
    try:
        return np.array([int("0x" + color[1:3], 16),
                         int("0x" + color[3:5], 16),
                         int("0x" + color[5:7], 16)])
    except ... as error:
        print("Error: required_functionalities->rgbStrToVec():", error)
        return np.array([0, 0, 0])


def rgbVecToStr(c_vec):
    """ Converts a numpy 3 vector within int variables into a rbg hex color
    code string.
    :param c_vec: Numpy 3 vector within int variables between 0 and 255.
    :return: Returns a string with a hex color code.
    """
    try:
        # Ensure for all components to be in range 0 up to 255.
        c_vec[0] = max(0, min(c_vec[0], 255))
        c_vec[1] = max(0, min(c_vec[1], 255))
        c_vec[2] = max(0, min(c_vec[2], 255))

        return "#" + str(hex(c_vec[0]))[2:4].zfill(2) + \
               str(hex(c_vec[1]))[2:4].zfill(2) + \
               str(hex(c_vec[2]))[2:4].zfill(2)
    except ... as error:
        print("Error: required_functionalities->rgbVecToStr()", error)
        return "#000000"


def colorRampPalette(colors, n):
    """ Interpolate colors linearly to create a color palette.
    :param colors:  List with color hex strings which is based on.
    :param n:       Number of required colors. That effects the greatness of return list.
    :return:        Gives a list with hex color strings.
    """
    result = []
    c_len = len(colors)
    if c_len < 1:
        return []
    if c_len == 1:
        return colors * n
    if n == 1:
        return [colors[0]]

    step = (len(colors) - 1) / (n - 1)
    for i in range(0, n):
        if math.floor(step * i) == math.ceil(step * i):
            result.append(colors[math.floor(step * i)])
        else:
            v_color_a = rgbStrToVec(colors[math.floor(step * i)])
            v_color_b = rgbStrToVec(colors[math.ceil(step * i)])

            v_color = (v_color_a + (v_color_b - v_color_a) *
                       (step * i % 1)).astype(int)
            result.append(rgbVecToStr(v_color))

    return result


def scatterDataExpress(x_2d=None, y_2d=None, x_3d=None, y_3d=None, z_3d=None, color=None, customdata=None, df=None,
                       dot_size=5, rev_size=0.5, opacity=0.9, border=None):
    """ Create three plots with in just one trace and the given custom options.
    All params until df be able to call a column name from dataframe.
    :param x_2d: List of x values for the 2-dimensional plot or df column name.
    :param y_2d: List of y values for the 2-dimensional plot or df column name.
    :param x_3d: List of x values for the 3-dimensional plot or df column name.
    :param y_3d: List of y values for the 3-dimensional plot or df column name.
    :param z_3d: List of z values for the 3-dimensional plot or df column name.
    :param color: List of colorcode strings supported by plotly or df column name which holds such list.
    :param customdata: List of objects excluding listed strings. These will be interpreted as list of column names.
    :param df: Pandas dataframe with any kind of plotting data.
    :param dot_size: Int for scatterplot dot size.
    :param rev_size: Float for sizing dots into scatterplot. For example legend dot size. View plotly doku.
    :param opacity: Float between 1.0 and 0.0 witch defines the dot opacity.
    :param border: List of colorcode strings supported by plotly or df column name which holds such list.
    :return: Returns a triplet of plotly-go-scatter data objects.
    """

    # Process column names argument handovers.
    if df is not None:
        if type(x_2d) == str and x_2d in df.columns:
            x_2d = list(df[x_2d])
        if type(y_2d) == str and y_2d in df.columns:
            y_2d = list(df[y_2d])
        if type(x_3d) == str and x_3d in df.columns:
            x_3d = list(df[x_3d])
        if type(y_3d) == str and y_3d in df.columns:
            y_3d = list(df[y_3d])
        if type(z_3d) == str and z_3d in df.columns:
            z_3d = list(df[z_3d])
        if type(color) == str:
            if color in df.columns:
                color = list(df[color])
            else:
                # Allow to use color string
                color = len(x_3d) * [color]
        if border and type(border) == str:
            if border in df.columns:
                border = list(df[border])
            else:
                # Allow to use color string
                border = len(x_3d) * [border]

        if type(customdata) == str:
            customdata = [customdata]
        if type(customdata) == list and all([type(it) == str and it in df.columns for it in customdata]):
            customdata = df[customdata].values.tolist()

    # Check if all data fit to each other. List of data have to hold the same amount of related data.
    if not x_2d or not y_2d or not x_3d or not y_3d or not z_3d or not color or not customdata:
        print("Error scatterExpress()! some arguments are none.")
        return []
    if not len(x_2d) == len(y_2d) == len(x_3d) == len(y_3d) == len(z_3d) == len(color) == len(customdata):
        print("Error scatterExpress()! given data do not have consistent size.")
        return []

    # The border exists only if set. Otherwise, none.
    if border and len(border) != len(x_2d):
        print("Error scatterExpress()! given data do not have consistent size.")
        return []

    # Create one trace data storage solution with go.Scattergl.
    data_2d = [go.Scattergl(
        customdata=np.array(customdata),
        hovertemplate="%{customdata[0]}<extra></extra>",
        x=x_2d,
        y=y_2d,
        name="main_trace",
        showlegend=False,
        mode='markers',
        marker=dict(
            symbol='circle',
            opacity=float(opacity),
            size=float(dot_size),
            sizemode='area',
            sizeref=float(rev_size),
            color=color,
            line=dict(
                width=dot_size * 0.25,
                color=border if border else [app_con.const.color_limpid] * len(x_2d)
            ),
        )
    )]

    # Create with go.Scatter3d a one trace data storage solution. Scatter3d is native supported by web gl.
    data_3d = [go.Scatter3d(
        customdata=np.array(customdata),
        hovertemplate="%{customdata[0]}<extra></extra>",
        x=x_3d,
        y=y_3d,
        z=z_3d,
        name="main_trace",
        showlegend=False,
        mode='markers',
        marker=dict(
            symbol='circle',
            opacity=float(opacity),
            size=float(dot_size) * 0.3,
            sizemode='area',
            sizeref=float(rev_size),
            color=color,
            line=dict(
                width=dot_size * 0.25,
                color=border if border else [app_con.const.color_limpid] * len(x_2d)
            ),
        )
    )]

    # Create with go.Splom a one trace data storage solution. Splom is native supported by web gl.
    data_matrix_3d = [go.Splom(
        customdata=np.array(customdata),
        hovertemplate="%{customdata[0]}<extra></extra>",
        dimensions=[{'axis': {'matches': True}, 'label': 'gene3d_x', 'values': x_3d},
                    {'axis': {'matches': True}, 'label': 'gene3d_y', 'values': y_3d},
                    {'axis': {'matches': True}, 'label': 'gene3d_z', 'values': z_3d}],
        name="main_trace",
        showlegend=False,
        marker=dict(
            symbol='circle',
            opacity=float(opacity),
            size=float(dot_size),
            sizemode='area',
            sizeref=float(rev_size),
            color=color,
            line=dict(
                width=dot_size * 0.25,
                color=border if border else [app_con.const.color_limpid] * len(x_2d)
            ),
        )
    )]

    return data_2d, data_3d, data_matrix_3d


def processOptionDf(df, data_sel, opt_arrange, opt_subset):
    """ Process dataframe reduction and arrangement based on the user selected options.
    :param df: Dataframe within all current phylogenetic data and as well all needed cluster data, for example gene1d_x.
    :param data_sel: A dictionary with data of selected genes as well there group name. Gene names are main keys.
    :param opt_arrange: String of the current selected option how to order genes. One of app_con.const.gene_arrange
    :param opt_subset: String of the current selected gene subset setting. One of app_con.const.gene_subset
    :return: Returns the dataframe with the applied options.
    """

    if opt_subset == "group-based" and data_sel:
        print("Waring processOptionDf() fail. No groups are selected.")
        # Filter by gene id´s form selection dict.
        df = df[df['geneID'].isin(list(data_sel.keys()))]

    # Just warn and continue with the full set to keep these code running.
    if "gene1d_x" not in df.columns:
        print("Waring processOptionDf() fail. Column gene1d_x should be available.")

    # Process the gene rearrangement.
    if opt_arrange == "1d-order" and "gene1d_x" in df.columns:
        df = df.sort_values('gene1d_x').reset_index(drop=True)
    elif opt_arrange == "group-based" or opt_arrange == "1d-group-based":
        groups = {}

        # All not grouped genes get ordered by origin or by 1d order. Only the grouped genes will put into first place.
        if opt_arrange == "1d-group-based" and "gene1d_x" in df.columns:
            gene_order = list(df.sort_values('gene1d_x')["geneID"].drop_duplicates())
        else:
            gene_order = list(df["geneID"].drop_duplicates())

        for it in data_sel.keys():
            gene_order.remove(it)
            if data_sel[it]['name'] in groups:
                groups[data_sel[it]['name']].append(it)
            else:
                groups[data_sel[it]['name']] = [it]

        # New gene oder is get finished. The list gene_order holds the distinct new order of genes.
        gene_order = sum([groups[itt] for itt in groups], []) + gene_order

        df['geneID'] = df['geneID'].astype('category')
        df['geneID'] = df['geneID'].cat.set_categories(gene_order)

        df = df.sort_values('geneID').reset_index(drop=True)

    # Just warn but give a dataframe back in any way.
    return df


def askSaveHelp(title, filetypes):
    """Help function to start tkinter save dialog into an own process. This avoids a thread based bug.
    :param title: Dialog box title string
    :param filetypes: Tuple with in pairs. One pair have to be (description, file_extension)
    :return: Path string to the selected location.
    """
    tk_root = tk.Tk()
    tk_root.withdraw()
    path = fd.asksaveasfilename(title=title, filetypes=filetypes)
    tk_root.destroy()

    return path


def askSave(title, filetypes):
    """Function to start tkinter save dialog into an own process.  This avoids a thread based bug.
        :param title: Dialog box title string
        :param filetypes: Tuple with in pairs. One pair have to be (description, file_extension)
        :return: Path string to the selected location."""
    try:

        tk_root = tk.Tk()
        tk_root.withdraw()
        path = fd.asksaveasfilename(title=title, filetypes=filetypes)
        tk_root.destroy()

        return path
    except Exception as error:
        print("Warning solve open dialog.", error)
        with concurrent.futures.ProcessPoolExecutor() as executor:
            process = executor.submit(askSaveHelp, title, filetypes)
        return process.result()


def askOpenHelp(title, filetypes):
    """Help function to start tkinter open file dialog into an own process. This avoids a thread based bug.
        :param title: Dialog box title string
        :param filetypes: Tuple with in pairs. One pair have to be (description, file_extension)
        :return: Path string to the selected location.
    """
    tk_root = tk.Tk()
    tk_root.withdraw()

    path = fd.askopenfilename(title=title, filetypes=filetypes)

    tk_root.destroy()

    return path


def askOpen(title, filetypes):
    """Function to start tkinter open file dialog into an own process. This avoids a thread based bug.
        :param title: Dialog box title string
        :param filetypes: Tuple with in pairs. One pair have to be (description, file_extension)
        :return: Path string to the selected location.
    """
    try:
        tk_root = tk.Tk()
        tk_root.withdraw()
        path = fd.askopenfilename(title=title, filetypes=filetypes)
        tk_root.destroy()

        return path
    except Exception as error:
        print("Warning solve open dialog.", error)
        with concurrent.futures.ProcessPoolExecutor() as executor:
            process = executor.submit(askOpenHelp, title, filetypes)
        return process.result()


def processAdvancedOptions(gene_task, taxa_task, str_adv_opt, task_arguments=[]):
    # Clean from whitespaces and prepare to split. Every new line, tab and ; become a new command.
    str_adv_opt = str_adv_opt.replace(" ", "").replace("\t", ";").replace("\n", ";")

    # First of the tuple hold task options and second keeps all job options.
    adv_opt = {'opt_general': ({}, {}), 'opt_gene': ({}, {}), 'opt_taxa': ({}, {})}
    p_opt = 'opt_general'

    # More or less parse the advance options input.
    for option in str_adv_opt.split(";"):
        if option[:len("gene:")] == "gene:":
            p_opt = 'opt_gene'
            option = option[len("gene:"):]
        if option[:len("taxa:")] == "taxa:":
            p_opt = 'opt_taxa'
            option = option[len("taxa:"):]

        if not option:
            continue

        option = option.split("=")
        if len(option) != 2 or len(option[0]) < 1 or len(option[1]) < 1:
            print("Error runCluster()! Can´t parse advance options.")
            print(option)
        else:
            try:
                if option[0] in task_arguments:
                    adv_opt[p_opt][0][str(option[0])] = json.loads(option[1].replace("'", '"'))
                else:
                    adv_opt[p_opt][1][str(option[0])] = json.loads(option[1].replace("'", '"'))
            except Exception as error:
                print("Error runCluster()! Can´t interpret value by json.loads().", error)
                print(option)

    # Update the 'gene_task'. Remember, the options for a task and the options for a job are different things.
    gene_task.update(adv_opt['opt_general'][0])
    gene_task.update(adv_opt['opt_gene'][0])
    for it in gene_task['jobs']:
        it.update(adv_opt['opt_general'][1])
        it.update(adv_opt['opt_gene'][1])
    print("Updated gene task:")
    print(gene_task)

    # Update the 'taxa_task'. Remember, the options for a task and the options for a job are different things.
    if taxa_task:
        taxa_task.update(adv_opt['opt_general'][0])
        taxa_task.update(adv_opt['opt_taxa'][0])
        for it in gene_task['jobs']:
            it.update(adv_opt['opt_general'][1])
            it.update(adv_opt['opt_taxa'][1])
        print("Updated taxa task:")
        print(taxa_task)

    return gene_task, taxa_task


def divideHelp(res, x):
    """ Help funktion which it used into updateAdditionalData(). We like to get res and res/x without calculate res a
    second time.
    :param res: Any float
    :param x: Int or float
    :return: Tuple with res and res divide by x
    """
    return res, res / x


def freePort():
    """ Function to retrieve an available port. This ensures that PhyloProfile launches on an unoccupied port.
    :return: Int representing the currently available port
    """
    with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as s:
        # Bind to port 0 will return a free port
        s.bind(('127.0.0.1', 0))
        ip, port = s.getsockname()
    return port


# -------------------------------------------- Class declarations ------------------------------------------------------
class DI(dict):
    """Better dictionary class to access by d.value_name"""
    __getattr__ = dict.get
    __setattr__ = dict.__setitem__
    __delattr__ = dict.__delitem__


# --------------------------------------- Global app container definition ----------------------------------------------
app_con = DI(
    ids=DI(

        coll_main_header="collapse_main_header",
        tabs_ugene_viewer="main_tabs_ugene_viewer",
        tab_ugene="main_tab_ugene",
        tab_viewer="main_tab_viewer",
        tab_download="main_tab_download",

        tabs_show_data="tabs_show_data",
        tab_show_phyloprofile="tab_show_phyloprofile",

        input_int_nn="input_int_n_neighbors",
        input_float_md="input_float_min_dist",
        dd_use_metric="dd_use_metric",
        bt_origin_order="bt_origin_order",
        store_order_id="store_order_id_main",
        bt_adv_opt="bt_advanced_opt",
        coll_adv_opt="coll_advanced_opt",
        text_adv_opt="text_advanced_opt",
        row_adv_label="row_advanced_label",

        bt_load_csv="bt_load_csv",
        bt_show_result="bt_show_result",
        bt_load_cluster="bt_load_cluster",
        bt_load_groups="bt_load_groups",

        spin_load_csv="spin_load_csv",
        spin_show_result="spin_show_results",
        spin_load_cluster="spin_load_cluster",
        spin_load_groups="spin_load_groups",

        user_location="user_location",
        input_csv="main_input_csv",
        input_cluster="main_input_cluster",
        filename_csv="filename_csv_file",
        filename_cluster="filename_cluster_file",
        df_selection="df_selection",
        plot_2d="main_plot_2d",
        plot_3d="main_plot_3d",
        scatter_3d="scatter_matrix_3d",

        column_gene_viewer="column_gene_viewer",
        column_stat="column_statistics",
        coll_stat="collapse_statistics",
        store_coll_stat="store_coll_stat",
        store_sel_stat="store_selection_stat",
        bt_unfold_stat="bt_unfold_statistics",
        bt_load_stat="bt_load_statistics",
        bar_plot_aaa="bar_plot_aaa",
        bar_plot_bbb="bar_plot_bbb",
        bar_plot_ccc="bar_plot_ccc",
        coll_bar_ccc="coll_bar_ccc",
        bar_plot_ddd="bar_plot_ddd",
        coll_bar_ddd="coll_bar_ddd",
        bar_plot_eee="bar_plot_eee",
        coll_bar_eee="coll_bar_eee",
        bar_plot_fff="bar_plot_fff",
        coll_bar_fff="coll_bar_fff",
        bar_plot_ggg="bar_plot_ggg",
        coll_bar_ggg="coll_bar_ggg",
        bar_plot_hhh="bar_plot_hhh",
        coll_bar_hhh="coll_bar_hhh",
        bar_plot_iii="bar_plot_iii",
        coll_bar_iii="coll_bar_iii",
        bar_plot_jjj="bar_plot_jjj",
        coll_bar_jjj="coll_bar_jjj",
        input_int_bar_n_best="input_int_bar_n_best",
        filename_stat="filename_stat_file",
        spin_load_stat="spin_load_stat",
        bt_hypergeometirc_stat="bt_hypergeometirc_stat",
        bt_bonferroni_cor="bt_bonferroni_correction",

        coll_options="collapse_option_bar",
        dd_color_pallette="dd_color_palette",
        active_color_palette="active_color_palette",
        slider_dot_size="slider_dot_size",
        slider_opacity="slider_opacity",
        card_phyloprofile="main_phyloprofile",
        phylo_iframe="phylo_iframe",

        di_ng="di_new_group",
        di_ng_name="di_ng_name",
        di_ng_color="di_ng_color",
        di_ng_cancel="di_ng_cancel",
        di_ng_create="di_ng_create",

        bt_di_ng="bt_new_group",
        bt_delete_group="bt_delete_group",
        dd_active_group="dd_active_group",
        dd_modify_tool="dd_modify_tool",

        bt_add_list="bt_add_list",
        di_add_list="di_add_list",
        di_add_list_input="di_add_list_input",
        di_add_list_cancel="di_add_list_cancel",
        di_add_list_apply="di_add_list_apply",

        sw_cluster_dir="sw_cluster_dir",

        deleted_group="store_deleted_group",
        dict_selection="store_selection",
        table_selection="main_data_table",

        dd_live_arrange="dd_live_arrange",
        dd_live_subset="dd_live_subset",
        bt_update_pyhlo="bt_update_pyhlo",
        coll_pyhlo_opt="coll_pyhlo_opt",

        debug_button="debug_button",
        dummy_restyle="dummy_restyle",
        dummy_message="dummy_message",
        user_message="user_message",

        bt_down_cluster="bt_down_cluster",
        bt_down_phylo="bt_down_phylo",
        bt_down_groups="bt_down_groups",
        bt_down_cat="bt_down_cat",
        dd_opt_subset="dd_opt_subset",
        dd_opt_arrange="dd_opt_arrange",

        spin_down_cluster="spin_down_cluster",
        spin_down_phylo="spin_down_phylo",
        spin_down_groups="spin_down_groups",
        spin_down_cat="spin_down_cat"

    ),
    data=DI(
        # Data like the main profiles will be stored here. File path of loaded data will become the dict keys.
    ),
    color_palettes=DI(
        Rainbow=['#DD0000', '#FFFF00', '#309000', '#00FF00', '#00DDDD', '#0202BB', '#FF80CC'],
        Black=['#000000', '#000000'],
        Grey=['#7F7F7F', '#7F7F7F'],
        Colorcycle=['#332288', '#88CCEE', '#44AA99', '#117733', '#999933', '#DDCC77', '#CC6677', '#882255', '#AA4499'],
        # The following color schematics are sourced online to enhance accessibility for color-blind individuals.
        # Author: Paul Tol
        # Email: p.j.j.tol@sron.nl
        # Title: "Introduction to Colour Schemes"
        # Accessed Date: 31.08.2023
        # Source URL: https://personal.sron.nl/~pault/#fig:scheme_bright
        Bright_Colorbilnd=['#4477AA', '#EE6677', '#228833', '#CCBB44', '#66CCEE', '#AA3377', '#BBBBBB'],
        Vibrant_Colorblind=['#EE7733', '#0077BB', '#33BBEE', '#EE3377', '#CC3311', '#009988', '#BBBBBB'],
        Dark_Colorblind=['#CC6677', '#332288', '#DDCC77', '#117733', '#88CCEE', '#882255', '#44AA99', '#999933',
                         '#AA4499'],
        Bight_Compromise=['#77AADD', '#EE8866', '#EEDD88', '#FFAABB', '#99DDFF', '#44BB99', '#BBCC33', '#AAAA00',
                          '#DDDDDD']
    ),
    const=DI(
        ugene_gui_version="0.5",
        color_limpid="rgba(235,235,250,0.0)",
        color_phylo_defauild="#FFFFFF",
        phylo_defaild_cat="__none__cat__",
        temp_folder="/uGeneTemp/",
        gene_col="geneID",
        taxa_col="ncbiID",
        phylo_col=['geneID', 'ncbiID', 'orthoID', 'FAS_F', 'FAS_B'],
        phylo_cat_col=['geneID', 'name', 'color'],
        umap_metrics=['euclidean', 'manhattan', 'chebyshev', 'minkowski'],
        cluster_col_names=['gene1d_x', 'gene2d_x', 'gene2d_y', 'gene3d_x', 'gene3d_y', 'gene3d_z'],
        gene_arrange=['1d-group-based', '1d-order', 'origin', 'group-based'],
        gene_subset=['full', 'group-based'],
        task_arguments=['dev_report', 'fill', 'pattern', 'drop', 'x_axis', 'y_axis', 'values']
    )
)

# ------------------------------------------- Create main dash app -----------------------------------------------------
app = dash.Dash(__name__, external_stylesheets=[dbc.themes.BOOTSTRAP, dbc.icons.FONT_AWESOME])
app.layout = dbc.Container(fluid=True, children=[
    dash.dcc.Location(id=app_con.ids.user_location),

    dash.dcc.Store(id=app_con.ids.filename_csv, storage_type='memory'),
    dash.dcc.Store(id=app_con.ids.filename_cluster, storage_type='memory'),
    dash.dcc.Store(id=app_con.ids.filename_stat, storage_type='memory'),

    dash.dcc.Store(id=app_con.ids.df_selection, storage_type='memory'),

    dash.dcc.Store(id=app_con.ids.active_color_palette, storage_type='memory'),
    dash.dcc.Store(id=app_con.ids.store_order_id, storage_type='memory'),

    dash.dcc.Store(id=app_con.ids.deleted_group, storage_type='memory', data=dict()),
    dash.dcc.Store(id=app_con.ids.dict_selection, storage_type='memory', data=dict()),

    dash.dcc.Store(id=app_con.ids.store_coll_stat, data=False),
    dash.dcc.Store(id=app_con.ids.store_sel_stat),

    # Design display true content.
    dbc.Collapse(id=app_con.ids.coll_main_header, is_open=True, children=[
        dbc.Card(color="secondary", inverse=True, children=dbc.CardHeader([
            dash.html.H1("uGene Dashboard"),
            dash.html.Div(className="position-absolute top-0 end-0", children=[
                dbc.Label("Version: " + str(app_con.const.ugene_gui_version) + "/" + str(uGene.UGENE_CORE_VERSION)),
            ])
        ])),
        dash.html.Br()
    ]),
    dbc.Tabs(id=app_con.ids.tabs_ugene_viewer, children=[
        # File selection and analyse tools.
        dbc.Tab(tab_id=app_con.ids.tab_ugene, label='uGene tool', children=[
            dbc.Row([
                dbc.Col([
                    dbc.Card([
                        dbc.CardHeader([
                            dash.html.H3("Create uMap cluster"),
                        ]),
                        dbc.CardBody([
                            dbc.Row([
                                dash.html.P(
                                    "Start your workflow here. Load a phyloprofile, do your clustering and at least " +
                                    "view the results into the gene viewer tab."
                                ),
                            ]),
                            dash.html.Hr(),
                            dbc.Row([
                                dbc.Col(dash.html.Label("Adjust n_neighbors"), width=3),
                                dbc.Col(
                                    dbc.Input(
                                        type="number",
                                        value=15,
                                        step="1",
                                        id=app_con.ids.input_int_nn
                                    ),
                                    width=9
                                ),
                            ]),
                            dbc.Row([
                                dbc.Col(dash.html.Label("Adjust min_dist"), width=3),
                                dbc.Col(
                                    dbc.Input(
                                        type="number",
                                        value=0.1,
                                        step="0.001",
                                        id=app_con.ids.input_float_md
                                    ),
                                    width=9
                                ),
                            ]),
                            dbc.Row([
                                dbc.Col(dash.html.Label("Select metric"), width=3),
                                dbc.Col(width=9, children=[
                                    dash.dcc.Dropdown(
                                        id=app_con.ids.dd_use_metric,
                                        clearable=False,
                                        options=[{'label': it, 'value': it} for it in app_con.const.umap_metrics],
                                        value=app_con.const.umap_metrics[0]
                                    )
                                ]),
                            ]),
                            dash.html.Br(),
                            dbc.Row(id=app_con.ids.row_adv_label, children=[
                                dbc.Col(dash.html.Label("Use advanced options."), width=3),
                                dbc.Col(dbc.Checkbox(id=app_con.ids.bt_adv_opt, value=False))
                            ]),
                            dbc.Tooltip(
                                "Here you can override all uMap parameters. Write Parameter = Value on each line. ",
                                "For more information on all parameters, please refer to the uGene documentation.",
                                target=app_con.ids.row_adv_label
                            ),
                            dbc.Row([
                                dbc.Collapse(id=app_con.ids.coll_adv_opt, is_open=False, children=[
                                    dbc.Textarea(id=app_con.ids.text_adv_opt, className="mb-5",
                                                 style={'height': '15vh'}, placeholder="Add here your adjustments"),
                                ])
                            ]),
                            dash.html.Hr(),
                            dbc.Row([
                                dbc.Col(children=[
                                    dash.html.Div(className="hstack gap-2", children=[
                                        dbc.Button("Select File", id=app_con.ids.bt_load_csv, color="success"),
                                        dbc.Spinner(dash.html.Div(id=app_con.ids.spin_load_csv, children="...."),
                                                    size="sm")
                                    ])
                                ]),
                                dbc.Col(children=[
                                    dash.html.Div(className="hstack gap-2", children=[
                                        dbc.Button("Start uGene", id=app_con.ids.bt_show_result, color="secondary"),
                                        dbc.Spinner(dash.html.Div(id=app_con.ids.spin_show_result, children="...."),
                                                    size="sm")
                                    ])
                                ])
                            ])
                        ])
                    ])
                ]),
                dbc.Col([
                    dbc.Card([
                        dbc.CardHeader([
                            dash.html.H3("Start over")
                        ]),
                        dbc.CardBody([
                            dbc.Row([
                                dash.html.P(
                                    "Here it is possible to load a already clustered file to view the clusters. " +
                                    "Only .cluster.csv files are supported."
                                ),
                            ]),
                            dash.html.Hr(),
                            dbc.Row([
                                dash.html.Div(className="hstack gap-2", children=[
                                    dbc.Button("Select File", id=app_con.ids.bt_load_cluster, color="success"),
                                    dbc.Spinner(dash.html.Div(id=app_con.ids.spin_load_cluster, children="...."),
                                                size="sm")
                                ])
                            ])
                        ])
                    ]),
                    dash.html.Br(),
                    dbc.Card([
                        dbc.CardHeader([
                            dash.html.H3("Load Previous Groups")
                        ]),
                        dbc.CardBody([
                            dbc.Row([
                                dash.html.P(
                                    "Here, you have the option to import your pre-defined handcrafted groups. " +
                                    "Please note that only CSV files of categories downloaded from uGene are supported."
                                )
                            ]),
                            dash.html.Hr(),
                            dbc.Row([
                                dash.html.Div(className="hstack gap-2", children=[
                                    dbc.Button("Select File", id=app_con.ids.bt_load_groups, color="success"),
                                    dbc.Spinner(dash.html.Div(id=app_con.ids.spin_load_groups, children="...."),
                                                size="sm")
                                ])
                            ])
                        ])
                    ])
                ])
            ])
        ]),
        # Display uGene Results
        dbc.Tab(tab_id=app_con.ids.tab_viewer, label='uGene viewer', children=[
            dbc.Row(className="hstack", children=[
                dash.html.Div(id=app_con.ids.column_gene_viewer,
                              style={'width': '95vw', "pattingInline": "0,5vw", 'flexGrow': 0, 'flexShrink': 0,
                                     'flexBasis': 'auto'},
                              children=[
                                  dbc.Tabs(id=app_con.ids.tabs_show_data, children=[
                                      dbc.Tab(label='2D View', children=[
                                          dbc.Card(children=[
                                              dbc.CardBody(children=[
                                                  dash.dcc.Graph(
                                                      style={'height': '70vh'},
                                                      className="plot",
                                                      id=app_con.ids.plot_2d,
                                                      config={"displayModeBar": True}
                                                  )
                                              ])
                                          ])
                                      ]),
                                      dbc.Tab(label='3D View', children=[
                                          dbc.Card(children=[
                                              dbc.CardBody(children=[
                                                  dash.dcc.Graph(
                                                      style={'height': '70vh'},
                                                      id=app_con.ids.plot_3d,
                                                      config={"displayModeBar": True},
                                                      className="plot"
                                                  )
                                              ])
                                          ])
                                      ]),
                                      dbc.Tab(label='3D Scatter', children=[
                                          dbc.Card(children=[
                                              dbc.CardBody(children=[
                                                  dash.dcc.Graph(
                                                      style={'height': '70vh'},
                                                      id=app_con.ids.scatter_3d,
                                                      config={"displayModeBar": True},
                                                      className="plot"
                                                  )
                                              ])
                                          ])
                                      ]),
                                      dbc.Tab(tab_id=app_con.ids.tab_show_phyloprofile, label='Phyloprofile', children=[
                                          dbc.Card(children=[
                                              dbc.CardBody(style={'height': '70vh'}, children=[
                                                  dash.html.Iframe(style={'height': '90%', 'width': '100%'},
                                                                   id=app_con.ids.phylo_iframe,
                                                                   srcDoc='<p>Load Phyloprofile ...</p>',
                                                                   src='about:blank')
                                              ])
                                          ])
                                      ]),
                                  ])
                              ]),
                # Additional statistic data.
                dash.html.Div(
                    id=app_con.ids.column_stat,
                    style={'width': '3vw',
                           "paddingInline": "0.5vw",
                           'flexGrow': 0,
                           'flexShrink': 0,
                           'flexBasis': 'auto'
                           },
                    children=[
                        dash.html.Br(),
                        dbc.Card(className="hstack", style={'flexShrink': 1}, children=[
                            dash.html.Div(style={'height': '73vh', 'width': '2vw'}, children=[
                                dbc.Button(">", id=app_con.ids.bt_unfold_stat, color="primary",
                                           style={'width': '100%', 'height': "100%"})
                            ]),
                            dbc.Collapse(id=app_con.ids.coll_stat, dimension="width", is_open=False, children=[
                                dbc.CardBody(style={"width": "25vw", "padding": "0.5vw", 'height': '73vh',
                                                    'overflowY': 'scroll'}, children=[
                                    dbc.Row([
                                        dbc.Col(dash.html.H4("Additional Statistics")),
                                        dbc.Col(children=[
                                            dash.html.Div(className="hstack gap-2", children=[
                                                dbc.Button("Load Data", id=app_con.ids.bt_load_stat, color="primary"),
                                                dbc.Spinner(
                                                    dash.html.Div(id=app_con.ids.spin_load_stat, children="...."),
                                                    size="sm")
                                            ])
                                        ])
                                    ]),
                                    dash.html.Br(),
                                    dbc.Row([
                                        dash.dcc.Graph(
                                            style={'height': '35vh'},
                                            id=app_con.ids.bar_plot_aaa,
                                            config={"displayModeBar": True},
                                            className="plot"
                                        )
                                    ]),
                                    dbc.Row([
                                        dash.dcc.Graph(
                                            style={'height': '35vh'},
                                            id=app_con.ids.bar_plot_bbb,
                                            config={"displayModeBar": True},
                                            className="plot"
                                        )
                                    ]),
                                    dbc.Collapse(id=app_con.ids.coll_bar_ccc, is_open=False, children=[
                                        dbc.Row([
                                            dash.dcc.Graph(
                                                style={'height': '35vh'},
                                                id=app_con.ids.bar_plot_ccc,
                                                config={"displayModeBar": True},
                                                className="plot"
                                            )
                                        ])
                                    ]),
                                    dbc.Collapse(id=app_con.ids.coll_bar_ddd, is_open=False, children=[
                                        dbc.Row([
                                            dash.dcc.Graph(
                                                style={'height': '35vh'},
                                                id=app_con.ids.bar_plot_ddd,
                                                config={"displayModeBar": True},
                                                className="plot"
                                            )
                                        ])
                                    ]),
                                    dbc.Collapse(id=app_con.ids.coll_bar_eee, is_open=False, children=[
                                        dbc.Row([
                                            dash.dcc.Graph(
                                                style={'height': '35vh'},
                                                id=app_con.ids.bar_plot_eee,
                                                config={"displayModeBar": True},
                                                className="plot"
                                            )
                                        ])
                                    ]),
                                    dbc.Collapse(id=app_con.ids.coll_bar_fff, is_open=False, children=[
                                        dbc.Row([
                                            dash.dcc.Graph(
                                                style={'height': '35vh'},
                                                id=app_con.ids.bar_plot_fff,
                                                config={"displayModeBar": True},
                                                className="plot"
                                            )
                                        ])
                                    ]),
                                    dbc.Collapse(id=app_con.ids.coll_bar_ggg, is_open=False, children=[
                                        dbc.Row([
                                            dash.dcc.Graph(
                                                style={'height': '35vh'},
                                                id=app_con.ids.bar_plot_ggg,
                                                config={"displayModeBar": True},
                                                className="plot"
                                            )
                                        ])
                                    ]),
                                    dbc.Collapse(id=app_con.ids.coll_bar_hhh, is_open=False, children=[
                                        dbc.Row([
                                            dash.dcc.Graph(
                                                style={'height': '35vh'},
                                                id=app_con.ids.bar_plot_hhh,
                                                config={"displayModeBar": True},
                                                className="plot"
                                            )
                                        ])
                                    ]),
                                    dbc.Collapse(id=app_con.ids.coll_bar_iii, is_open=False, children=[
                                        dbc.Row([
                                            dash.dcc.Graph(
                                                style={'height': '35vh'},
                                                id=app_con.ids.bar_plot_iii,
                                                config={"displayModeBar": True},
                                                className="plot"
                                            )
                                        ])
                                    ]),
                                    dbc.Collapse(id=app_con.ids.coll_bar_jjj, is_open=False, children=[
                                        dbc.Row([
                                            dash.dcc.Graph(
                                                style={'height': '35vh'},
                                                id=app_con.ids.bar_plot_jjj,
                                                config={"displayModeBar": True},
                                                className="plot"
                                            )
                                        ])
                                    ]),
                                    dbc.Row([
                                        dbc.Col(dash.html.Label("Setup number of displayed bars:"), width=9),
                                        dbc.Col(dbc.Input(type="number", value=5, step="1",
                                                          id=app_con.ids.input_int_bar_n_best), width=3)
                                    ]),
                                    dbc.Row([
                                        dbc.Col(dash.html.Label("Use hyper-geometric statistic:"), width=9),
                                        dbc.Col(dbc.Checkbox(id=app_con.ids.bt_hypergeometirc_stat, value=False),
                                                width=3)
                                    ]),
                                    dbc.Row([
                                        dbc.Col(dash.html.Label("Use Bonferroni correction (p = 0.05):"), width=9),
                                        dbc.Col(dbc.Checkbox(id=app_con.ids.bt_bonferroni_cor, value=False), width=3)
                                    ])
                                ])
                            ])
                        ])
                    ])
            ]),
            # Modify plot settings
            dbc.Collapse(id=app_con.ids.coll_options, is_open=True, children=[
                dbc.Row([
                    dbc.Col([
                        dash.html.Span("Pick a color palette"),
                        dash.dcc.Dropdown(
                            id=app_con.ids.dd_color_pallette,
                            clearable=False,
                            options=[{'label': it.replace('_', ' '), 'value': it} for it in app_con.color_palettes],
                            placeholder="Default palette"
                        )
                    ]),
                    dbc.Col(children=[
                        dash.html.Span("Pick a dot size"),
                        dash.dcc.Slider(id=app_con.ids.slider_dot_size, min=1, max=20, step=1, value=10)
                    ]),
                    dbc.Col(children=[
                        dash.html.Span("Pick a opacity"),
                        dash.dcc.Slider(id=app_con.ids.slider_opacity, min=0.2, max=1.0, step=0.1, value=0.9)
                    ]),
                    dbc.Col(children=[
                        dash.html.Span("Gene color by origin order"),
                        dbc.Checkbox(id=app_con.ids.bt_origin_order, value=False)
                    ])
                ]),
                dbc.Row([
                    dash.html.Hr()
                ]),
                # Group and selection stuff.
                dbc.Row([
                    dbc.Col(width=2, children=[
                        dbc.Row([
                            dbc.Button("Create new group", id=app_con.ids.bt_di_ng, color="primary",
                                       style={'margin': '3px'}),
                        ]),
                        dbc.Row([
                            dbc.Button("Delete group", id=app_con.ids.bt_delete_group, color="danger",
                                       style={'margin': '3px'}),
                        ]),
                        dbc.Modal(id=app_con.ids.di_ng, children=[
                            dbc.ModalHeader("Create new group"),
                            dbc.ModalBody([
                                dash.html.Br(),
                                dash.html.Span("Choose a group name"),
                                dbc.Input(id=app_con.ids.di_ng_name, placeholder="Group name"),
                                dash.html.Br(),
                                dash.html.Span("Choose a group color"),
                                dbc.Input(type="color", id=app_con.ids.di_ng_color, value="#000000"),
                                dash.html.Br()
                            ]),
                            dbc.ModalFooter([
                                dbc.Button("Cancel", id=app_con.ids.di_ng_cancel, className="ml-auto"),
                                dbc.Button("Create", id=app_con.ids.di_ng_create, className="ml-auto"),
                            ])
                        ])
                    ]),
                    dbc.Col(width=3, children=[
                        dash.html.Span("Current working group"),
                        dash.dcc.Dropdown(
                            id=app_con.ids.dd_active_group,
                            clearable=True,
                            placeholder="Default group",
                            value=None
                        )
                    ]),
                    dbc.Col(width=3, children=[
                        dash.html.Span("Select modify tool"),
                        dash.dcc.Dropdown(
                            id=app_con.ids.dd_modify_tool,
                            clearable=False,
                            options=[{'label': 'Add gene', 'value': 'add'}, {'label': 'No changes', 'value': 'pause'},
                                     {'label': 'Remove gene', 'value': 'remove'}],
                            value="add"
                        )
                    ]),
                    dbc.Col(width={"size": 2}, children=[
                        dbc.Row(dash.html.Span("Select by custom list")),
                        dbc.Row(dbc.Button("Add list", id=app_con.ids.bt_add_list, color="primary",
                                           style={'margin': '3px'}))
                    ]),
                    dbc.Modal(id=app_con.ids.di_add_list, children=[
                        dbc.ModalHeader("Select by id names"),
                        dbc.ModalBody([
                            dash.html.Br(),
                            dash.html.Span("Select items by there ID. You be able to hand over a list of items like:" +
                                           "'ID1,ID2,... '"),
                            dbc.Textarea(id=app_con.ids.di_add_list_input, className="mb-3",
                                         placeholder="Type here ..."),
                            dash.html.Br(),
                        ]),
                        dbc.ModalFooter([
                            dbc.Button("Cancel", id=app_con.ids.di_add_list_cancel, className="ml-auto"),
                            dbc.Button("Create", id=app_con.ids.di_add_list_apply, className="ml-auto"),
                        ])
                    ]),
                    dbc.Col(width={"size": 2}, children=[
                        dbc.Row(dash.html.Span("Cluster direction")),
                        dbc.Row(style={'paddingTop': '0.5em'}, children=[
                            dbc.Col(dash.html.Span("Gene", className="float-end")),
                            dbc.Col(dbc.Switch(id=app_con.ids.sw_cluster_dir,
                                               value=False, style={'margin': '0px', 'padding': '0px', 'width': "100%",
                                                                   'height': '100%'},
                                               input_style={'margin': '0px', 'padding': '0px', 'width': "100%",
                                                            'height': '100%'}
                                               ), style={'maxWidth': '5em'}),
                            dbc.Col(dash.html.Span("Taxa")),
                        ])
                    ]),
                ]),
            ]),
            # PhyloProfile pre settings
            dbc.Collapse(id=app_con.ids.coll_pyhlo_opt, is_open=False, children=[
                dash.html.Br(),
                dbc.Row([
                    dbc.Col([
                        dash.html.Span("Pick subset preset"),
                        dash.dcc.Dropdown(
                            id=app_con.ids.dd_live_subset,
                            clearable=False,
                            options=[{'label': it, 'value': it} for it in app_con.const.gene_subset],
                            value=app_con.const.gene_subset[0]
                        )
                    ]),
                    dbc.Col([
                        dash.html.Span("Pick a arrangement preset"),
                        dash.dcc.Dropdown(
                            id=app_con.ids.dd_live_arrange,
                            clearable=False,
                            options=[{'label': it, 'value': it} for it in app_con.const.gene_arrange],
                            value=app_con.const.gene_arrange[0]
                        )
                    ]),
                    dbc.Col(className="position-relative", children=[
                        dash.html.Div(className="position-absolute bottom-0 start-10", children=[
                            dbc.Button("Update preset", id=app_con.ids.bt_update_pyhlo, color="success")
                        ])
                    ])
                ]),
                dash.html.Br()
            ]),
            # Group Table
            dbc.Row([
                dash.html.Hr(),
                dash.dash_table.DataTable(
                    id=app_con.ids.table_selection,
                    page_size=50,
                    style_header={'textAlign': 'left'},
                    style_table={'overflowX': 'auto', 'height': 'auto'},
                    style_cell={'textAlign': 'left'},
                    columns=[{"name": "Group name", "id": "name"},
                             {"name": "Group color", "id": "color"},
                             {"name": "Group member", "id": "geneID"}
                             ]
                )
            ])
        ]),
        # Download results section.
        dbc.Tab(tab_id=app_con.ids.tab_download, label='Downloads', children=[
            dbc.Row([
                dbc.Col([
                    dbc.Card([
                        dbc.CardHeader([
                            dash.html.H3("Download uGene groups."),
                        ]),
                        dbc.CardBody([
                            dbc.Row([
                                dash.html.P(
                                    "Download selection table content. Here you ara be able to download your " +
                                    "handcrafted groups, based on the cluster selection. Two different file formats " +
                                    "are supported."),
                            ]),
                            dash.html.Hr(),
                            dbc.Row([
                                dbc.Col(dash.html.Label("Phyloprofile compatible categories")),
                                dbc.Col([
                                    dash.html.Div(className="hstack gap-2", children=[
                                        dbc.Button("Download .cat", id=app_con.ids.bt_down_cat, color="secondary"),
                                        dbc.Spinner(dash.html.Div(id=app_con.ids.spin_down_cat, children="...."),
                                                    size="sm")
                                    ])
                                ])
                            ]),
                            dash.html.Br(),
                            dbc.Row([
                                dbc.Col(dash.html.Label("Gene groups export.")),
                                dbc.Col([
                                    dash.html.Div(className="hstack gap-2", children=[
                                        dbc.Button("Download .csv", id=app_con.ids.bt_down_groups, color="secondary"),
                                        dbc.Spinner(dash.html.Div(id=app_con.ids.spin_down_groups, children="...."),
                                                    size="sm")
                                    ])
                                ])
                            ])
                        ])
                    ])
                ]),
                dbc.Col([
                    dbc.Card([
                        dbc.CardHeader([
                            dash.html.H3("Download main Profiles")
                        ]),
                        dbc.CardBody([
                            dbc.Row([
                                dash.html.P(
                                    "Here you can get the full profile in .cluster.csv or .phyloprofile format. " +
                                    "These files contain the full molten data of the phylogenetic profile."),
                            ]),
                            dash.html.Hr(),
                            dbc.Row([
                                dbc.Col(dash.html.Label("Select gene arrangement"), width=3),
                                dbc.Col(width=9, children=[
                                    dash.dcc.Dropdown(
                                        id=app_con.ids.dd_opt_arrange,
                                        clearable=False,
                                        options=[{'label': it, 'value': it} for it in app_con.const.gene_arrange],
                                        value=app_con.const.gene_arrange[0]
                                    )
                                ])
                            ]),
                            dash.html.Br(),
                            dbc.Row([
                                dbc.Col(dash.html.Label("Select gene subset"), width=3),
                                dbc.Col(width=9, children=[
                                    dash.dcc.Dropdown(
                                        id=app_con.ids.dd_opt_subset,
                                        clearable=False,
                                        options=[{'label': it, 'value': it} for it in app_con.const.gene_subset],
                                        value=app_con.const.gene_subset[0]
                                    )
                                ])
                            ]),
                            dash.html.Hr(),
                            dbc.Row([
                                dbc.Col(dash.html.Span("Get the phyloprofile file")),
                                dbc.Col([
                                    dash.html.Div(className="hstack gap-2", children=[
                                        dbc.Button("Download .phyloprofile", id=app_con.ids.bt_down_phylo,
                                                   color="secondary"),
                                        dbc.Spinner(dash.html.Div(id=app_con.ids.spin_down_phylo, children="...."),
                                                    size="sm")
                                    ])
                                ])
                            ]),
                            dash.html.Br(),
                            dbc.Row([
                                dbc.Col(dash.html.Span("Get the full profile as .cluster.csv")),
                                dbc.Col([
                                    dash.html.Div(className="hstack gap-2", children=[
                                        dbc.Button("Download Cluster", id=app_con.ids.bt_down_cluster,
                                                   color="secondary"),
                                        dbc.Spinner(dash.html.Div(id=app_con.ids.spin_down_cluster, children="...."),
                                                    size="sm")
                                    ])
                                ])
                            ])
                        ])
                    ])
                ])
            ])
        ]),
    ]),

    # These are output dummies and user feedback dummies.
    dash.dcc.Store(id=app_con.ids.user_message, storage_type='memory'),
    dash.dcc.Store(id=app_con.ids.dummy_message, storage_type='memory'),
    dash.dcc.Store(id=app_con.ids.dummy_restyle, storage_type='memory')
])
app.title = "uGene"


# -------------------------------------- Dash backend callback definitions ---------------------------------------------
@app.callback(
    dash.dependencies.Output(app_con.ids.filename_csv, 'data'),
    dash.dependencies.Output(app_con.ids.bt_show_result, 'color'),
    dash.dependencies.Output(app_con.ids.spin_load_csv, 'children'),
    dash.dependencies.Output(app_con.ids.user_message, 'data', allow_duplicate=True),
    dash.dependencies.Input(app_con.ids.bt_load_csv, 'n_clicks'),
    dash.dependencies.State(app_con.ids.user_location, 'href'),
    prevent_initial_call=True
)
def load_csv_phyloprofile(click, url):
    """ Load a csv file or a phyloprofile file into the app_con.data container.
    :param click: Just the number of n_clicks to trigger the callback.
    :param url: String with current browser url to check localhost ues. Non localhost usage is not supported jet.
    :return: A tuple of three strings. Fist filename, second enable button statement, third short display filename.
    """
    if "127.0.0.1" not in url:
        er_ms = "Only local host use is jet supported, because dash upload can´t handel big files jet."
        print(er_ms)
        return dash.no_update, dash.no_update, dash.no_update, er_ms

    file_type = (("PhyloProfile", ".phyloprofile"), ("csv-file", ".csv"))

    filename = askOpen("Open profile", file_type)

    if filename[-len('.cluster.csv'):] == ".cluster.csv" or not filename:
        er_ms = "Error ! Input file is not valid."
        print(er_ms)
        return dash.no_update, dash.no_update, dash.no_update, er_ms

    # Handle PhyloProfile or csv input.
    if filename[-len('.phyloprofile'):] == '.phyloprofile':
        temp_df = pd.read_csv(filename, sep="\t")
        filename = filename.replace(".phyloprofile", ".csv")
        app_con.data[filename] = temp_df
    else:
        app_con.data[filename] = pd.read_csv(filename)

    return filename, 'success', "Current file " + str(filename.replace('\\', '/').split('/')[-1]), dash.no_update


@app.callback(
    dash.dependencies.Output(app_con.ids.filename_cluster, 'data', allow_duplicate=True),
    dash.dependencies.Output(app_con.ids.spin_load_cluster, 'children'),
    dash.dependencies.Output(app_con.ids.user_message, 'data', allow_duplicate=True),
    dash.dependencies.Input(app_con.ids.bt_load_cluster, 'n_clicks'),
    dash.dependencies.State(app_con.ids.user_location, 'href'),
    prevent_initial_call=True
)
def loadCluster(n_click_load, url):
    """ Load a .cluster.csv file into the app_con.data container. Load uGeneCore.py output files.
        These files have to hold the cluster columns gene1d_x, gene2d_x, gene2d_y, gene3d_x, gene3d_y, gene3d_z and as
        well all standard PhyloProfile-File columns like geneID, orthoID, ncbiID.
    :param n_click_load: Just the number of n_clicks to trigger the callback.
    :param url: String with current browser url to check localhost ues. Non localhost usage is not supported jet.
    :return: A tuple of two strings. Fist filename, second short display filename.
    """
    if "127.0.0.1" not in url:
        er_ms = "Only local host use is jet supported, because dash upload can´t handel big files jet."
        print(er_ms)
        return dash.no_update, dash.no_update, er_ms

    filename = askOpen("Open cluster file", (("Cluster file", ".cluster.csv"),))

    if not filename:
        er_ms = "Error ! Input file is not valid."
        print(er_ms)
        return dash.no_update, dash.no_update, er_ms

    app_con.data[filename] = pd.read_csv(filename)

    return filename, "Current file " + str(filename.replace('\\', '/').split('/')[-1]), dash.no_update


@app.callback(
    dash.dependencies.Output(app_con.ids.dict_selection, 'data', allow_duplicate=True),
    dash.dependencies.Output(app_con.ids.spin_load_groups, 'children'),
    dash.dependencies.Output(app_con.ids.user_message, 'data', allow_duplicate=True),
    dash.dependencies.Input(app_con.ids.bt_load_groups, 'n_clicks'),
    dash.dependencies.State(app_con.ids.user_location, 'href'),
    prevent_initial_call=True
)
def loadGroups(n_click_load, url):
    """ Load a .cluster.csv file into the app_con.data container.
        These files have to hold the cluster columns gene1d_x, gene2d_x, gene2d_y, gene3d_x, gene3d_y, gene3d_z and as
        well all standard PhyloProfile-File columns like geneID, orthoID, ncbiID.
    :param n_click_load: Just the number of n_clicks to trigger the callback.
    :param url: String with current browser url to check localhost ues. Non localhost usage is not supported jet.
    :return: A tuple of two strings. Fist filename, second short display filename.
    """
    if "127.0.0.1" not in url:
        er_ms = "Only local host use is jet supported, because dash upload can´t handel big files jet."
        print(er_ms)
        return dash.no_update, dash.no_update, er_ms

    filename = askOpen("Open groups csv", (("Cluster file", ".cat.txt"), ("Cluster file ", "*")))

    # df_sel = pd.DataFrame.from_dict(data_sel, orient='index').reset_index(names=["geneID"])
    # df_sel = df_sel[[it for it in app_con.const.phylo_cat_col if it in df_sel.columns]]

    if not filename:
        er_ms = "Error ! Input file is not valid."
        print(er_ms)
        return dash.no_update, dash.no_update, er_ms

    df_sel = pd.read_csv(filename)
    if not all([it in app_con.const.phylo_cat_col for it in df_sel.columns]):
        er_ms = "Error ! Category file do not holds expected column names"
        print(er_ms)
        return dash.no_update, dash.no_update, er_ms

    return df_sel.set_index('geneID').to_dict(orient='index'), "Current file " + str(
        filename.replace('\\', '/').split('/')[-1]), dash.no_update


@app.callback(
    dash.dependencies.Output(app_con.ids.filename_cluster, 'data', allow_duplicate=True),
    dash.dependencies.Output(app_con.ids.spin_show_result, 'children'),
    dash.dependencies.Output(app_con.ids.user_message, 'data', allow_duplicate=True),
    dash.dependencies.Input(app_con.ids.bt_show_result, 'n_clicks'),
    dash.dependencies.State(app_con.ids.filename_csv, 'data'),
    dash.dependencies.State(app_con.ids.input_int_nn, 'value'),
    dash.dependencies.State(app_con.ids.input_float_md, 'value'),
    dash.dependencies.State(app_con.ids.dd_use_metric, 'value'),
    dash.dependencies.State(app_con.ids.bt_adv_opt, 'value'),
    dash.dependencies.State(app_con.ids.text_adv_opt, 'value'),
    dash.dependencies.State(app_con.ids.sw_cluster_dir, 'value'),
    prevent_initial_call=True
)
def runCluster(n_click, filename_csv, n_neighbors, min_dist, metric, use_adv_opt, adv_opt, taxa_main):
    """ Function to submit a cluster job.
    :param n_click: Not used Int of number button hits. Needed to trigger the function.
    :param filename_csv: A string which belongs to the filename with including path and refers to an app_con.data key.
    :param n_neighbors: Int for the n_neighbors param from umap cluster algorithm.
    :param min_dist: Float for the min_dist param from umap cluster algorithm.
    :param metric: String to set up the used metric by the umap cluster algorithm.
    :param bool use_adv_opt: Flac to enable use advance options.
    :param str adv_opt: String with all given additional options. Standard settings will be overridden.
    :param taxa_main: Flac to switch to taxa ids. Effects the cluster of taxa will be displayed.
    :return: Tuple of strings. Fist, new filename, to plot the current clustered data. Second, placeholder.
    """

    if not filename_csv or filename_csv not in app_con.data or app_con.data[filename_csv] is None:
        er_ms = "Error runCluster()! No input file selected."
        print(er_ms)
        return dash.no_update, dash.no_update, er_ms

    # Create the current task with in the standard jobs. For uGeneCore.py external pyhon use, take this serious.
    gene_task = {'y_axis': 'geneID',
                 'x_axis': 'ncbiID',
                 'values': ['FAS_F', 'FAS_B'],
                 'jobs': [{'job_name': 'gene', 'n_components': 1, 'n_neighbors': n_neighbors, 'min_dist': min_dist,
                           'metric': metric},
                          {'job_name': 'gene', 'n_components': 2, 'n_neighbors': n_neighbors, 'min_dist': min_dist,
                           'metric': metric},
                          {'job_name': 'gene', 'n_components': 3, 'n_neighbors': n_neighbors, 'min_dist': min_dist,
                           'metric': metric}]}
    taxa_task = None if not taxa_main else {
        'x_axis': 'geneID',
        'y_axis': 'ncbiID',
        'values': ['FAS_F', 'FAS_B'],
        'jobs': [
            {'job_name': 'taxa', 'n_components': 1, 'n_neighbors': n_neighbors, 'min_dist': min_dist, 'metric': metric},
            {'job_name': 'taxa', 'n_components': 2, 'n_neighbors': n_neighbors, 'min_dist': min_dist, 'metric': metric},
            {'job_name': 'taxa', 'n_components': 3, 'n_neighbors': n_neighbors, 'min_dist': min_dist, 'metric': metric}
        ]}

    # Process advanced options
    if use_adv_opt:
        gene_task, taxa_task = processAdvancedOptions(gene_task, taxa_task, adv_opt, app_con.const.task_arguments)
    # Do the cluster job with the external main file.
    df = uGene.mainAnalytics(app_con.data[filename_csv], **gene_task)

    if taxa_task:
        df = uGene.mainAnalytics(df, **taxa_task)

    # Check if all required cluster data have being produced.
    if not all([it in df.columns for it in app_con.const.cluster_col_names]):
        er_ms = "Error runCluster() ! Not all required cluster data being produced.\n" + \
                "Pleas use uGeneCore.py manual and check for column names " + ",".join(app_con.const.cluster_col_names)
        print(er_ms)
        return dash.no_update, dash.no_update, er_ms

    # Change file name extension, because now it is a full functional .cluster.csv file.
    filename = filename_csv.replace(".csv", ".cluster.csv")
    app_con.data[filename] = df
    return filename, "", dash.no_update


@app.callback(
    dash.dependencies.Output(app_con.ids.plot_2d, 'figure', allow_duplicate=True),
    dash.dependencies.Output(app_con.ids.plot_2d, 'config'),
    dash.dependencies.Output(app_con.ids.plot_3d, 'figure', allow_duplicate=True),
    dash.dependencies.Output(app_con.ids.plot_3d, 'config'),
    dash.dependencies.Output(app_con.ids.scatter_3d, 'figure', allow_duplicate=True),
    dash.dependencies.Output(app_con.ids.scatter_3d, 'config'),
    dash.dependencies.Output(app_con.ids.tabs_ugene_viewer, 'active_tab', allow_duplicate=True),
    dash.dependencies.Output(app_con.ids.sw_cluster_dir, 'value'),
    dash.dependencies.Output(app_con.ids.store_order_id, 'data'),
    dash.dependencies.Input(app_con.ids.filename_cluster, 'data'),
    dash.dependencies.Input(app_con.ids.sw_cluster_dir, 'value'),
    dash.dependencies.Input(app_con.ids.bt_origin_order, 'value'),
    dash.dependencies.State(app_con.ids.slider_dot_size, 'value'),
    dash.dependencies.State(app_con.ids.slider_opacity, 'value'),
    dash.dependencies.State(app_con.ids.dd_color_pallette, 'value'),
    dash.dependencies.State(app_con.ids.input_int_nn, 'value'),
    dash.dependencies.State(app_con.ids.input_float_md, 'value'),
    dash.dependencies.State(app_con.ids.dd_use_metric, 'value'),
    dash.dependencies.State(app_con.ids.bt_adv_opt, 'value'),
    dash.dependencies.State(app_con.ids.text_adv_opt, 'value'),
    dash.dependencies.Input(app_con.ids.dict_selection, 'data'),
    prevent_initial_call=True
)
def createPrePlot(p_data, taxa_main, origin_order, dot_size, dot_opacity, dot_palette, n_neighbors, min_dist, metric,
                  use_adv_opt, adv_opt, data_sel):
    """ Create all three main plots for 2d , 3d and scatter matrix cluster view.
    :param p_data: A string that points to the data in app_con.data. Basically the filename with including path.
    :param taxa_main: Flac to switch to taxa ids. Effects the cluster of taxa will be displayed.
    :param origin_order: Flac to enable origin order. Option how to color and order ids.
    :param dot_size: Int of the current dot size of scatter points.
    :param dot_opacity: Float of dot opacity.
    :param dot_palette: A string with the name of a color palette. Have to be one key of app_con.color_palettes.
    :param n_neighbors: Int for the n_neighbors param from umap cluster algorithm.
    :param min_dist: Float for the min_dist param from umap cluster algorithm.
    :param metric: String to set up the used metric by the umap cluster algorithm.
    :param bool use_adv_opt: Flac to enable use advance options.
    :param str adv_opt: String with all given additional options. Standard settings will be overridden.
    :param data_sel: A dictionary with data of selected genes as well there group name. Gene names are main keys.
    :return: A tuple with four items. First all created plots and at least a current tap setting to change the view.
    """
    if p_data not in app_con.data or "gene2d_x" not in app_con.data[p_data].columns or \
            "gene2d_y" not in app_con.data[p_data].columns or "gene3d_x" not in app_con.data[p_data].columns or \
            "gene3d_y" not in app_con.data[p_data].columns or "gene3d_z" not in app_con.data[p_data].columns:
        raise dash.exceptions.PreventUpdate()

    # Taxa cluster data have not to been in to the dataset. They will be added in this case.
    if taxa_main:
        if not all([it in app_con.data[p_data].columns for it in
                    ['taxa2d_x', 'taxa2d_y', 'taxa3d_x', 'taxa3d_y', 'taxa3d_z']]):
            print(" Start taxa cluster !!!!!", app_con.data[p_data].columns)
            taxa_task = {
                'x_axis': 'geneID',
                'y_axis': 'ncbiID',
                'values': ['FAS_F', 'FAS_B'],
                'jobs': [
                    {'job_name': 'taxa', 'n_components': 1, 'n_neighbors': n_neighbors, 'min_dist': min_dist,
                     'metric': metric},
                    {'job_name': 'taxa', 'n_components': 2, 'n_neighbors': n_neighbors, 'min_dist': min_dist,
                     'metric': metric},
                    {'job_name': 'taxa', 'n_components': 3, 'n_neighbors': n_neighbors, 'min_dist': min_dist,
                     'metric': metric}
                ]}

            if use_adv_opt:
                dummy, taxa_task = processAdvancedOptions({'jobs': []}, taxa_task, adv_opt,
                                                          app_con.const.task_arguments)

            # Run taxa cluster job.
            app_con.data[p_data] = uGene.mainAnalytics(app_con.data[p_data], **taxa_task)

    # Manage column names. All column names are strings.
    if taxa_main:
        col_id = app_con.const.taxa_col
        col_x_1d = "taxa1d_x"
        col_x_2d = "taxa2d_x"
        col_y_2d = "taxa2d_y"
        col_x_3d = "taxa3d_x"
        col_y_3d = "taxa3d_y"
        col_z_3d = "taxa3d_z"
    else:
        col_id = app_con.const.gene_col
        col_x_1d = "gene1d_x"
        col_x_2d = "gene2d_x"
        col_y_2d = "gene2d_y"
        col_x_3d = "gene3d_x"
        col_y_3d = "gene3d_y"
        col_z_3d = "gene3d_z"

    print("Info createPrePlot(): Start plotting.")

    # Get the profile data form app_con.data source.
    df = app_con.data[p_data].drop_duplicates(subset=['geneID'])

    if not origin_order:
        df = df.sort_values([col_x_1d]).reset_index(drop=True)
    # Create custom colors for each gene.
    colors = colorRampPalette(app_con.color_palettes[dot_palette if dot_palette else 'Rainbow'], len(df[col_id]))
    order_main_ids = df[col_id].to_list()

    border = None
    if data_sel:
        border = [data_sel[it]['color'] if it in data_sel else app_con.const.color_limpid for it in list(df[col_id])]

    data_2d, data_3d, data_matrix_3d = scatterDataExpress(df=df, x_2d=col_x_2d, y_2d=col_y_2d, x_3d=col_x_3d,
                                                          y_3d=col_y_3d, z_3d=col_z_3d, color=colors, border=border,
                                                          dot_size=dot_size, opacity=dot_opacity, customdata=[col_id])
    layout_2d = go.Layout(
        title='2D Scatter Plot',
        xaxis=dict(title='X-Axis'),
        yaxis=dict(title='Y-Axis'),
        template='plotly_white'
    )
    layout_3d = go.Layout(
        legend={'title': {'text': col_id}, 'tracegroupgap': 0},
        scene={'domain': {'x': [0.0, 1.0], 'y': [0.0, 1.0]},
               'xaxis': {'title': {'text': col_x_3d}},
               'yaxis': {'title': {'text': col_y_2d}},
               'zaxis': {'title': {'text': col_z_3d}}},
        template='plotly_white'
    )
    layout_matrix_3d = go.Layout(
        legend={'title': {'text': 'geneID'}, 'tracegroupgap': 0},
    )

    fig_2d = go.Figure(data=data_2d, layout=layout_2d)
    fig_3d = go.Figure(data=data_3d, layout=layout_3d)
    fig_matrix_3d = go.Figure(data=data_matrix_3d, layout=layout_matrix_3d)

    print("Info createPrePlot(): Done plotting.")

    plot_config = {
        'toImageButtonOptions': {
            'format': 'jpeg',  # jpeg, png, svg, webp
            'filename': 'uGeneScreenshot',
            'height': 1080,
            'width': 1920,
            'scale': 10
        }
    }
    return fig_2d, plot_config, fig_3d, plot_config, fig_matrix_3d, plot_config, app_con.ids.tab_viewer, taxa_main, \
           order_main_ids


@app.callback(
    dash.dependencies.Output(app_con.ids.active_color_palette, 'data'),
    dash.dependencies.Input(app_con.ids.dd_color_pallette, 'value'),
    dash.dependencies.State(app_con.ids.filename_cluster, 'data'),
    dash.dependencies.State(app_con.ids.store_order_id, 'data'),
    prevent_initial_call=True
)
def updateColorPalette(palette, filename, order_main_ids):
    """ Funktion to set a new custom color vector to based on the selected color palette.
    :param str palette: A string with the name of a color palette. Have to be one key of app_con.color_palettes.
    :param str filename: A string which belongs to the filename with including path and refers to an app_con.data key.
    :param list order_main_ids: Ordered list of ids. These ids are the in use ids form the main dataframe.
    :return: A list with containing strings, holding hex color codes.
    """
    if not palette or not filename or not order_main_ids:
        print("Error updateColorPalette(). No valid arguments.")
        raise dash.exceptions.PreventUpdate()

    # Create an individual color for each gene.
    colors = colorRampPalette(app_con.color_palettes[palette], len(order_main_ids))
    return colors


@app.callback(
    dash.dependencies.Output(app_con.ids.table_selection, 'data'),
    dash.dependencies.Output(app_con.ids.table_selection, 'style_data_conditional'),
    dash.dependencies.Input(app_con.ids.dict_selection, 'data'),
    dash.dependencies.State(app_con.ids.filename_cluster, 'data'),
    prevent_initial_call=True
)
def updateTable(data_sel, filename):
    """ Update the current selection table content.
    :param data_sel: A dictionary with data of selected genes as well there group name. Gene names are main keys.
    :param filename: A string which belongs to the filename with including path and refers to an app_con.data key.
    :return: A tuple with first table content data list and second a list of style constrains to color the table.
    """
    if not filename:
        raise dash.exceptions.PreventUpdate()
    if not data_sel:
        return [], []

    df_sel = pd.DataFrame.from_dict(data_sel, orient='index').reset_index(names=["geneID"])

    data = df_sel.groupby('name').apply(lambda x: {col: ",".join(set(x[col])) for col in x.columns})

    # Color the background of the group color field.
    style_data = [{
        'if': {
            'filter_query': '{name} =' + str(it['name']),
            'column_id': 'color'
        },
        'backgroundColor': str(it['color'])
    } for it in data]

    return [it for it in data], style_data


@app.callback(
    dash.dependencies.Output(app_con.ids.phylo_iframe, 'srcDoc'),
    dash.dependencies.Output(app_con.ids.phylo_iframe, 'src'),
    dash.dependencies.Output(app_con.ids.user_message, 'data', allow_duplicate=True),
    dash.dependencies.Input(app_con.ids.bt_update_pyhlo, 'n_clicks'),
    dash.dependencies.Input(app_con.ids.tabs_show_data, 'active_tab'),
    dash.dependencies.State(app_con.ids.dict_selection, 'data'),
    dash.dependencies.State(app_con.ids.filename_cluster, 'data'),
    dash.dependencies.State(app_con.ids.dd_live_arrange, 'value'),
    dash.dependencies.State(app_con.ids.dd_live_subset, 'value'),
    prevent_initial_call=True
)
def updatePhylo(clicks, active_tab, data_sel, filename, opt_arrange, opt_subset):
    """ This function starts and stop phyloprofile to make it realtime accessible to the dashboard.
    :param clicks: Not used Int of number button hits. Needed to trigger the function manual.
    :param active_tab: String with the id of the active tab.
    :param data_sel: A dictionary with data of selected genes as well there group name. Gene names are main keys.
    :param filename: A string which belongs to the filename with including path and refers to an app_con.data key.
    :param opt_arrange: String of the current selected option how to order genes. One of app_con.const.gene_arrange
    :param opt_subset: String of the current selected gene subset setting. One of app_con.const.gene_subset
    :return: Tuple with two string, that controls the iframe src and srcDoc.
    """
    if not filename or filename not in app_con.data:
        raise dash.exceptions.PreventUpdate()

    filename_core = filename.replace('\\', '/').split('/')[-1].replace('.cluster.csv', '').replace('.csv', '')
    if not filename_core:
        raise dash.exceptions.PreventUpdate()

    path_dir = os.path.dirname(os.path.abspath(__file__)).replace('\\', '/')
    # Check location and temp_folder
    if not (os.path.exists(path_dir + str(app_con.const.temp_folder))):
        er_ms = "Error updatePhylo()! Folder " + str(app_con.const.temp_folder) + " not found"
        print(er_ms)
        return dash.no_update, dash.no_update, er_ms

    # The port must be available and different from the previously used port. If not, the iframe won't change.
    port = freePort()
    # Check for running server.
    if filename + ".phyloprocess" in app_con.data:
        if app_con.data[filename + ".phyloprocess"] is not None:
            app_con.data[filename + ".phyloprocess"].kill()
            del app_con.data[filename + ".phyloprocess"]

    if active_tab != app_con.ids.tab_show_phyloprofile:
        print("Info: Phyloprofile off.")
        # Return First -> iframe:srcDoc, Second -> iframe:src, Third -> collapse:is_open
        return "<p> Load Phyloprofile ...</p>", 'about:blank', dash.no_update

    # Create all needed file locations
    path_t_data = path_dir + app_con.const.temp_folder + filename_core + '.phyloprofile'
    path_t_cat = path_dir + app_con.const.temp_folder + filename_core + '.phylocategories'
    path_t_config = path_dir + app_con.const.temp_folder + filename_core + '.phyloconfig.yml'

    if data_sel:
        df_sel = pd.DataFrame.from_dict(data_sel, orient='index').reset_index(names=["geneID"])
        df_sel = df_sel[[it for it in app_con.const.phylo_cat_col if it in df_sel.columns]]

        # Fix PhyloProfile uses categories like filters.
        ser_id = app_con.data[filename]['geneID'].drop_duplicates()

        # Fix only gene ids are allowed. Maybe data_sel could contain ncbiID´s.
        df_sel = df_sel[df_sel['geneID'].isin(ser_id)]

        # Continue with PhyloProfile fix
        ser_id = ser_id[~ser_id.isin(data_sel.keys())]

        df_fix = pd.DataFrame(columns=df_sel.columns)
        df_fix[app_con.const.phylo_cat_col[0]] = ser_id.reset_index(drop=True)
        df_fix[app_con.const.phylo_cat_col[1]] = len(ser_id) * [app_con.const.phylo_defaild_cat]
        df_fix[app_con.const.phylo_cat_col[2]] = len(ser_id) * [app_con.const.color_phylo_defauild]

        # Save new category dataframe
        pd.concat([df_sel, df_fix], ignore_index=True).to_csv(path_t_cat, header=False, index=False, sep='\t')

    # Produce .phyloprofile data.
    df = processOptionDf(app_con.data[filename], data_sel, opt_arrange, opt_subset)
    df = df[[it for it in app_con.const.phylo_col if it in app_con.data[filename].columns]]
    df.to_csv(path_t_data, index=False, sep='\t')

    # Create phylo config file
    with open(path_t_config, 'w', encoding='utf-8') as config:
        config.write('---\n')
        config.write('mainInput: ' + str(path_t_data) + '\n')
        config.write('domainInput: NULL\n')
        config.write('fastaInput: NULL\n')
        config.write('treeInput: NULL\n')
        config.write('ordering: 0\n')
        if data_sel:
            config.write('geneCategory: ' + str(path_t_cat) + '\n')
            config.write('colorByGroup: 1\n')
        # New after update.
        config.write('launchBrowser: False\n')
        config.write('host: 127.0.0.1\n')
        config.write('port: ' + str(port) + '\n')
        config.write('...\n')

    # Execute phyloprofile
    phylo_process = subprocess.Popen(
        ['R', '-e', 'library(PhyloProfile); runPhyloProfile(configFile ="' + path_t_config + '");'],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE)

    ip = None
    timeout = 120
    t_start = time.time()

    # We need to wait until the phyloprofile is started. No ip mean no successful running phyloprofile.
    while not ip:
        stderr = phylo_process.stderr.readline().decode('utf-8')
        status_code = phylo_process.poll()

        print(stderr)
        if stderr.find("Listening on") != -1:
            ip = stderr.replace("Listening on ", "")
            break
        elif time.time() - t_start > timeout:
            # We set a timeout to prevent this callback to run forever.
            print("Fail to start phyloprofile! Make sure phyloprofile is installed and callable by commandline.")
            phylo_process.kill()
            phylo_process = None
            ip = None
            break
        elif status_code:
            # A status_code means the function is finished for some reasons and tells back his finish code.
            print("Fail to run PhyloProfile !", status_code)
            phylo_process = None
            ip = None
            break

    # Check if PhyloProfile runs successful.
    if not phylo_process or not ip:
        return "<p> Fail to use Phyloprofile!</p><p>Check Phyloprofile installation, as well server logs.</p>", \
               'about:blank', dash.no_update
    else:
        app_con.data[filename + ".phyloprocess"] = phylo_process
        print("Info: Phyloprofile is running.")
        # Return First -> iframe:srcDoc, Second -> iframe:src, Third -> collapse:is_open
        return None, ip, dash.no_update


@app.callback(
    dash.dependencies.Output(app_con.ids.coll_options, 'is_open'),
    dash.dependencies.Output(app_con.ids.coll_pyhlo_opt, 'is_open'),
    dash.dependencies.Input(app_con.ids.tabs_show_data, 'active_tab'),
    prevent_initial_call=True)
def collapseOptions(active_tab):
    """ Callback to hide and show the current relevant options below the main plots.
    :param active_tab: String with the id of the current open tab.
    :return: Tuple of two bools, which controls the display options.
    """
    if active_tab == app_con.ids.tab_show_phyloprofile:
        return False, True
    else:
        return True, False


@app.callback(
    dash.dependencies.Output(app_con.ids.coll_main_header, 'is_open'),
    dash.dependencies.Input(app_con.ids.tabs_ugene_viewer, 'active_tab'),
    prevent_initial_call=True)
def collapseHeader(active_tab):
    """ Callback to hide and show the page header. Into the view tab, there is no benefit to have such header.
    :param active_tab: String with the id of the current open tab.
    :return: Bool, which control the display option.
    """
    if active_tab == app_con.ids.tab_viewer:
        return False
    else:
        return True


@app.callback(
    dash.dependencies.Output(app_con.ids.bt_down_cat, 'color'),
    dash.dependencies.Output(app_con.ids.bt_down_groups, 'color'),
    dash.dependencies.Output(app_con.ids.bt_down_phylo, 'color'),
    dash.dependencies.Output(app_con.ids.bt_down_cluster, 'color'),
    dash.dependencies.Input(app_con.ids.filename_cluster, 'data'),
    dash.dependencies.Input(app_con.ids.dict_selection, 'data'),
    prevent_initial_call=True)
def enableDownButton(filename, selection):
    """ Set the color of download buttons to show if a button is usable or not.
    :param filename: A string which belongs to the filename with including path and refers to an app_con.data key.
    :param selection: A dictionary with data of selected genes as well there group name. Gene names are main keys.
    :return: Tuple within bootstrap dash button color codes.
    """
    if filename and selection and filename in app_con.data:
        return "primary", "primary", "primary", "primary"
    elif filename and filename in app_con.data:
        return "secondary", "secondary", "primary", "primary"
    else:
        return "secondary", "secondary", "secondary", "secondary"


@app.callback(
    dash.dependencies.Output(app_con.ids.spin_down_cat, 'children'),
    dash.dependencies.Input(app_con.ids.bt_down_cat, 'n_clicks'),
    dash.dependencies.State(app_con.ids.filename_cluster, 'data'),
    dash.dependencies.State(app_con.ids.dict_selection, 'data'),
    prevent_initial_call=True)
def downPhyloProfileCat(clicks, filename, data_sel):
    """ Function to download groups into PhyloProfile category format.
    :param clicks: Not used Int of number button hits. Needed to trigger the function manual.
    :param filename: A string which belongs to the filename with including path and refers to an app_con.data key.
    :param data_sel: A dictionary with data of selected genes as well there group name. Gene names are main keys.
    :return: String with download feedback.
    """
    if not filename or not data_sel or filename not in app_con.data:
        return "(no data)"

    # Prepare export data
    # df_sel = pd.DataFrame.from_dict(data_sel, orient='index').reset_index(names=["geneID"])
    # df_sel = df_sel[[it for it in app_con.const.phylo_cat_col if it in df_sel.columns]]

    df_sel = pd.DataFrame.from_dict(data_sel, orient='index').reset_index(names=["geneID"])
    df_sel = df_sel[[it for it in app_con.const.phylo_cat_col if it in df_sel.columns]]

    # Fix PhyloProfile uses categories like filters.
    ser_id = app_con.data[filename]['geneID'].drop_duplicates()

    # Fix only gene ids are allowed. Maybe data_sel could contain ncbiID´s
    df_sel = df_sel[df_sel['geneID'].isin(ser_id)]

    # Continue with PhyloProfile fix
    ser_id = ser_id[~ser_id.isin(data_sel.keys())]

    df_fix = pd.DataFrame(columns=df_sel.columns)
    df_fix[app_con.const.phylo_cat_col[0]] = ser_id.reset_index(drop=True)
    df_fix[app_con.const.phylo_cat_col[1]] = len(ser_id) * [app_con.const.phylo_defaild_cat]
    df_fix[app_con.const.phylo_cat_col[2]] = len(ser_id) * [app_con.const.color_phylo_defauild]

    export_path = askSave('Save Phyloprofile categories ', (("Cat file", ".cat.txt"), ("File", "*")))

    # Save new category dataframe
    # df_sel.to_csv(export_path, header=False, index=False, sep='\t')
    pd.concat([df_sel, df_fix], ignore_index=True).to_csv(export_path, header=False, index=False, sep='\t')

    return "(done)"


@app.callback(
    dash.dependencies.Output(app_con.ids.spin_down_groups, 'children'),
    dash.dependencies.Input(app_con.ids.bt_down_groups, 'n_clicks'),
    dash.dependencies.State(app_con.ids.filename_cluster, 'data'),
    dash.dependencies.State(app_con.ids.dict_selection, 'data'),
    prevent_initial_call=True)
def downSelectedGroups(clicks, filename, data_sel):
    """ Function to download groups into a standard csv file format, to make it accessible to and other software.
    :param clicks: Not used Int of number button hits. Needed to trigger the function manual.
    :param filename: A string which belongs to the filename with including path and refers to an app_con.data key.
    :param data_sel: A dictionary with data of selected genes as well there group name. Gene names are main keys.
    :return: String with download feedback.
    """
    if not filename or not data_sel or filename not in app_con.data:
        return "(no data)"

    # Prepare export data
    df_sel = pd.DataFrame.from_dict(data_sel, orient='index').reset_index(names=["geneID"])
    df_sel = df_sel[[it for it in app_con.const.phylo_cat_col if it in df_sel.columns]]

    export_path = askSave('Save Gene groups ', (("File", "*"),))

    df_sel.to_csv(export_path, index=False)

    return "(done)"


@app.callback(
    dash.dependencies.Output(app_con.ids.spin_down_phylo, 'children'),
    dash.dependencies.Input(app_con.ids.bt_down_phylo, 'n_clicks'),
    dash.dependencies.State(app_con.ids.filename_cluster, 'data'),
    dash.dependencies.State(app_con.ids.dict_selection, 'data'),
    dash.dependencies.State(app_con.ids.dd_opt_arrange, 'value'),
    dash.dependencies.State(app_con.ids.dd_opt_subset, 'value'),
    prevent_initial_call=True)
def downPhyloProfile(clicks, filename, data_sel, opt_arrange, opt_subset):
    """ Function to download the resulting phyloprofile according to the download settings.
    :param clicks: Not used Int of number button hits. Needed to trigger the function manual.
    :param filename: A string which belongs to the filename with including path and refers to an app_con.data key.
    :param data_sel: A dictionary with data of selected genes as well there group name. Gene names are main keys.
    :param opt_arrange: String of the current selected option how to order genes. One of app_con.const.gene_arrange
    :param opt_subset: String of the current selected gene subset setting. One of app_con.const.gene_subset
    :return: String with download feedback.
    """
    if not filename or filename not in app_con.data:
        return "(no data)"

    df = app_con.data[filename]
    # Check because processOptionDf() is not secure for this case.
    if opt_subset == "group-based" and not data_sel:
        print("Waring downClusterProfile() fail. No groups are selected.")
        return "(fail)"

    # Check because processOptionDf() is not secure for this case.
    if opt_arrange == "1d-order" and "gene1d_x" not in df.columns:
        print("Error downClusterProfile(). Column gene1d_x not in cluster dataframe.")
        return "(fail)"

    df = processOptionDf(df, data_sel, opt_arrange, opt_subset)

    # Filter by phyloprofile compatible column names.
    df = df[[it for it in app_con.const.phylo_col if it in app_con.data[filename].columns]]

    export_path = askSave('Save Phyloprofile', (("Phyloprofile", ".phyloprofile"), ("file", "*")))

    if not export_path:
        return "(abort)"

    if export_path[-len(".phyloprofile"):] != ".phyloprofile":
        export_path = export_path + ".phyloprofile"

    df.to_csv(export_path, index=False, sep='\t')

    return "(done)"


@app.callback(
    dash.dependencies.Output(app_con.ids.spin_down_cluster, 'children'),
    dash.dependencies.Input(app_con.ids.bt_down_cluster, 'n_clicks'),
    dash.dependencies.State(app_con.ids.filename_cluster, 'data'),
    dash.dependencies.State(app_con.ids.dict_selection, 'data'),
    dash.dependencies.State(app_con.ids.dd_opt_arrange, 'value'),
    dash.dependencies.State(app_con.ids.dd_opt_subset, 'value'),
    prevent_initial_call=True)
def downClusterProfile(clicks, filename, data_sel, opt_arrange, opt_subset):
    """ Function to download the resulting cluster csv file according to the download settings.
    :param clicks: Not used Int of number button hits. Needed to trigger the function manual.
    :param filename: A string which belongs to the filename with including path and refers to an app_con.data key.
    :param data_sel: A dictionary with data of selected genes as well there group name. Gene names are main keys.
    :param opt_arrange: String of the current selected option how to order genes. One of app_con.const.gene_arrange
    :param opt_subset: String of the current selected gene subset setting. One of app_con.const.gene_subset
    :return: String with download feedback.
    """
    if not filename or filename not in app_con.data:
        return "(no data)"

    df = app_con.data[filename]

    # Check because processOptionDf() is not secure for this case.
    if opt_subset == "group-based" and not data_sel:
        print("Waring downClusterProfile() fail. No groups are selected.")
        return "(fail)"

    # Check because processOptionDf() is not secure for this case.
    if opt_arrange == "1d-order" and "gene1d_x" not in df.columns:
        print("Error downClusterProfile(). Column gene1d_x not in cluster dataframe.")
        return "(fail)"

    df = processOptionDf(df, data_sel, opt_arrange, opt_subset)

    export_path = askSave('Save cluster profile', (("Cluster Profile", ".cluster.csv"), ("file", "*")))
    if not export_path:
        return "(abort)"

    # Allow .csv or .cluster.csv. Force to .cluster.csv if no .csv file is detected.
    if export_path[-len(".cluster.csv"):] != ".cluster.csv":
        if export_path[-len(".csv"):] != ".csv":
            export_path = export_path + ".cluster.csv"

    df.to_csv(export_path, index=False)
    return "(done)"


@app.callback(
    dash.dependencies.Output(app_con.ids.filename_stat, 'data', allow_duplicate=True),
    dash.dependencies.Output(app_con.ids.spin_load_stat, 'children'),
    dash.dependencies.Output(app_con.ids.user_message, 'data', allow_duplicate=True),
    dash.dependencies.Input(app_con.ids.bt_load_stat, 'n_clicks'),
    dash.dependencies.State(app_con.ids.user_location, 'href'),
    prevent_initial_call=True
)
def loadStat(n_click_load, url):
    """ Load additional statistic file into the app_con.data container.
        These files have to hold the cluster columns gene1d_x, gene2d_x, gene2d_y, gene3d_x, gene3d_y, gene3d_z and as
        well all standard phyloprofile columns like geneID, orthoID, ncbiID.
    :param n_click_load: Just the number of n_clicks to trigger the callback.
    :param url: String with current browser url to check localhost ues. Non localhost usage is not supported jet.
    :return: A tuple of two strings. First filename, second short display filename.
    """
    if "127.0.0.1" not in url:
        er_ms = "Only local host use is jet supported, because dash upload can´t handel big files jet."
        print(er_ms)
        return dash.no_update, dash.no_update, er_ms

    filename = askOpen('Add additional information', (("File", ".csv"), ("File *", "*")))

    if not filename:
        er_ms = "Error ! Input file is not valid."
        print(er_ms)
        return dash.no_update, dash.no_update, er_ms

    # Read the dataframe with the additional statistic.
    df = pd.read_csv(filename)

    if len(df) < 1 or len(df.columns) != 3:
        er_ms = "Error ! This additional data do not follow the required constrains. Review uGene doku."
        print(er_ms)
        return dash.no_update, dash.no_update, er_ms

    app_con.data[filename] = df

    return filename, "Current file " + str(filename.replace('\\', '/').split('/')[-1]), dash.no_update


@app.callback(
    dash.dependencies.Output(app_con.ids.bar_plot_aaa, 'figure'),
    dash.dependencies.Output(app_con.ids.bar_plot_bbb, 'figure'),
    dash.dependencies.Output(app_con.ids.bar_plot_ccc, 'figure'),
    dash.dependencies.Output(app_con.ids.bar_plot_ddd, 'figure'),
    dash.dependencies.Output(app_con.ids.bar_plot_eee, 'figure'),
    dash.dependencies.Output(app_con.ids.bar_plot_fff, 'figure'),
    dash.dependencies.Output(app_con.ids.bar_plot_ggg, 'figure'),
    dash.dependencies.Output(app_con.ids.bar_plot_hhh, 'figure'),
    dash.dependencies.Output(app_con.ids.bar_plot_iii, 'figure'),
    dash.dependencies.Output(app_con.ids.bar_plot_jjj, 'figure'),
    dash.dependencies.Output(app_con.ids.coll_bar_ccc, 'is_open'),
    dash.dependencies.Output(app_con.ids.coll_bar_ddd, 'is_open'),
    dash.dependencies.Output(app_con.ids.coll_bar_eee, 'is_open'),
    dash.dependencies.Output(app_con.ids.coll_bar_fff, 'is_open'),
    dash.dependencies.Output(app_con.ids.coll_bar_ggg, 'is_open'),
    dash.dependencies.Output(app_con.ids.coll_bar_hhh, 'is_open'),
    dash.dependencies.Output(app_con.ids.coll_bar_iii, 'is_open'),
    dash.dependencies.Output(app_con.ids.coll_bar_jjj, 'is_open'),
    dash.dependencies.Input(app_con.ids.store_sel_stat, 'data'),
    dash.dependencies.Input(app_con.ids.store_order_id, 'data'),
    dash.dependencies.Input(app_con.ids.active_color_palette, 'data'),
    dash.dependencies.Input(app_con.ids.bt_hypergeometirc_stat, 'value'),
    dash.dependencies.Input(app_con.ids.bt_bonferroni_cor, 'value'),
    dash.dependencies.Input(app_con.ids.input_int_bar_n_best, 'value'),
    dash.dependencies.State(app_con.ids.filename_cluster, 'data'),
    dash.dependencies.State(app_con.ids.filename_stat, 'data'),
    dash.dependencies.State(app_con.ids.dd_color_pallette, 'value'),
    prevent_initial_call=True)
def updateAdditionalData(select_ids, order_main_ids, active_colors, use_hypergeo_test, use_bonf_cor, bar_n_show,
                         filename_cluster, filename_stat, dot_palette):
    """ Update the additional statistics based on the current data selection.
    :param select_ids: List of selected items.
    :param order_main_ids: Ordered list of ids. These ids are the in use ids form the main dataframe.
    :param active_colors: A string with the name of the current active color palette. One of app_con.color_palettes.
    :param use_hypergeo_test: Flac if the hyper-geometric test is in use.
    :param use_bonf_cor: Flac to display the significant border by Bonferroni
    :param bar_n_show: Int between 1 until n. Effects how many best hists into the statistic file will be displayed.
    :param filename_cluster: A string which belongs to the filename with including path and refers to app_con.data.
    :param filename_stat: A string which belongs to the filename of the additional statistic file.
    :param dot_palette: A string with the name of a color palette. Have to be one key of app_con.color_palettes.
    :return: Tuple with 1 up to 10 plotly figures and 0 up to 8 bools with the plot visibly status.
    """

    if not filename_stat or not filename_cluster or filename_stat not in app_con.data or filename_cluster \
            not in app_con.data or not select_ids or not order_main_ids:
        raise dash.exceptions.PreventUpdate()

    df_main = app_con.data[filename_cluster]
    df_stat = app_con.data[filename_stat]

    if df_stat.columns[0] not in df_main.columns:
        print("Error updateAdditionalData()! There are no matching columns.")
        raise dash.exceptions.PreventUpdate()

    if not active_colors or len(active_colors) != len(order_main_ids):
        # Prepare colors. Create custom colors for each gene.
        active_colors = colorRampPalette(app_con.color_palettes[dot_palette if dot_palette else 'Rainbow'],
                                         len(order_main_ids))

    dict_colors = {it[0]: it[1] for it in zip(order_main_ids, active_colors)}

    # List of bar plots. Max 10 plots are usable.
    li_fig = []

    if use_hypergeo_test:
        # Filter to the current generated plot
        stat_grouped = df_stat.groupby(by=df_stat.columns[1])
        for plot_name, df_sub in stat_grouped:

            # Size of the full set for the hyper geometric test. Actually M.
            total_sub_m = len(df_sub.iloc[:, 2])

            # Filter subset of selected genes. df_sub.columns[0] == "geneID"
            df_set = df_sub[df_sub[df_sub.columns[0]].isin(select_ids)].reset_index(drop=True)

            # Size of the subset for hyper geometric test. Actually N.
            total_set_n = len(df_set.iloc[:, 2])

            # Counting the IDs and saving the statistics as a dataframe.
            set_count = df_set.groupby(df_sub.columns[2]).count()
            sub_count = df_sub[df_sub.iloc[:, 2].isin(set_count.index)].reset_index(drop=True).groupby(
                df_sub.columns[2]).count()

            # Perform the hyper-geometric test and save results into first column. Column name geneID is an artifact.
            set_count[['pValue', 'pValueAdjust']] = [
                divideHelp(1 - sst.hypergeom.sf(it[0] - 1, total_sub_m, it[1], total_set_n), it[0]) for it in
                zip(list(set_count.iloc[:, 1]), list(sub_count.iloc[:, 1]))]

            # Get number of independent test for Bonferroni correction.
            k_independent_tests = len(set_count)

            # Show at least one result
            n = bar_n_show if bar_n_show > 0 else 1
            # Filter to the most represented significant items.
            set_count = set_count.sort_values(by=set_count.columns[2], ascending=False).iloc[:n]
            filter_n_x = set_count.index.tolist()

            # Apply filter to the sub dataframe
            df_set = df_set[df_set.iloc[:, 2].isin(filter_n_x)].reset_index(drop=True)
            df_set = df_set.merge(set_count.iloc[:, 2:4], left_on=df_set.columns[2], right_index=True, how="left")

            # Create the figure.
            if not df_set.empty:
                fig = px.bar(df_set, y=df_set.columns[4],
                             x=df_set.columns[2],
                             color=df_set.columns[0],
                             title=str(plot_name),
                             custom_data=[df_set.columns[0], 'pValue'],
                             category_orders={df_set.columns[2]: filter_n_x})
            else:
                raise dash.exceptions.PreventUpdate()

            # Update to plot color scale. Using not common use of 'legendgroup'
            for itt in fig.data:

                # Update Color
                if itt['legendgroup'] in dict_colors:
                    itt['marker']['color'] = dict_colors[itt['legendgroup']]

                # Update hovertemplate
                itt['hovertemplate'] = "PointID=%{customdata[0]}<br>AnnotatedID=%{x}<br>pValue=%{customdata[1]}" + \
                                       "<extra></extra>"

            # Add a p= 0.05 confidence interval corrected by Bonferroni rule.
            if use_bonf_cor:
                interval_y = 1 - 0.05 / k_independent_tests
                fig.add_shape(type='line', x0=-0.5, y0=interval_y, x1=len(filter_n_x) - 0.5, y1=interval_y,
                              line=dict(color='Black', width=2, dash="dashdot"), xref='x', yref='y')

                # Add to tile the current confidence interval.
                fig.update_layout(title=dict(
                    text=str(plot_name) + "\t(" + str(k_independent_tests) + "," + str(round(interval_y, 5)) + ")"))

                # Mark all results which pass the significance border.
                for it in range(len(set_count)):
                    # Check for significant result.
                    if set_count.iloc[it, 2] > interval_y:
                        fig.add_shape(type='line', x0=it - 0.4, y0=interval_y, x1=it + 0.4, y1=interval_y,
                                      line=dict(color='Green', width=3), xref='x', yref='y')

            li_fig.append(fig)
    else:
        # Filter to the current selected genes. df_sub.columns[0] == "geneID"
        df_stat = df_stat[df_stat[df_stat.columns[0]].isin(select_ids)].reset_index(drop=True)
        for it in list(df_stat.iloc[:, 1].drop_duplicates()):
            # Filter to plot duty
            df_sub = df_stat[df_stat.iloc[:, 1].isin([it])].reset_index(drop=True)
            df_sub.iloc[:, 1] = df_sub.iloc[:, 1].apply(lambda x: 1)

            # Show at least one result
            n = bar_n_show if bar_n_show > 0 else 1
            # Filter to the most represented n items
            filter_n_x = df_sub.iloc[:, [1, 2]].groupby(by=df_sub.columns[2]).count() \
                             .sort_values(by=df_sub.columns[1], ascending=False).iloc[:n].index.tolist()

            # Apply filter to the sub dataframe
            df_sub = df_sub[df_sub.iloc[:, 2].isin(filter_n_x)].reset_index(drop=True)

            if not df_sub.empty:
                fig = px.bar(df_sub, x=df_sub.columns[2], y=df_sub.columns[1], color=df_sub.columns[0], title=str(it),
                             custom_data=df_sub.columns[0], category_orders={df_sub.columns[2]: filter_n_x})

                # Update to plot color scale. Using not common use of 'legendgroup'
                for itt in fig.data:
                    if itt['legendgroup'] in dict_colors:
                        itt['marker']['color'] = dict_colors[itt['legendgroup']]

                    # Update hovertemplate
                    itt['hovertemplate'] = "PointID=%{customdata[0]}<br>AnnotatedID=%{x}<extra></extra>"
                li_fig.append(fig)

    # Manage plot visibilities for a dynamic number of plots.
    max_plots = 10
    enable_view = len(li_fig) * [True]

    # Fill up less plots. Disable display for not present plots in enable_view list.
    while len(li_fig) < max_plots:
        li_fig.append(dash.no_update)
        enable_view.append(False)

    return li_fig[:max_plots] + enable_view[2:]


@app.callback(
    dash.dependencies.Output(app_con.ids.dict_selection, 'data', allow_duplicate=True),
    dash.dependencies.Output(app_con.ids.di_add_list, 'is_open', allow_duplicate=True),
    dash.dependencies.Input(app_con.ids.di_add_list_apply, 'n_clicks'),
    dash.dependencies.State(app_con.ids.filename_cluster, 'data'),
    dash.dependencies.State(app_con.ids.dict_selection, 'data'),
    dash.dependencies.State(app_con.ids.dd_active_group, 'value'),
    dash.dependencies.State(app_con.ids.di_add_list_input, 'value'),
    prevent_initial_call=True
)
def addListToGroup(click, filename, d_data, active_group, add_list):
    """ Handle add by id selections. This is useful if you know a group and want to invest these ids.
    :param click: Not used Int of number button hits. Needed to trigger the function manual.
    :param filename: A string which belongs to the filename with including path and refers to an app_con.data key.
    :param d_data: A dictionary with data of selected genes as well there group name. Gene names are main keys.
    :param active_group: The current active editable group. String format with in "group-name,colorcode" information.
    :param add_list: String with id separated by ',' which will add to the selected group.
    :return: Dictionary with updated d_data based on add_list modifications.
    """
    if not filename or filename not in app_con.data or not active_group or not add_list:
        return dash.no_update, False

    if not d_data:
        d_data = {}

    old_data = str(d_data)
    add_list = add_list.replace(" ", "").replace("\t", ",").replace("\n", ",").replace(",,,", ",").replace(",,", ",")
    add_list = add_list[1:] if add_list[0] == "," else add_list
    add_list = add_list[:-1] if add_list[-1] == "," else add_list

    add_list = add_list.split(",")
    active_group = active_group.split(",")

    df = app_con.data[filename]['geneID'].drop_duplicates()
    df = df[df.isin(add_list)]

    for it in list(df):
        d_data[it] = {'name': active_group[0], 'color': active_group[1]}

    if old_data == str(d_data):
        # We got no changes
        return dash.no_update, False

    # Set new Update
    return d_data, False


# ------------------------------------ Dash client side callback definitions -------------------------------------------
app.clientside_callback(
    """function restylePlot( dot_size, dot_opacity, dot_palette, data_sel, fig){
    /**
    *restylePlot() is needed to update just in time display settings.
    *:param dot_size: Int for scatter plot dot size.
    *:param dot_opacity: Float between 1.0 and 0.0 which defines the dot opacity.
    *:param dot_palette: A sting with the name of a color palette. The name have to be in app_con.color_palettes.
    *:param data_sel: A dictionary with data of selected genes as well there group name. Gene names are main keys.
    *:param fig: Plotly figure of the current 2d plot.
    *:return: A Sting with no function.
    */
    if(! fig){
        throw window.dash_clientside.PreventUpdate;
        return "";
    }

    let color_group = []
    fig.data[0].customdata.forEach(it =>{
        if(data_sel[it[0]] !== undefined){
            color_group.push(data_sel[it[0]].color);
        }else{
            color_group.push("rgba(235,235,250,0.0)");
        }
    })

    if(dot_palette){
        var update = {'marker.size': dot_size,
            'marker.opacity': dot_opacity,
            'marker.line.width': dot_size * 0.25,
            'marker.line.color': [color_group],
            'marker.color': [dot_palette]};
    }else{
        var update = {'marker.size': dot_size,
            'marker.opacity': dot_opacity,
            'marker.line.width': dot_size * 0.25,
            'marker.line.color': [color_group],};
    }

    let scatDiv = document.getElementById('main_plot_2d')
    if(scatDiv !== undefined && scatDiv.children !== undefined && scatDiv.children.length >= 2){
        window.Plotly.restyle(scatDiv.children[1], update);
    }

    scatDiv = document.getElementById('scatter_matrix_3d')
    if(scatDiv !== undefined && scatDiv.children !== undefined && scatDiv.children.length >= 2){
        window.Plotly.restyle(scatDiv.children[1], update);
    }

    update['marker.size'] = dot_size * 0.3;
    scatDiv = document.getElementById('main_plot_3d')
    if(scatDiv !== undefined && scatDiv.children !== undefined && scatDiv.children.length >= 2){
        window.Plotly.restyle(scatDiv.children[1], update);
    }



    throw window.dash_clientside.PreventUpdate;
    return "";
    }""",
    dash.dependencies.Output(app_con.ids.dummy_restyle, 'data', allow_duplicate=True),
    dash.dependencies.Input(app_con.ids.slider_dot_size, 'value'),
    dash.dependencies.Input(app_con.ids.slider_opacity, 'value'),
    dash.dependencies.Input(app_con.ids.active_color_palette, 'data'),
    dash.dependencies.Input(app_con.ids.dict_selection, 'data'),
    dash.dependencies.State(app_con.ids.plot_2d, 'figure'),
    prevent_initial_call=True)

app.clientside_callback("""function openGroupDialog(clicks){
    /**
    *openGroupDialog() open dialog box to create new group.
    *:param clicks: Not used Int  of number button hits. Needed to trigger the function manual.
    *:return: Bool to open the dialog box.
    */    
    return true; 
    }""",
                        dash.dependencies.Output(app_con.ids.di_ng, 'is_open', allow_duplicate=True),
                        dash.dependencies.Input(app_con.ids.bt_di_ng, 'n_clicks'),
                        prevent_initial_call=True)

app.clientside_callback("""function closeGroupDialog(clicks){
    /**
    *closeGroupDialog() close dialog box without create a new group.
    *:param clicks: Not used Int  of number button hits. Needed to trigger the function manual.
    *:return: Bool to close the dialog box.
    */
    return false; 
    }""",
                        dash.dependencies.Output(app_con.ids.di_ng, 'is_open', allow_duplicate=True),
                        dash.dependencies.Input(app_con.ids.di_ng_cancel, 'n_clicks'),
                        prevent_initial_call=True)

app.clientside_callback("""function newGroup(clicks, group_name, group_color, old_groups, old_value){
    /**
    *newGroup() to create new group. Groups are used to sort genes according to a user selection.
    *:param clicks: Not used Int  of number button hits. Needed to trigger the function manual.
    *:group_name: Sting with the new group name.
    *:param group_color: Sting with a hex color code to mark the group members.
    *:param old_groups: Dash dropdown options object with all known groups.
    *:param old_value: The current active editable group.
    *:return: Tuple with three entries. 1st Bool show dialog, 2nd new Dash dropdown options, 3rd new editable group.
    */
    if (! group_name || !group_color){
        alert("Fail to create group. Is is necessary to define a group name and color.");
        return [true, old_groups, old_value];
    }
    let new_group = {'label': group_name, 'value': [group_name, group_color].join(",")};
    if (!old_groups) {
        return [false, [new_group], new_group.value];
    }else{
        for (var i = 0; i < old_groups.length; i++) {
            if(old_groups[i].label === group_name){
                alert("Fail to create group. Every group needs a unique name.");
                return [true, old_groups, old_value];
            }
        }
        old_groups.push(new_group);
        return [false, old_groups, new_group.value];
    }
    alert("Something went wrong! Abort create group.")
    return [false, old_groups, old_value];
    }""",
                        dash.dependencies.Output(app_con.ids.di_ng, 'is_open', allow_duplicate=True),
                        dash.dependencies.Output(app_con.ids.dd_active_group, 'options', allow_duplicate=True),
                        dash.dependencies.Output(app_con.ids.dd_active_group, 'value', allow_duplicate=True),
                        dash.dependencies.Input(app_con.ids.di_ng_create, 'n_clicks'),
                        dash.dependencies.State(app_con.ids.di_ng_name, 'value'),
                        dash.dependencies.State(app_con.ids.di_ng_color, 'value'),
                        dash.dependencies.State(app_con.ids.dd_active_group, 'options'),
                        dash.dependencies.State(app_con.ids.dd_active_group, 'value'),
                        prevent_initial_call=True)

app.clientside_callback("""function deleteGroup(clicks, old_groups, old_value, data_sel){
    /**
    *deleteGroup() to delete the current active Group. Groups are used to sort genes according to a user selection.
    *:param clicks: Not used Int  of number button hits. Needed to trigger the function manual.
    *:param old_groups: Dash dropdown options object with all known groups.
    *:param old_value: The current active editable group.
    *:param data_sel: A dictionary with data of selected genes as well there group name. Gene names are main keys.
    *:return: Tuple with three entries. 1st New dropdown options, 2nd new active editable group, 3rd updated data_sel.
    */
    if (! old_groups || ! old_value) {
        throw window.dash_clientside.PreventUpdate
    }

    let new_groups = [];
    for (var i = 0; i < old_groups.length; i++) {
        if(old_groups[i].value !== old_value){
            new_groups.push(old_groups[i]);
        }
    }

    name_del = old_value.split(",")[0];
    for (const [key, value] of Object.entries(data_sel)) {
        if (value.name == name_del) {
            delete data_sel[key];
        }
    }

    if(new_groups.length > 0){
        return [new_groups, new_groups[0].value, data_sel];
    }
    return [undefined, undefined, data_sel];
    }""",
                        dash.dependencies.Output(app_con.ids.dd_active_group, 'options', allow_duplicate=True),
                        dash.dependencies.Output(app_con.ids.dd_active_group, 'value', allow_duplicate=True),
                        dash.dependencies.Output(app_con.ids.dict_selection, 'data', allow_duplicate=True),
                        dash.dependencies.Input(app_con.ids.bt_delete_group, 'n_clicks'),
                        dash.dependencies.State(app_con.ids.dd_active_group, 'options'),
                        dash.dependencies.State(app_con.ids.dd_active_group, 'value'),
                        dash.dependencies.State(app_con.ids.dict_selection, 'data'),
                        prevent_initial_call=True)

app.clientside_callback(
    """function clickPlot(click_2d, select_2d, click_3d, click_3d_sc, select_3d_sc, active_group, d_data, current_tool){
    /**
    *clickPlot() handles click and select events for all scatter plots. Write the dict_selection.
    *:param click_2d: Object with current clicked data into the 2d plot.
    *:param select_2d: Object with all currently selected data into the 2d scatter plot.
    *:param click_3d: Object with current clicked data into the 3d plot. Select it´s not supported by plotly.
    *:param click_3d_sc: Object with current clicked data into the 3d scatter matrix.
    *:param select_3d_sc: Object with all currently selected data into the 3d scatter matrix.
    *:param active_group: The current active editable group. Sting format with in "group-name,colorcode" information.
    *:param d_data: A dictionary with data of selected genes as well there group name. Gene names are main keys.
    *:param current_tool: The current selected edit tool. A sting with on of the options ['add', 'remove', 'pause'].
    *:return: Dictionary with updated d_data based on modifications.
    */
    if(! click_2d && ! select_2d && ! click_3d && ! click_3d_sc && ! select_3d_sc){
        throw window.dash_clientside.PreventUpdate; return d_data;
    }

    var ctx = dash_clientside.callback_context;
    if(! ctx.triggered[0].value){
        throw window.dash_clientside.PreventUpdate;
    }

    var points = ctx.triggered[0].value.points;

    if (! d_data){d_data = {};}

    let old_data = JSON.stringify(d_data)
    if(current_tool === "pause"){
        throw window.dash_clientside.PreventUpdate;
    }
    else if(current_tool === "add"){
        if (! active_group || ! current_tool){throw window.dash_clientside.PreventUpdate; return d_data;}
        active_group = active_group.split(",")

        for (let i = 0; i < points.length; i++) {
            d_data[points[i].customdata[0]] = {'name': active_group[0], 'color' : active_group[1]};
        }
    }
    else if(current_tool === "remove"){
        if (! active_group || ! current_tool){throw window.dash_clientside.PreventUpdate; return d_data;}
        active_group = active_group.split(",")

        for (let i = 0; i < points.length; i++) {
            if (d_data[points[i].customdata[0]] !== undefined){
                delete d_data[points[i].customdata[0]];
            }
        }
    }

    if (old_data === JSON.stringify(d_data)){
        throw window.dash_clientside.PreventUpdate;
    }

    return d_data;}""",
    dash.dependencies.Output(app_con.ids.dict_selection, 'data'),
    dash.dependencies.Input(app_con.ids.plot_2d, 'clickData'),
    dash.dependencies.Input(app_con.ids.plot_2d, 'selectedData'),
    dash.dependencies.Input(app_con.ids.plot_3d, 'clickData'),
    dash.dependencies.Input(app_con.ids.scatter_3d, 'clickData'),
    dash.dependencies.Input(app_con.ids.scatter_3d, 'selectedData'),
    dash.dependencies.State(app_con.ids.dd_active_group, 'value'),
    dash.dependencies.State(app_con.ids.dict_selection, 'data'),
    dash.dependencies.State(app_con.ids.dd_modify_tool, 'value'),
    prevent_initial_call=True)

app.clientside_callback(
    """function clickPlotAdditional(click_aaa, select_aaa, click_bbb, select_bbb, click_ccc, select_ccc, 
    click_ddd, select_ddd, click_eee, select_eee, click_fff, select_fff, click_ggg, select_ggg,click_hhh, 
    select_hhh,click_iii, select_iii,click_jjj, select_jjj, active_group, d_data, current_tool){
    /**
    *clickPlotAdditional() handles click and select additional bar plots. Write the dict_selection.
    *:param click_***: Object with current clicked data of a additional bar plot.
    *:param select_***: Object with current selected data points of a additional bar plot.
    *:param active_group: The current active editable group. Sting format with in "group-name,colorcode" information.
    *:param d_data: A dictionary with data of selected genes as well there group name. Gene names are main keys.
    *:param current_tool: The current selected edit tool. A sting with on of the options ['add', 'remove', 'pause'].
    *:return: Dictionary with updated d_data based on modifications.
    */
    if(click_aaa && ! select_aaa && ! click_bbb && ! select_bbb && !click_ccc && ! select_ccc && 
    ! click_ddd && ! select_ddd && ! click_eee && ! select_eee && ! click_fff && ! select_fff && 
    ! click_ggg && ! select_ggg && ! click_hhh && ! select_hhh && ! click_iii && ! select_iii && 
    ! click_jjj && ! select_jjj){throw window.dash_clientside.PreventUpdate; return d_data;}

    var ctx = dash_clientside.callback_context;

    if(! ctx.triggered[0].value){
        throw window.dash_clientside.PreventUpdate;
    }

    var points = ctx.triggered[0].value.points;

    if (! d_data){d_data = {};}

    let old_data = JSON.stringify(d_data)
    if(current_tool === "pause"){
        throw window.dash_clientside.PreventUpdate;
    }
    else if(current_tool === "add"){
        if (! active_group || ! current_tool){throw window.dash_clientside.PreventUpdate; return d_data;}
        active_group = active_group.split(",")

        for (let i = 0; i < points.length; i++) {
            d_data[points[i].customdata[0]] = {'name': active_group[0], 'color' : active_group[1]};
        }
    }
    else if(current_tool === "remove"){
        if (! active_group || ! current_tool){throw window.dash_clientside.PreventUpdate; return d_data;}
        active_group = active_group.split(",")

        for (let i = 0; i < points.length; i++) {
            if (d_data[points[i].customdata[0]] !== undefined){
                delete d_data[points[i].customdata[0]];
            }
        }
    }

    if (old_data === JSON.stringify(d_data)){
        throw window.dash_clientside.PreventUpdate;
    }

    return d_data;}""",
    dash.dependencies.Output(app_con.ids.dict_selection, 'data', allow_duplicate=True),
    dash.dependencies.Input(app_con.ids.bar_plot_aaa, 'clickData'),
    dash.dependencies.Input(app_con.ids.bar_plot_aaa, 'selectedData'),
    dash.dependencies.Input(app_con.ids.bar_plot_bbb, 'clickData'),
    dash.dependencies.Input(app_con.ids.bar_plot_bbb, 'selectedData'),
    dash.dependencies.Input(app_con.ids.bar_plot_ccc, 'clickData'),
    dash.dependencies.Input(app_con.ids.bar_plot_ccc, 'selectedData'),
    dash.dependencies.Input(app_con.ids.bar_plot_ddd, 'clickData'),
    dash.dependencies.Input(app_con.ids.bar_plot_ddd, 'selectedData'),
    dash.dependencies.Input(app_con.ids.bar_plot_eee, 'clickData'),
    dash.dependencies.Input(app_con.ids.bar_plot_eee, 'selectedData'),
    dash.dependencies.Input(app_con.ids.bar_plot_fff, 'clickData'),
    dash.dependencies.Input(app_con.ids.bar_plot_fff, 'selectedData'),
    dash.dependencies.Input(app_con.ids.bar_plot_ggg, 'clickData'),
    dash.dependencies.Input(app_con.ids.bar_plot_ggg, 'selectedData'),
    dash.dependencies.Input(app_con.ids.bar_plot_hhh, 'clickData'),
    dash.dependencies.Input(app_con.ids.bar_plot_hhh, 'selectedData'),
    dash.dependencies.Input(app_con.ids.bar_plot_iii, 'clickData'),
    dash.dependencies.Input(app_con.ids.bar_plot_iii, 'selectedData'),
    dash.dependencies.Input(app_con.ids.bar_plot_jjj, 'clickData'),
    dash.dependencies.Input(app_con.ids.bar_plot_jjj, 'selectedData'),
    dash.dependencies.State(app_con.ids.dd_active_group, 'value'),
    dash.dependencies.State(app_con.ids.dict_selection, 'data'),
    dash.dependencies.State(app_con.ids.dd_modify_tool, 'value'),
    prevent_initial_call=True)

app.clientside_callback("""function closeAdditionalStat(clicks, active_tab, open_stat){
    /**
    *closeAdditionalStat() manage sidebar display and hide option.
    *:param click: Not used Int  of number button hits. Needed to trigger the function manual.
    *:param active_tab: String with name of the current active tab. Native hide for phyloprofile.
    *:param open_stat: Bool if the user select view or hide sidebar.
    *:return: Tuple with 5 Elements. First Bool to save the user selection, Second Bool to hide or show sidebar.
    *         3rd and 4th restyle dictionaries. 5th string to restyle the trigger button.
    */
    const getStyle = (w) => {return {'width': w, 'flexGrow': 0, 'flexShrink': 0, 'flexBasis': 'auto', 'paddingInline': '0.5vw'};}
    let ctx = dash_clientside.callback_context;
    console.log(ctx)

    if(ctx.triggered[0].prop_id == "tabs_show_data.active_tab"){
        if (active_tab == 'tab_show_phyloprofile'){
            return [false, open_stat, getStyle("95vw"), getStyle("3vw"), "<"]; 
        }
        else{
            if (open_stat){
                return [true, true, getStyle("70vw"), getStyle("28vw"), ">"];

            }
            else{
                return [false, false, getStyle("95vw"), getStyle("3vw"), "<"]; 
            }
        }    
    }

    if (ctx.triggered[0].prop_id == "bt_unfold_statistics.n_clicks")
    {
        if (open_stat){
            return [false, false, getStyle("95vw"), getStyle("3vw"), "<"]; 
        }
        else{
            return [true, true, getStyle("70vw"), getStyle("28vw"), ">"]; 
        }
    }
    console.log("Error closeAdditionalStat()! Unexpected trigger.")
    throw window.dash_clientside.PreventUpdate;
    }""",
                        dash.dependencies.Output(app_con.ids.coll_stat, 'is_open', allow_duplicate=True),
                        dash.dependencies.Output(app_con.ids.store_coll_stat, 'data'),
                        dash.dependencies.Output(app_con.ids.column_gene_viewer, "style"),
                        dash.dependencies.Output(app_con.ids.column_stat, "style"),
                        dash.dependencies.Output(app_con.ids.bt_unfold_stat, 'children'),
                        dash.dependencies.Input(app_con.ids.bt_unfold_stat, 'n_clicks'),
                        dash.dependencies.Input(app_con.ids.tabs_show_data, 'active_tab'),
                        dash.dependencies.State(app_con.ids.store_coll_stat, 'data'),
                        prevent_initial_call=True)

app.clientside_callback("""function openAddListDialog(clicks){
    /**
    *openAddListDialog() open dialog box to add new group members form a list.
    *:param clicks: Not used Int  of number button hits. Needed to trigger the function manual.
    *:return: Bool to open the dialog box.
    */
    return [true, null]; 
    }""",
                        dash.dependencies.Output(app_con.ids.di_add_list, 'is_open', allow_duplicate=True),
                        dash.dependencies.Output(app_con.ids.di_add_list_input, 'value'),
                        dash.dependencies.Input(app_con.ids.bt_add_list, 'n_clicks'),
                        prevent_initial_call=True)

app.clientside_callback("""function closeAddListDialog(clicks){
    /**
    *closeAddListDialog() close dialog box without create adding new group members.
    *:param clicks: Not used Int  of number button hits. Needed to trigger the function manual.
    *:return: Bool to close the dialog box.
    */
    return false; 
    }""",
                        dash.dependencies.Output(app_con.ids.di_add_list, 'is_open', allow_duplicate=True),
                        dash.dependencies.Input(app_con.ids.di_add_list_cancel, 'n_clicks'),
                        prevent_initial_call=True)

app.clientside_callback("""function openAdvancedOptions(value){
    /**
    *openAdvancedOptions() open and close advanced option input.
    *:param clicks: Not used Int  of number button hits. Needed to trigger the function manual.
    *:return: Bool to close the dialog box.
    */
    return value; 
    }""",
                        dash.dependencies.Output(app_con.ids.coll_adv_opt, 'is_open'),
                        dash.dependencies.Input(app_con.ids.bt_adv_opt, 'value'),
                        prevent_initial_call=True)

app.clientside_callback("""function collectStatSelection(click_sel_table, click_2d, select_2d, click_3d, click_3d_sc, 
    select_3d_sc,data_sel_table, filename_cluster){
    /** 
    *collectStatSelection() Collect data points to evaluate the additional statistic.
    *:param click_sel_table: Object with current clicked cell into lower data table. 
    *:param click_2d: Object with current clicked data into the 2d plot.
    *:param select_2d: Object with all currently selected data into the 2d scatter plot.
    *:param click_3d: Object with current clicked data into the 3d plot. Select it´s not supported by plotly.
    *:param click_3d_sc: Object with current clicked data into the 3d scatter matrix.
    *:param select_3d_sc: Object with all currently selected data into the 3d scatter matrix.
    *:param data_sel_table: Main data of the plotly table.
    *:param filename_cluster: A string which belongs to the filename with including path and refers to app_con.data.
    *:return: List with all ids of the selected items.
    */
    if(! click_sel_table && ! click_2d && ! select_2d && ! click_3d && ! click_3d_sc && ! select_3d_sc && 
        ! data_sel_table){
        throw window.dash_clientside.PreventUpdate;
    }

    var ctx = dash_clientside.callback_context;

    if(! ctx.triggered[0] || ! ctx.triggered[0].value || ! filename_cluster){
        throw window.dash_clientside.PreventUpdate;
    }

    let select_ids = null; 

    if(ctx.triggered[0].prop_id  === 'main_data_table.active_cell'){
        var cell = ctx.triggered[0].value;
        if (!cell.hasOwnProperty('column_id') || !cell.hasOwnProperty('row') || cell.column_id !== 'geneID' || 
        !data_sel_table) {
            throw window.dash_clientside.PreventUpdate;
        }
        select_ids = data_sel_table[parseInt(cell.row)][cell.column_id].split(",");
    }
    else if (ctx.triggered[0].value.points){
        select_ids = [];
        ctx.triggered[0].value.points.forEach( it => {
            if (it.customdata[0]) {
                select_ids.push(it.customdata[0]);
            }
        });   
    }
    else {
        throw window.dash_clientside.PreventUpdate;
    }

    return select_ids;
}""",
                        dash.dependencies.Output(app_con.ids.store_sel_stat, 'data'),
                        dash.dependencies.Input(app_con.ids.table_selection, 'active_cell'),
                        dash.dependencies.Input(app_con.ids.plot_2d, 'clickData'),
                        dash.dependencies.Input(app_con.ids.plot_2d, 'selectedData'),
                        dash.dependencies.Input(app_con.ids.plot_3d, 'clickData'),
                        dash.dependencies.Input(app_con.ids.scatter_3d, 'clickData'),
                        dash.dependencies.Input(app_con.ids.scatter_3d, 'selectedData'),
                        dash.dependencies.State(app_con.ids.table_selection, 'data'),
                        dash.dependencies.State(app_con.ids.filename_cluster, 'data'),
                        prevent_initial_call=True)

app.clientside_callback("""function alertUserMessage(user_message){
    /**
    *alertUserMessage() Print alert the given string.
    *:param user_message: Message which will be displayed.
    *:return: null.
    */
    alert(user_message);
    return null; 
    }""",
                        dash.dependencies.Output(app_con.ids.dummy_message, 'data'),
                        dash.dependencies.Input(app_con.ids.user_message, 'data'),
                        prevent_initial_call=True)
# -------------------------------------- Dash backend callback definitions ---------------------------------------------
if __name__ == "__main__":
    # Clean last temp files. This is a workaround, to avoid not deleted temp files fill the storage.
    temp_dir = os.path.dirname(os.path.abspath(__file__)).replace('\\', '/') + app_con.const.temp_folder
    for it in os.listdir(temp_dir):
        os.remove(temp_dir + '/' + it)

    app.run_server(host='127.0.0.1', port='8050', debug=True)
