"""Code that evaluates the quality of the MVM products generated by the drizzlepac package.

Visualization of these Pandas DataFrames with Bokeh can follow the example
from:

https://programminghistorian.org/en/lessons/visualizing-with-bokeh

PUBLIC INTERFACE FOR THIS MODULE
 - build_mvm_plots(HDF5_FILE, output_basename='', display_plot=False):

python mvm_quality_graphics.py mvm_qa_dataframe.h5

"""

# Standard library imports
import argparse
import json
import os
import random
import sys

# Related third party imports
import pandas as pd
from bokeh.layouts import row, column, gridplot
from bokeh.plotting import figure, output_file, save, show
from bokeh.models import ColumnDataSource, Label, CDSView, Div

# Local application imports
from drizzlepac.haputils.graph_utils import HAPFigure, build_tooltips
from drizzlepac.haputils.pandas_utils import get_pandas_data
from stsci.tools import logutil

__taskname__ = 'mvm_quality_graphics'

MSG_DATEFMT = '%Y%j%H%M%S'
SPLUNK_MSG_FORMAT = '%(asctime)s %(levelname)s src=%(name)s- %(message)s'
log = logutil.create_logger(__name__, level=logutil.logging.NOTSET, stream=sys.stdout,
                            format=SPLUNK_MSG_FORMAT, datefmt=MSG_DATEFMT)


# ------------------------------------------------------------------------------------------------------------

def build_mvm_plots(data_source, output_basename='mvm_qa', display_plot=False, log_level=logutil.logging.INFO):
    """Create all the plots for the results generated by these comparisons

    Parameters
    ----------
    data_source : str
        name of the .h5 file produced by diagnostic_json_harvester.py that holds the information to plot.

    output_basename : str
        text string that will be used as the basis for all .html files generated by this script. Unless
        explicitly specified, the default value is 'mvm_qa'. If a value is explicitly specified, the text
        string '_mvm_qa' will be automatically appended to the end of the user-specified value.

    display_plot : bool, optional
        If set to Boolean 'True', plots .html files will be automatically opened in the default web browser
        as they are generated. Unless explicitly specified, the default value is Boolean 'False', meaning
        that plots will only be written to output .html files but not displayed on-screen.

    log_level : int, optional
        The desired level of verboseness in the log statements displayed on the screen and written to the
        .log file. Default value is 20, or 'info'.

    Returns
    -------
    Nothing.
    """
    log.setLevel(log_level)
    if output_basename == '':
        output_basename = "mvm_qa"
    else:
        output_basename = "{}_mvm_qa".format(output_basename)

    # Generate overlap crossmatch plots
    try:
        build_overlap_crossmatch_plots(data_source,
                                       display_plot,
                                       output_basename=output_basename,
                                       log_level=log_level)
    except Exception:
        log.warning("Overlap crossmatch plot generation encountered a problem.")
        log.exception("message")
        log.warning("Continuing to next plot...")
#     -      -     -      -     -      -     -      -     -      -     -      -     -      -     -      -


# ------------------------------------------------------------------------------------------------------------


def build_overlap_crossmatch_plots(data_source, display_plot=False, output_basename='mvm_qa', log_level=logutil.logging.INFO):
    """retrieve required data and reformat the dataframe in preperation for the generation of overlap
    crossmatch plots.

    Parameters
    ----------
    data_source : str
        name of the .h5 file produced by diagnostic_json_harvester.py that holds the information to plot.

    display_plot : bool, optional.
        If set to Boolean 'True', plots .html files will be automatically opened in the default web browser
        as they are generated. Unless explicitly specified, the default value is Boolean 'False', meaning
        that plots will only be written to output .html files but not displayed on-screen.

    output_basename : str, optional.
        text string that will be used as the basis for all .html files generated by this script. Unless
        explicitly specified, the default value is 'mvm_qa'.

    log_level : int, optional
        The desired level of verboseness in the log statements displayed on the screen and written to the
        .log file. Default value is 20, or 'info'.

    Returns
    -------
    Nothing
    """
    log.setLevel(log_level)
    # lists of column titles that will be used.
    details_column_basenames = ["overlap_region_size",
                                "reference_catalog_name",
                                "comparison_catalog_name",
                                "total_reference_catalog_size",
                                "total_comparison_catalog_size",
                                "number_of_reference_sources_available_for_crossmatch",
                                "number_of_comparison_sources_available_for_crossmatch",
                                "number_of_crossmatched_sources"]

    difference_column_basenames = ["X-axis_differences",
                                   "Y-axis_differences",
                                   "On-sky_separation_(X-Y)",
                                   "On-sky_separation_(RA-Dec)"]

    stats_column_basenames = ["Non-clipped_minimum",
                              "Non-clipped_maximum",
                              "Non-clipped_mean",
                              "Non-clipped_median",
                              "Non-clipped_standard_deviation",
                              "3x3-sigma_clipped_mean",
                              "3x3-sigma_clipped_median",
                              "3x3-sigma_clipped_standard_deviation",
                              "Percent_of_all_diff_values_within_1-sigma_of_the_clipped_mean",
                              "Percent_of_all_diff_values_within_2-sigma_of_the_clipped_mean",
                              "Percent_of_all_diff_values_within_3-sigma_of_the_clipped_mean"]

    data_table_column_basename = "crossmatched_reference_X,_Y,_RA,_Dec,_and_crossmatched_comparison_-_reference_difference_values"
    data_table_colnames = ["X-Skycell",
                           "Y-Skycell",
                           "RA",
                           "DEC",
                           "X-axis differences",
                           "Y-axis differences",
                           "On-sky separation (X-Y)",
                           "On-sky separation (RA-Dec)"]
    color_list = ["black", "blue", "brown", "fuchsia", "gold", "green", "olive", "orange", "purple",
                  "rebeccapurple", "red", "rosybrown", "royalblue", "saddlebrown", "salmon", "sandybrown",
                  "seagreen", "springgreen", "steelblue", "tan", "teal", "thistle", "tomato", "turquoise",
                  "violet", "wheat"]

    n_layers_colname = 'gen_info.number of overlap regions present'
    num_layers = get_pandas_data(data_source, [n_layers_colname])[n_layers_colname]

    # retrieve relevant data and "restack" and rename dataframe so that all the information for each overlap
    # region is stored in discrete rows
    restacked_overlap_dataframe = pd.DataFrame()
    for df_indexname, layer_val in zip(num_layers.index.values, num_layers.values):
        for layer_ctr in range(1, layer_val + 1):
            columns_to_retrieve = []
            column_basename = "overlap_region_#{}".format(layer_ctr)
            # add "overlap details" columns
            for details_colname in details_column_basenames:
                columns_to_retrieve.append("{}_details.{}".format(column_basename, details_colname))
            # add stats columns for each difference type
            for diff_type in difference_column_basenames:
                for stats_colname in stats_column_basenames:
                    columns_to_retrieve.append("{}_{}.{}".format(column_basename, diff_type, stats_colname))
            # add all the data table columns
            for data_table_colname in data_table_colnames:
                columns_to_retrieve.append(
                    "{}_{}.{}".format(column_basename, data_table_column_basename, data_table_colname))
            overlap_dataframe = get_pandas_data(data_source, columns_to_retrieve)
            overlap_dataframe = overlap_dataframe[overlap_dataframe['gen_info.dataframe_index'] == df_indexname]
            overlap_dataframe['gen_info.dataframe_index'] = "{}_overlap_region_{}".format(overlap_dataframe['gen_info.dataframe_index'][0], layer_ctr)
            col_rename_dict = {}
            for colname in columns_to_retrieve:
                col_rename_dict[colname] = colname.replace(column_basename, "overlap_region")
            overlap_dataframe = overlap_dataframe.rename(columns=col_rename_dict)
            overlap_dataframe['colormap'][0] = random.sample(color_list, 1)[0]  # assign each row a random color for the plot
            restacked_overlap_dataframe = restacked_overlap_dataframe.append(overlap_dataframe)
    # Sort columns alphabetically to make it more human-friendly
    restacked_overlap_dataframe = restacked_overlap_dataframe[overlap_dataframe.columns.sort_values()]
    # optionally write dataframe to .csv file.
    if log_level == logutil.logging.DEBUG:
        output_csv_filename = "testout.csv"
        if os.path.exists(output_csv_filename):
            os.remove(output_csv_filename)
        restacked_overlap_dataframe.to_csv(output_csv_filename)
        log.debug("Wrote restacked dataframe to csv file {}".format(output_csv_filename))

    # generate plots!
    generate_overlap_crossmatch_graphics(restacked_overlap_dataframe,
                                         display_plot=display_plot,
                                         output_basename=output_basename,
                                         log_level=log_level)


# ------------------------------------------------------------------------------------------------------------

def generate_overlap_crossmatch_graphics(dataframe, display_plot=False, output_basename='mvm_qa', log_level=logutil.logging.INFO):
    """Generate plots to statistically quantify the quality of the alignment of crossmatched sources in
    regions where observations from different proposal/visits overlap.

    Parameters
    ----------
    dataframe : pandas dataframe
        dataframe containing results from the overlap crossmatch(s) to plot.

    display_plot : bool, optional.
        If set to Boolean 'True', plots .html files will be automatically opened in the default web browser
        as they are generated. Unless explicitly specified, the default value is Boolean 'False', meaning
        that plots will only be written to output .html files but not displayed on-screen.

    output_basename : str, optional.
        text string that will be used as the basis for all .html files generated by this script. Unless
        explicitly specified, the default value is 'mvm_qa'.

    log_level : int, optional
        The desired level of verboseness in the log statements displayed on the screen and written to the
        .log file. Default value is 20, or 'info'.

    Returns
    -------
    Nothing
    """
    log.setLevel(log_level)
    xmatch_cds = ColumnDataSource(dataframe)
    # generate plots of x vs. y components of various stat. measures for each difference
    output = "{}_overlap_crossmatch_stats_plots".format(output_basename)
    if not output.endswith('.html'):
        output = output + '.html'
    # Set the output file immediately as advised by Bokeh.
    output_file(output)

    # Generate the graphic-specific tooltips - be mindful of the default tooltips defined in graph_utils.py
    hover_columns = ["gen_info.dataframe_index",
                     "overlap_region_details.overlap_region_size",
                     "overlap_region_details.reference_catalog_name",
                     "overlap_region_details.comparison_catalog_name",
                     "overlap_region_details.number_of_crossmatched_sources"]
    tooltips_list = ['Overlap region name',
                     'Overlap region size (pixels)',
                     'Ref. SVM catalog name',
                     'Comp. SVM catalog name',
                     'Number of crossmatched sources']
    tooltips_list_dynamic = ["X value", "Y value"]
    # hover_tips = build_tooltips(tooltips_list, hover_columns, list(range(0, len(hover_columns))))
    # Define the graphics
    plot_list = []
    # Create title text at the top of the html file
    html_title_text = Div(text="""<h1>Distribution characteristics of crossmatched sources identified in regions of overlapping observations in the MVM product</h1>""")
    plot_list.append(html_title_text)
    # Scatter plots!
    hover_columns_dynamic = ['overlap_region_X-axis_differences.Non-clipped_minimum',
                             'overlap_region_Y-axis_differences.Non-clipped_minimum']
    hover_tips = build_tooltips(tooltips_list+tooltips_list_dynamic, hover_columns+hover_columns_dynamic,
                                list(range(0, len(hover_columns)+2)))
    p0 = HAPFigure(title='Minimum difference value',
                   x_label='X minimum difference (pixels)',
                   y_label='Y minimum difference (pixels)', hover_tips=hover_tips)
    p0.build_glyph('circle',
                   x='overlap_region_X-axis_differences.Non-clipped_minimum',
                   y='overlap_region_Y-axis_differences.Non-clipped_minimum',
                   sourceCDS=xmatch_cds,
                   glyph_color='colormap',
                   legend_group='gen_info.dataframe_index')

    hover_columns_dynamic = ['overlap_region_X-axis_differences.Non-clipped_maximum',
                             'overlap_region_Y-axis_differences.Non-clipped_maximum']
    hover_tips = build_tooltips(tooltips_list+tooltips_list_dynamic, hover_columns+hover_columns_dynamic,
                                list(range(0, len(hover_columns)+2)))
    p1 = HAPFigure(title='Maximum difference value',
                   x_label='X maximum difference (pixels)',
                   y_label='Y maximum difference (pixels)', hover_tips=hover_tips)
    p1.build_glyph('circle',
                   x='overlap_region_X-axis_differences.Non-clipped_maximum',
                   y='overlap_region_Y-axis_differences.Non-clipped_maximum',
                   sourceCDS=xmatch_cds,
                   glyph_color='colormap',
                   legend_group='gen_info.dataframe_index')
    row1 = row(p0.fig, p1.fig)
    plot_list.append(row1)

    hover_columns_dynamic = ['overlap_region_X-axis_differences.Non-clipped_median',
                             'overlap_region_Y-axis_differences.Non-clipped_median']
    hover_tips = build_tooltips(tooltips_list+tooltips_list_dynamic, hover_columns+hover_columns_dynamic,
                                list(range(0, len(hover_columns)+2)))
    p2 = HAPFigure(title='Median difference value',
                   x_label='X median difference (pixels)',
                   y_label='Y median difference (pixels)', hover_tips=hover_tips)
    p2.build_glyph('circle',
                   x='overlap_region_X-axis_differences.Non-clipped_median',
                   y='overlap_region_Y-axis_differences.Non-clipped_median',
                   sourceCDS=xmatch_cds,
                   glyph_color='colormap',
                   legend_group='gen_info.dataframe_index')

    hover_columns_dynamic = ['overlap_region_X-axis_differences.Non-clipped_mean',
                             'overlap_region_Y-axis_differences.Non-clipped_mean']
    hover_tips = build_tooltips(tooltips_list+tooltips_list_dynamic, hover_columns+hover_columns_dynamic,
                                list(range(0, len(hover_columns)+2)))
    p3 = HAPFigure(title='Mean difference value',
                   x_label='X mean difference (pixels)',
                   y_label='Y mean difference (pixels)', hover_tips=hover_tips)
    p3.build_glyph('circle',
                   x='overlap_region_X-axis_differences.Non-clipped_mean',
                   y='overlap_region_Y-axis_differences.Non-clipped_mean',
                   sourceCDS=xmatch_cds,
                   glyph_color='colormap',
                   legend_group='gen_info.dataframe_index')
    row2 = row(p2.fig, p3.fig)
    plot_list.append(row2)

    hover_columns_dynamic = ['overlap_region_X-axis_differences.Non-clipped_standard_deviation',
                             'overlap_region_Y-axis_differences.Non-clipped_standard_deviation']
    hover_tips = build_tooltips(tooltips_list+tooltips_list_dynamic, hover_columns+hover_columns_dynamic,
                                list(range(0, len(hover_columns)+2)))
    p4 = HAPFigure(title='Difference standard deviation value',
                   x_label='X standard deviation difference (pixels)',
                   y_label='Y standard deviation difference (pixels)', hover_tips=hover_tips)
    p4.build_glyph('circle',
                   x='overlap_region_X-axis_differences.Non-clipped_standard_deviation',
                   y='overlap_region_Y-axis_differences.Non-clipped_standard_deviation',
                   sourceCDS=xmatch_cds,
                   glyph_color='colormap',
                   legend_group='gen_info.dataframe_index')

    hover_columns_dynamic = ['overlap_region_X-axis_differences.3x3-sigma_clipped_median',
                             'overlap_region_Y-axis_differences.3x3-sigma_clipped_median']
    hover_tips = build_tooltips(tooltips_list+tooltips_list_dynamic, hover_columns+hover_columns_dynamic,
                                list(range(0, len(hover_columns)+2)))
    p5 = HAPFigure(title='3x3 sigma-clipped median difference value',
                   x_label='X sigma-clipped median difference (pixels)',
                   y_label='Y sigma-clipped median difference (pixels)', hover_tips=hover_tips)
    p5.build_glyph('circle',
                   x='overlap_region_X-axis_differences.3x3-sigma_clipped_median',
                   y='overlap_region_Y-axis_differences.3x3-sigma_clipped_median',
                   sourceCDS=xmatch_cds,
                   glyph_color='colormap',
                   legend_group='gen_info.dataframe_index')
    row3 = row(p4.fig, p5.fig)
    plot_list.append(row3)

    hover_columns_dynamic = ['overlap_region_X-axis_differences.3x3-sigma_clipped_mean',
                             'overlap_region_Y-axis_differences.3x3-sigma_clipped_mean']
    hover_tips = build_tooltips(tooltips_list+tooltips_list_dynamic, hover_columns+hover_columns_dynamic,
                                list(range(0, len(hover_columns)+2)))
    p6 = HAPFigure(title='3x3 sigma-clipped mean difference value',
                   x_label='X sigma-clipped mean difference (pixels)',
                   y_label='Y sigma-clipped mean difference (pixels)', hover_tips=hover_tips)
    p6.build_glyph('circle',
                   x='overlap_region_X-axis_differences.3x3-sigma_clipped_mean',
                   y='overlap_region_Y-axis_differences.3x3-sigma_clipped_mean',
                   sourceCDS=xmatch_cds,
                   glyph_color='colormap',
                   legend_group='gen_info.dataframe_index')

    hover_columns_dynamic = ['overlap_region_X-axis_differences.3x3-sigma_clipped_standard_deviation',
                             'overlap_region_Y-axis_differences.3x3-sigma_clipped_standard_deviation']
    hover_tips = build_tooltips(tooltips_list+tooltips_list_dynamic, hover_columns+hover_columns_dynamic,
                                list(range(0, len(hover_columns)+2)))

    p7 = HAPFigure(title='3x3 sigma-clipped Difference standard deviation value',
                   x_label='X sigma-clipped difference standard_deviation (pixels)',
                   y_label='Y sigma-clipped difference standard_deviation (pixels)', hover_tips=hover_tips)
    p7.build_glyph('circle',
                   x='overlap_region_X-axis_differences.3x3-sigma_clipped_standard_deviation',
                   y='overlap_region_Y-axis_differences.3x3-sigma_clipped_standard_deviation',
                   sourceCDS=xmatch_cds,
                   glyph_color='colormap',
                   legend_group='gen_info.dataframe_index')
    row4 = row(p6.fig, p7.fig)
    plot_list.append(row4)

    if display_plot:
        show(column(plot_list))
    # Just save
    else:
        save(column(plot_list))
    log.info("Output HTML graphic file {} has been written.\n".format(output))

    # generate quad resid (x vs. dx, y vs. dx, x vs. dy, y vs. dy) plots for each DF row

    for dfindex, dfrow in dataframe.iterrows():
        qr_df = pd.DataFrame()
        qr_df = qr_df.append(dfrow)
        qr_cds = ColumnDataSource(qr_df)
        qr_dict = {"x": qr_cds.data['overlap_region_crossmatched_reference_X,_Y,_RA,_Dec,_and_crossmatched_comparison_-_reference_difference_values.X-Skycell'][0],
                   "y": qr_cds.data['overlap_region_crossmatched_reference_X,_Y,_RA,_Dec,_and_crossmatched_comparison_-_reference_difference_values.Y-Skycell'][0],
                   "dx": qr_cds.data['overlap_region_crossmatched_reference_X,_Y,_RA,_Dec,_and_crossmatched_comparison_-_reference_difference_values.X-axis differences'][0],
                   "dy": qr_cds.data['overlap_region_crossmatched_reference_X,_Y,_RA,_Dec,_and_crossmatched_comparison_-_reference_difference_values.Y-axis differences'][0]}
        qr_cds = ColumnDataSource(qr_dict)
        output = "{} {}_overlap_crossmatch_residuals_plots".format(output_basename, qr_df['gen_info.dataframe_index'].values[0])
        if not output.endswith('.html'):
            output = output + '.html'
        # Set the output file immediately as advised by Bokeh.
        output_file(output)
        # Define plots
        plot_list = []
        html_title_text = Div(text="""<h1>Crossmatched comparison - reference residuals:<br>{}</h1>""".format(qr_df['gen_info.dataframe_index'].values[0]))
        plot_list.append(html_title_text)
        # add descriptive info
        for detail_title, detail_value in zip(tooltips_list[1:], hover_columns[1:]):
            if detail_title in ["Overlap region size (pixels)", "Number of crossmatched sources"]:
                dv = int(qr_df[detail_value].values[0])
            else:
                dv = qr_df[detail_value].values[0]
            detail_html_text = Div(text="""<h3>{}: {}</h3>""".format(detail_title, dv))
            plot_list.append(detail_html_text)
        p1 = HAPFigure(title='X vs DX',
                       x_label="X (pixels)",
                       y_label='Delta[X] (pixels)',
                       use_hover_tips=False)
        p1.build_glyph('circle',
                       x='x',
                       y='dx',
                       sourceCDS=qr_cds)

        p2 = HAPFigure(title='X vs DY',
                       x_label="X (pixels)",
                       y_label='Delta[Y] (pixels)',
                       use_hover_tips=False)
        p2.build_glyph('circle',
                       x='x',
                       y='dy',
                       sourceCDS=qr_cds)
        row1 = row(p1.fig, p2.fig)
        plot_list.append(row1)

        p3 = HAPFigure(title='Y vs DX',
                       x_label="Y (pixels)",
                       y_label='Delta[X] (pixels)',
                       use_hover_tips=False)
        p3.build_glyph('circle',
                       x='y',
                       y='dx',
                       sourceCDS=qr_cds)

        p4 = HAPFigure(title='Y vs DY',
                       x_label="Y (pixels)",
                       y_label='Delta[Y] (pixels)',
                       use_hover_tips=False)
        p4.build_glyph('circle',
                       x='y',
                       y='dy',
                       sourceCDS=qr_cds)
        row2 = row(p3.fig, p4.fig)
        plot_list.append(row2)

        # Display and save
        if display_plot:
            show(column(plot_list))
        # Just save
        else:
            save(column(plot_list))
        log.info("Output HTML graphic file {} has been written.\n".format(output))


# ------------------------------------------------------------------------------------------------------------

if __name__ == "__main__":
    """Simple command-line interface. That is all.
    """
    # process command-line inputs with argparse
    parser = argparse.ArgumentParser(description='Generate MVM quality assessment plots based on information'
                                                 ' stored in the user-specified .h5 file.')
    parser.add_argument('input_filename',
                        help='Input .h5 file produced by diagnostic_json_harvester.py that holds the '
                             'information to plot.')
    parser.add_argument('-d', '--display_plot', required=False, action='store_true',
                        help='If specified, plots will be automatically opened in the default web browser as '
                             'they are generated. Otherwise, .html plot files will be generated but not '
                             'opened.')
    parser.add_argument('-l', '--log_level', required=False, default='info',
                        choices=["critical", "error", "warning", "info", "debug", "notset"],
                        help='The desired level of verboseness in the log statements displayed on the screen '
                             'and written to the .log file.')
    parser.add_argument('-o', '--output_basename', required=False, default='',
                        help='Base name for the HMTL file generated by Bokeh.')
    user_args = parser.parse_args()

    # verify that input file exists
    if not os.path.exists(user_args.input_filename):
        err_msg = "File {} doesn't exist.".format(user_args.input_filename)
        log.critical(err_msg)
        sys.exit(err_msg)

    # set up logging
    log_dict = {"critical": logutil.logging.CRITICAL,
                "error": logutil.logging.ERROR,
                "warning": logutil.logging.WARNING,
                "info": logutil.logging.INFO,
                "debug": logutil.logging.DEBUG,
                "notset": logutil.logging.NOTSET}
    log_level = log_dict[user_args.log_level]
    log.setLevel(log_level)

    # execute plot generation!
    build_mvm_plots(user_args.input_filename, output_basename=user_args.output_basename,
                    display_plot=user_args.display_plot, log_level=log_level)
