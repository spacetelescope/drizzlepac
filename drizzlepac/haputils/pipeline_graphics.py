"""Module for generating Q&A plots for pipeline processed data"""
import os
import sys

from bokeh.layouts import row, column
from bokeh.plotting import figure, output_file, save
from bokeh.models import ColumnDataSource, Label, Range1d, Paragraph

import numpy as np

from stsci.tools import logutil

from .pandas_utils import PandasDFReader
from . import graph_utils
from .graph_utils import HAPFigure



MSG_DATEFMT = '%Y%j%H%M%S'
SPLUNK_MSG_FORMAT = '%(asctime)s %(levelname)s src=%(name)s- %(message)s'
log = logutil.create_logger(__name__, level=logutil.logging.NOTSET, stream=sys.stdout,
                            format=SPLUNK_MSG_FORMAT, datefmt=MSG_DATEFMT)
                            
__taskname__ = 'pipeline_graphics'

"""
NOTES
-----
Code for generating vectors in Bokeh was found in:

 https://github.com/holoviz/holoviews/blob/master/holoviews/plotting/bokeh/chart.py#L309

This code shows how they define the arrow heads for the vector plots.  A
'non-invasive' way to use the Bokeh Segment Glyph to create a vector plot would
be to define the arrow heads as 2 additional segments starting at the endpoint
of the vector. The new segments would then be added to the data points to be
plotted by Segment to appear as an arrow head on each segment.



# Auto-scaling computation for Quiver plot from matplotlib.
# https://github.com/matplotlib/matplotlib/blob/master/lib/matplotlib/quiver.py#L666-673
 
# self.N : number of points to be plotted
# amean : mean of vector lengths
# self.span : width of plot

            sn = max(10, math.sqrt(self.N))
            if self.Umask is not ma.nomask:
                amean = a[~self.Umask].mean()
            else:
                amean = a.mean()
            # crude auto-scaling
            # scale is typical arrow length as a multiple of the arrow width
            scale = 1.8 * amean * sn / self.span

Example from bokeh.org on how to create a tabbed set of plots:

.. code:: python

    from bokeh.io import output_file, show
    from bokeh.models import Panel, Tabs
    from bokeh.plotting import figure

    output_file("slider.html")

    p1 = figure(plot_width=300, plot_height=300)
    p1.circle([1, 2, 3, 4, 5], [6, 7, 2, 4, 5], size=20, color="navy", alpha=0.5)
    tab1 = Panel(child=p1, title="circle")

    p2 = figure(plot_width=300, plot_height=300)
    p2.line([1, 2, 3, 4, 5], [6, 7, 2, 4, 5], line_width=3, color="navy", alpha=0.5)
    tab2 = Panel(child=p2, title="line")

    tabs = Tabs(tabs=[ tab1, tab2 ])

    show(tabs)


"""

# -------------------------------------------------------------------------------
# Generate the Bokeh plot for the pipeline astrometric data.
#
HOVER_COLUMNS = ['gen_info.imgname',
                 'gen_info.instrument',
                 'gen_info.detector',
                 'gen_info.filter',                
                 'header.DATE-OBS',
                 'header.RA_TARG',
                 'header.DEC_TARG',
                 'header.GYROMODE',
                 'header.EXPTIME',
                 'fit_results.aligned_to']
                   
TOOLTIPS_LIST = ['EXPNAME', 'INSTRUMENT', 'DET', 'FILTER', 
                  'DATE', 'RA', 'DEC', 'GYRO', 'EXPTIME',
                  'ALIGNED_TO']
INSTRUMENT_COLUMN = 'full_instrument'

RESULTS_COLUMNS = ['fit_results.rms_x',
                   'fit_results.rms_y',
                   'fit_results.xsh',
                   'fit_results.ysh',
                   'fit_results.rot',
                   'fit_results.scale',
                   'fit_results.rot_fit',
                   'fit_results.scale_fit',
                   'fit_results.nmatches',
                   'fit_results.skew']

RESIDS_COLUMNS = ['residuals.x',
                  'residuals.y',
                  'residuals.ref_x',
                  'residuals.ref_y']
TOOLSEP_START = '{'
TOOLSEP_END = '}'
    
FIGURE_TOOLS = graph_utils.FIGURE_TOOLS_BASE

def build_vector_plot(sourceCDS, **plot_dict):
    """Create figure object for plotting desired columns as a scatter plot with circles

    Parameters
    ----------
    sourceCDS : Pandas ColumnDataSource
        Object with all the input data

    x_label, y_label : str, optional
        Labels to use for the X and Y axes (respectively)

    title : str, optional
        Title of the plot

    color : string, optional
        Single color to use for data points in the plot if `colormap` is not used

    """
    # Check for optional elements
    x_label = plot_dict.get('x_label', '')
    y_label = plot_dict.get('y_label', '')
    title = plot_dict.get('title', '')
    color = plot_dict.get('color', 'blue')

    # Define a figure object with square aspect ratio
    p1 = HAPFigure(title=title,
                   x_label=x_label,
                   y_label=y_label,
                   use_hover_tips=False)

    # Add the glyphs
    p1.build_vector_glyph(sourceCDS,
                          color=color)

    return p1


def get_pandas_data(pandas_filename):
    """Load the harvested data, stored in a CSV file, into local arrays.

    Parameters
    ==========
    pandas_filename: str
        Name of the CSV file created by the harvester.

    Returns
    =======
    phot_data: Pandas dataframe
        Dataframe which is a subset of the input Pandas dataframe written out as
        a CSV file.  The subset dataframe consists of only the requested columns
        and rows where all of the requested columns did not contain NaNs.

    """
    
    # Instantiate a Pandas Dataframe Reader (lazy instantiation)
    # df_handle = PandasDFReader_CSV("svm_qa_dataframe.csv")
    df_handle = PandasDFReader(pandas_filename, log_level=logutil.logging.NOTSET)

    # In this particular case, the names of the desired columns do not
    # have to be further manipulated, for example, to add dataset specific
    # names.
    # 
    # Get the relevant column data, eliminating all rows which have NaNs
    # in any of the relevant columns.
    if pandas_filename.endswith('.h5'):
        fit_data = df_handle.get_columns_HDF5(HOVER_COLUMNS + RESULTS_COLUMNS)
        resids_data = df_handle.get_columns_HDF5(RESIDS_COLUMNS)
    else:
        fit_data = df_handle.get_columns_CSV(HOVER_COLUMNS + RESULTS_COLUMNS)
        resids_data = df_handle.get_columns_CSV(RESIDS_COLUMNS)

    return fit_data, resids_data



def generate_summary_plots(fitCDS, output='cal_qa_results.html'):
    """Generate the graphics associated with this particular type of data.

    Parameters
    ==========
    fitCDS : ColumnDataSource object
        Dataframe consisting of the relative alignment astrometric fit results
    
    output : str, optional
        Filename for output file with generated plot
         
    Returns
    -------
    output : str
        Name of HTML file where the plot was saved.

    """
    # TODO: include the date from the input data as part of the html filename
    # Set the output file immediately as advised by Bokeh.
    output_file(output)

    num_of_datasets = len(fitCDS.data['index'])
    print('Number of datasets: {}'.format(num_of_datasets))
        
    plot_list = []

    # Data point figures
    p1 = HAPFigure(title='RMS Values',
                   x_label="RMS_X (pixels)",
                   y_label="RMS_Y (pixels)",
                   tools=FIGURE_TOOLS)

    p1.build_glyph('circle', 
                   x=RESULTS_COLUMNS[0], 
                   y=RESULTS_COLUMNS[1],
                   sourceCDS=fitCDS,  
                   glyph_color='colormap',
                   legend_group='inst_det')


    # Data point figures
    p2 = HAPFigure(title='Offsets',
                   x_label="Shift X (pixels)",
                   y_label="Shift Y (pixels)",
                   tools=FIGURE_TOOLS)

    p2.build_glyph('circle', 
                   x=RESULTS_COLUMNS[2], 
                   y=RESULTS_COLUMNS[3],
                   sourceCDS=fitCDS,  
                   glyph_color='colormap',
                   legend_group='inst_det')

    p3 = HAPFigure(title='Rotation',
                   x_label="Number of matched sources",
                   y_label="Rotation (degrees)",
                   tools=FIGURE_TOOLS)

    p3.build_glyph('circle', 
                   x=RESULTS_COLUMNS[8], 
                   y=RESULTS_COLUMNS[4],
                   sourceCDS=fitCDS,  
                   glyph_color='colormap',
                   legend_group='inst_det')

    p4 = HAPFigure(title='Scale',
                   x_label="Number of matched sources",
                   y_label="Scale",
                   tools=FIGURE_TOOLS)

    p4.build_glyph('circle', 
                   x=RESULTS_COLUMNS[8], 
                   y=RESULTS_COLUMNS[5],
                   sourceCDS=fitCDS,  
                   glyph_color='colormap',
                   legend_group='inst_det')

    p5 = HAPFigure(title='Skew',
                   x_label="Number of matched sources",
                   y_label="Skew (degrees)",
                   tools=FIGURE_TOOLS)

    p5.build_glyph('circle', 
                   x=RESULTS_COLUMNS[8], 
                   y=RESULTS_COLUMNS[9],
                   sourceCDS=fitCDS,  
                   glyph_color='colormap',
                   legend_group='inst_det')

    # Save the generated plots to an HTML file define using 'output_file()'
    save(column(p1.fig, p2.fig, p3.fig, p4.fig, p5.fig))
                     
    return output
    
def generate_residual_plots(residsCDS, filename, rms=(0.,0.),
                            display_plot=False, output_dir=None, output=''):
    rootname = '_'.join(filename.split("_")[:-1])
    output = '{}_vectors_{}'.format(rootname, output)
    if output_dir is not None:
        output = os.path.join(output_dir, output)
    output_file(output)

    delta_x = np.array(residsCDS.data['x']) - np.array(residsCDS.data['xr'])
    delta_y = np.array(residsCDS.data['y']) - np.array(residsCDS.data['yr'])
    residsCDS.data['dx'] = delta_x
    residsCDS.data['dy'] = delta_y
    
    npoints = len(delta_x)
    if npoints < 2:
        return None

    residsCDS.data['x'] = np.array(residsCDS.data['x']) - min(residsCDS.data['x'])
    residsCDS.data['y'] = np.array(residsCDS.data['y']) - min(residsCDS.data['y'])

    # Define Y scale for plots that will be used for all residual plots
    
    xmax = np.abs(delta_x).max()
    ymax = np.abs(delta_y).max()
    max_range = round(max(xmax, ymax, 0.4), 1) + 0.1
    plot_range = Range1d(-max_range, max_range)

    rms_x,rms_y = rms
    title_start = 'Residuals for {} '.format(filename)
    title_rms = 'RMS(X)={:.3f}, RMS(Y)={:.3f}'.format(rms_x, rms_y)
    title_base = '{}[{} sources, {}]'.format(title_start, npoints, title_rms)

    p1 = HAPFigure(title='X vs DX',
                   x_label="X (pixels)",
                   y_label='Delta[X] (pixels)',
                   y_range=plot_range,
                   use_hover_tips=False)
    p1.build_glyph('circle',
                   x='x',
                   y='dx',
                   sourceCDS=residsCDS)

    p2 = HAPFigure(title='X vs DY',
                   x_label="X (pixels)",
                   y_label='Delta[Y] (pixels)',
                   y_range=plot_range,
                   use_hover_tips=False)
    p2.build_glyph('circle',
                   x='x',
                   y='dy',
                   sourceCDS=residsCDS)

    p3 = HAPFigure(title='Y vs DX',
                   x_label="Y (pixels)",
                   y_label='Delta[X] (pixels)',
                   y_range=plot_range,
                   use_hover_tips=False)
    p3.build_glyph('circle',
                   x='y',
                   y='dx',
                   sourceCDS=residsCDS)

    p4 = HAPFigure(title='Y vs DY',
                   x_label="Y (pixels)",
                   y_label='Delta[Y] (pixels)',
                   y_range=plot_range,
                   use_hover_tips=False)
    p4.build_glyph('circle',
                   x='y',
                   y='dy',
                   sourceCDS=residsCDS)

    row1 = row(p1.fig, p3.fig)
    
    title_text = Paragraph(text=title_base, width=600, height=40)

    row2 = row(p2.fig, p4.fig)

    pv = build_vector_plot(residsCDS,
                           title='Vector plot of Image - Reference for {}'.format(filename),
                           x_label='X (pixels)',
                           y_label='Y (pixels)',
                           color='blue')

    # Display and save
    if display_plot:
        show(column(row1, title_text, row2, pv.fig))
    else:
        print('Saving: {}'.format(output))
        save(column(row1, title_text, row2, pv.fig))

    return output


def build_astrometry_plots(pandas_file, 
                           output_dir=None, 
                           output='cal_qa_results.html'):

    fit_data, resids_data = get_pandas_data(pandas_file)
    fitCDS = ColumnDataSource(fit_data)
    residsCDS = ColumnDataSource(resids_data)
    

    if output_dir is not None:
        summary_filename = os.path.join(output_dir, output)
    else:
        summary_filename = output
        
    # Generate the astrometric plots
    astrometry_plot_name = generate_summary_plots(fitCDS, output=summary_filename)
    
    resids_plot_names = [astrometry_plot_name]
    
    for filename,rootname in zip(fitCDS.data[HOVER_COLUMNS[0]], fitCDS.data['index']):
        i = residsCDS.data['index'].tolist().index(rootname)
        resids_dict = {'x': residsCDS.data[RESIDS_COLUMNS[0]][i],
                       'y': residsCDS.data[RESIDS_COLUMNS[1]][i],
                       'xr': residsCDS.data[RESIDS_COLUMNS[2]][i],
                       'yr': residsCDS.data[RESIDS_COLUMNS[3]][i]}
        cds = ColumnDataSource(resids_dict)

        rms_x = float(fitCDS.data[RESULTS_COLUMNS[0]][i])
        rms_y = float(fitCDS.data[RESULTS_COLUMNS[1]][i])
        
        resids_plot_name = generate_residual_plots(cds, filename, 
                                                   rms=(rms_x,rms_y),
                                                   output_dir=output_dir,
                                                   output=output)
        resids_plot_names.append(resids_plot_name)

    resids_plot_names = list(filter(lambda a: a != None, resids_plot_names))
    
    return resids_plot_names

