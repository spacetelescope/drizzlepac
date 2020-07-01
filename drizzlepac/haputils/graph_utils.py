""" Definition of Class to handle HAP SVM Bokeh Graphics

"""
import logging
import sys
import os
import traceback
import shutil

from stsci.tools import logutil
from astropy.io import fits
import numpy as np

from bokeh.plotting import figure, output_file, show, save
from bokeh.models import ColumnDataSource, Label, Circle
from bokeh.models.tools import HoverTool

MSG_DATEFMT = '%Y%j%H%M%S'
SPLUNK_MSG_FORMAT = '%(asctime)s %(levelname)s src=%(name)s- %(message)s'
log = logutil.create_logger(__name__, level=logutil.logging.NOTSET, stream=sys.stdout,
                            format=SPLUNK_MSG_FORMAT, datefmt=MSG_DATEFMT)

# Data columns from the Pandas dataframe which are to be used as common display value
# for all hover tooltips.
HOVER_BASIC_TIPS = [('Inst/Det', '@{inst_det}'),
                    ('Dataset', '@{gen_info.dataset}'),
                    ('Filter', '@{gen_info.filter}'),
                    ('ImageName', '@{gen_info.imgname}')]

# Default for Bokeh is (‘pan,wheel_zoom,box_zoom,save,reset,help’)
FIGURE_TOOLS_BASE = 'box_zoom, wheel_zoom, box_select, lasso_select, reset, save'

class HAPFigure:
    """ HAPFigure is a base class whose attributes are common values applicable to all
        HAP figures.

        Careful: grid plot will collect all the tools onto one toolbar, but a tool 
        assciated with only one figure, will only work for that one figure.

        Parameters
        ----------
        x_label, y_label : str
            Labels to use for the X and Y axes (respectively)

        title : str
            Title of the plot

        output_filename: str, optional
            Base name for the output HTML file generated and written to the current directory
            Default: hap_graphic.html

        tools : str, optional
            List of figure tools/controls associated with the plot

        toolbar_location : str, optional
            Location of toolbar

        #hover_tips : list of int, optional
        #    List of indices for the columns from `source` to use as hints 
        #    in the HoverTool

        hover_tips : list of tuples, optional
            List of tuples which are ('DisplayName', ColumnDataSource value)
            This list will be appended to the base set of hover tooltips

        log_level : int
            The desired level of verboseness in the log statements displayed on the screen and written to the .log file.
        
    """
    def __init__(self, **figure_dict):
        # set logging level to user-specified level
        #log.setLevel(log_level)
        #self.log_level = log_level
 
        # Append any user requested tools to the base set 
        user_tools = figure_dict.get('tools', '')
        fig_tools = FIGURE_TOOLS_BASE + ', ' + user_tools
        fig_tool_loc = figure_dict.get('toolbar', 'right')
        self.fig = figure(tools=fig_tools, toolbar_location=fig_tool_loc)
        self.fig.legend.click_policy = figure_dict.get('click_policy', 'hide')
        self.fig.legend.location = figure_dict.get('legend_location', 'center-right')

        # Append any user requested tooltips to the base set 
        user_tips = figure_dict.get('hover_tips', [])
        hover_fig = HoverTool()
        hover_fig.tooltips = HOVER_BASIC_TIPS + user_tips
        self.fig.add_tools(hover_fig)

        # Basic figure attributes
        self.fig.title.text = figure_dict.get('title', '')
        self.fig.xaxis.axis_label = figure_dict.get('x_label', '')
        self.fig.yaxis.axis_label = figure_dict.get('y_label', '')
        self.fig.background_fill_color = figure_dict.get('background_fill_color', 'gainsboro')
        #self.fig.background_fill_color = figure_dict.get('background_fill_color', None)

        self.fig.x_range.start = figure_dict.get('xstart', None)
        self.fig.y_range.start = figure_dict.get('ystart', None)
        self.fig.grid.grid_line_color = figure_dict.get('grid_line_color', None) 

        # These attributes are not really on the figure, but are styling attributes
        # used for the "shape" glyphs.  They may be set to non-default values when 
        # the build_glyph() routine is invoked.
        self.marker_color = 'blue'
        self.marker_size = 10
        self.colormap = True
        #self.legend_group = 
        self.legend_label = None
        self.fill_alpha = 0.5 
        self.line_alpha = 0.5 


    def build_glyph(self, glyph_name, x, y, sourceCDS, **data_dict):
        """Generate a glyph representation of X and Y data.
    
            Parameters
            ----------
            sourceCDS : ColumnDataSource object
                The sourceCDS is an object based upon the Pandas dataframe holding the 
                data of interest where the data can be referenced by column name
       
            x, y : str
                Names of X and Y columns of data in sourceCDS.  The columns represent 
                the X- and Y-axis coordinates for the center of the circle marker.

            name : str, optional
                Indentification of the glyph - useful when using HoverTool where the
                value of the name is a dictionary key as the value itself is variable.
    
            marker_color : string, optional
                Single color to use for data points in the plot if `colormap` is not used
        
            colormap : bool, optional
                Specify whether or not to use a pre-defined set of colors for the points
                derived from the `colormap` column in the input data `source`
        
            marker_size : int, optional
                Size of each marker in the plot
        
            legend_group : str, optional
                If `colormap` is used, this is the name of the column from the input
                data `source` to use for defining the legend for the colors used.  The
                same colors should be assigned to all the same values of data from the 
                column, for example, all 'ACS/WFC' data points from the `instrument`
                column should have a `colormap` column value of `blue`.

            glyph_name : str
                When the actual column being used as the data changes dynamically and
                the assciated value. XXX
        """

        # Set the required attributes
        self.glyph_name = glyph_name
        self.x = x
        self.y = y
        self.sourceCDS = sourceCDS
    
        # This will use the 'colormap' column from 'source' for the colors of 
        # each point.  This column should have been populated by the calling
        # routine. 
        #if colormap:
        #    marker_color = 'colormap'

        # Check for optional attributes and use defaults as necessary.  The
        # avaiable attributes are declared in the constructor.
        self.marker_color = data_dict.get('marker_color', self.marker_color)
        self.marker_size = data_dict.get('marker_size', self.marker_size) 
        self.colormap = data_dict.get('colormap', self.colormap)
        #self.legend_group = data_dict.get('legend_group', self.legend_group)
        self.legend_label = data_dict.get('legend_label', self.legend_label)
        self.fill_alpha = data_dict.get('fill_transparency', self.fill_alpha)
        self.line_alpha = data_dict.get('line_transparency', self.line_alpha)

        #renderer = self.fig.select(name=glyph_name)
        #renderer.selection_glyph = Circle(fill_color = self.marker_color)

        # Dictionary of supported "shape" glyphs.  These are really references to
        # the associated method names.
        glyph_types = {'circle': HAPFigure.build_circle_glyph, 
                       'square': HAPFigure.build_square_glyph,
                       'triangle': HAPFigure.build_triangle_glyph}

        glyph_types[glyph_name](self)


    # "Shape" glyphs
    def build_circle_glyph(self):

        self.fig.circle(x = self.x, 
                        y = self.y, 
                        source = self.sourceCDS, 
                        size = self.marker_size, 
                        color = self.marker_color,
                        #legend_label = self.legend_label,
                        fill_alpha = self.fill_alpha,
                        line_alpha = self.line_alpha,
                        name = self.glyph_name)


    def build_square_glyph(self):

        self.fig.square(x = self.x, 
                        y = self.y, 
                        source = self.sourceCDS, 
                        size = self.marker_size, 
                        color = self.marker_color,
                        #legend_label = self.legend_label,
                        fill_alpha = self.fill_alpha,
                        line_alpha = self.line_alpha,
                        name = self.glyph_name)


    def build_triangle_glyph(self):

        self.fig.triangle(x = self.x, 
                        y = self.y, 
                        source = self.sourceCDS, 
                        size = self.marker_size, 
                        color = self.marker_color,
                        legend_label = self.legend_label,
                        fill_alpha = self.fill_alpha,
                        line_alpha = self.line_alpha,
                        name = self.glyph_name)


    def build_histogram(self, left, right, top, bottom, **data_dict):
        """Generate a glyph representation of X and Y data.

        """

        fill_color = data_dict.get('fill_color', 'gray')
        line_color = data_dict.get('line_color', 'black')
        alpha = data_dict.get('fill_transparency', 1.0)

        self.fig.quad(left = left,
                      right = right,
                      top = top,
                      bottom = bottom,
                      fill_color = fill_color,
                      line_color = line_color, 
                      alpha = alpha)
    

    #def build_vector_plot(self, left, right, top, bottom, **data_dict):
