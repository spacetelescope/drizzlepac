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

from bokeh.plotting import figure
from bokeh.models import ColumnDataSource, Label, Circle
from bokeh.models.tools import HoverTool

MSG_DATEFMT = '%Y%j%H%M%S'
SPLUNK_MSG_FORMAT = '%(asctime)s %(levelname)s src=%(name)s- %(message)s'
log = logutil.create_logger(__name__, level=logutil.logging.NOTSET, stream=sys.stdout,
                            format=SPLUNK_MSG_FORMAT, datefmt=MSG_DATEFMT)

# Data columns from the Pandas dataframe which are to be used as common display value
# for all hover tooltips.
# TO DO: Nice to have WCSNAME for all images.  This has to be added to the JSON writer
# and then the harvester too.
HOVER_BASIC_TIPS = [('ImageName', '@{gen_info.imgname}'),
                    ('Proposal ID', '@{gen_info.proposal_id}'),
                    ('ASN ID', '@{header.ASN_ID}'),
                    ('Filter', '@{gen_info.filter}'),
                    ('Inst/Det', '@{inst_det}')]

# Default figure tools.  Bokeh default is (‘pan,wheel_zoom,box_zoom,save,reset,help’)
FIGURE_TOOLS_BASE = 'box_zoom,wheel_zoom,pan,box_select,lasso_select,reset,save'


def build_tooltips(tips):
    """Return list of tuples for tooltips to use in hover tool.

    Parameters
    ----------
    tips : list
        List of indices for the HOVER_COLUMNS entries to be used as tooltips
        to be included in the hover tool.

    """
    tools = [(TOOLTIPS_LIST[i], '@{}{}{}'.format(
                                TOOLSEP_START,
                                HOVER_COLUMNS[i],
                                TOOLSEP_END)) for i in tips]

    return tools


class HAPFigure:
    """ HAPFigure is a class whose attributes are common values applicable to all
        HAP figures.

        The constructor is used to declare and define the basic "figure" attributes.
        Glyph attributes are also declared and set to default values.

        Careful: grid plot will collect all the tools onto one toolbar, but a tool
        assciated with only one figure, will only work for that one figure.

        Parameters
        ----------
        log_level : int
            The desired level of verboseness in the log statements displayed on the screen and written to the .log file.

        x_label, y_label : str, optional
            Labels to use for the X and Y axes (respectively)

        title : str, optional
            Title of the plot

        background_fill_color : str, optional
            Background color of figure area

        tools : str, optional  DEPRECATED 09 July 2020
            List of figure tools/controls associated with the plot

        toolbar_location : str, optional
            Location of toolbar

        hover_tips : list of tuples, optional
            List of tuples which are ('DisplayName', ColumnDataSource value)
            This list will be appended to the base set of hover tooltips

        xstart, ystart : int, optional
            Starting value of respective figure axis

        x_range, y_range : int, optional
            Range value of respective figure axis

        grid_line_color : str, optional
            Color of grid lines

        These are the characteristics which the developers have found useful thus far.  It
        is expected more arguments may need to be added to support desired functionality.

    """
    def __init__(self, log_level=logutil.logging.NOTSET, **figure_dict):
        # set logging level to user-specified level
        log.setLevel(log_level)

        # Append any user requested tools to the base set - eliminate
        # duplicate entries
        user_tools = figure_dict.get('tools', '')
        all_tools = FIGURE_TOOLS_BASE + ',' + user_tools
        tool_tokens = all_tools.split(',')
        strip_tool_tokens = []
        for t in tool_tokens:
            strip_tool_tokens.append(t.strip())
        fig_tools = list(set(strip_tool_tokens))

        fig_tool_loc = figure_dict.get('toolbar_location', 'right')

        # Generate the figure instance
        self.fig = figure(tools=fig_tools, toolbar_location=fig_tool_loc)

        # Basic figure attributes
        self.fig.xaxis.axis_label = figure_dict.get('x_label', '')
        self.fig.yaxis.axis_label = figure_dict.get('y_label', '')
        self.fig.title.text = figure_dict.get('title', '')
        self.fig.background_fill_color = figure_dict.get('background_fill_color', 'gainsboro')

        # Append any user requested tooltips to the base set
        user_tips = figure_dict.get('hover_tips', [])
        hover_fig = HoverTool()
        hover_fig.tooltips = HOVER_BASIC_TIPS + user_tips
        self.fig.add_tools(hover_fig)

        self.fig.x_range = figure_dict.get('x_range', None)
        self.fig.y_range = figure_dict.get('y_range', None)
        self.fig.x_range.start = figure_dict.get('xstart', None)
        self.fig.y_range.start = figure_dict.get('ystart', None)
        self.fig.grid.grid_line_color = figure_dict.get('grid_line_color', 'white')

        # These attributes are not for the figure, but are styling attributes
        # used for the "shape" glyphs.  They may be set to non-default values when
        # the build_glyph() routine is invoked.
        self.glyph_color = 'colormap'
        self.color = 'colormap'
        self.size = 10
        self.legend_group = ''
        self.legend_label = ''
        self.fill_alpha = 0.5
        self.line_alpha = 0.5
        self.line_width = 1
        self.fig.legend.location = 'top-right'
        self.fill_color = 'gainsboro'
        self.angle = 0.0
        self.text_size = '10px'

        # This value is set so the functionality is always active - this allows
        # the user click on an item in the legend and turn on/off the display of
        # that particular dataset
        self.fig.legend.click_policy = 'hide'

    def build_glyph(self, glyph_name, x, y, sourceCDS, **data_dict):
        """Generate a glyph representation of X and Y data.

            This method is for generating figures with "shape" glyphs where
            (x,y) coordinates and ColumnDataSource are in use.
6+
            Parameters
            ----------
            glyph_name : str
                Name of the supported "shape" glyph to use (eg., 'circle')

            x, y : str
                Names of X and Y columns of data in sourceCDS.  The columns represent


                the X- and Y-axis coordinates for the center of the circle marker.

            sourceCDS : ColumnDataSource object
                The sourceCDS is an object based upon the Pandas dataframe holding the
                data of interest where the data can be referenced by column name

            name : str, optional
                Indentification of the glyph - useful when using HoverTool where the
                value of the name is a dictionary key as the value itself is variable.

            glyph_color : str, optional
                The glyph_color of the "shape" glyphs can be specified via an input string
                (e.g. 'red', 'blue', 'colormap').  The term `colormap` specifies the
                use of a pre-defined set of colors for the points derived from the
                `colormap` column in the input data `source`.  The default value for
                this value is "colormap".
            legend_label : str, optional
                Legend label to assign the particular glyph ONLY when the "glyph_color"
                argument is not set to "colormap".  If the "glyph_color" is "colormap",
                the "legend_label" argument is ignored.

            These are the characteristics which the developers have found useful thus far.  It
            is expected more arguments may need to be added to support desired functionality.

        """

        # Set the required attributes
        self.glyph_name = glyph_name
        self.x = x
        self.y = y
        self.sourceCDS = sourceCDS

        # Check for optional attributes and use defaults as necessary.  The
        # available attributes are declared in the constructor.
        self.size = data_dict.get('marker_size', self.size)
        self.fill_alpha = data_dict.get('fill_transparency', self.fill_alpha)
        self.fill_color = data_dict.get('fill_color', 'gainsboro')
        self.line_alpha = data_dict.get('line_transparency', self.line_alpha)
        self.line_width = data_dict.get('line_width', self.line_width)
        self.angle = data_dict.get('angle', self.angle)

        # If "glyph_color" is a color, then the specified color and legend label
        # are used.
        # If "glyph_color" is "colormap", then "colormap" and "legend_group" columns
        # from "sourceCDS" are used for the corresponding glyph to set the color
        # and legend text.
        self.glyph_color = data_dict.get('glyph_color', self.glyph_color)
        if self.glyph_color is 'colormap':
            self.legend_group = 'inst_det'
            self.color = 'colormap'
        else:
            self.color = self.glyph_color

        # The "legend_label" is ignored when "glyph_color" is set to "colormap".
        self.legend_label = data_dict.get('legend_label', self.legend_label)

        # Dictionary of supported "shape" glyphs.  These are really references to
        # the associated method names.
        glyph_types = {'circle': HAPFigure.__build_circle_glyph,
                       'square': HAPFigure.__build_square_glyph,
                       'triangle': HAPFigure.__build_triangle_glyph}

        glyph_types[glyph_name](self)

    # "Shape" glyphs
    def __build_circle_glyph(self):

        # hacky if/else
        if self.glyph_color is 'colormap':
            self.fig.circle(x=self.x,
                            y=self.y,
                            source=self.sourceCDS,
                            size=self.size,
                            color=self.color,
                            legend_group=self.legend_group,
                            fill_alpha=self.fill_alpha,
                            line_alpha=self.line_alpha,
                            hover_color='green',
                            name=self.glyph_name)
        else:
            self.fig.circle(x=self.x,
                            y=self.y,
                            source=self.sourceCDS,
                            size=self.size,
                            color=self.color,
                            legend_label=self.legend_label,
                            fill_alpha=self.fill_alpha,
                            line_alpha=self.line_alpha,
                            hover_color='green',
                            name=self.glyph_name)

    def __build_square_glyph(self):

        if self.glyph_color is 'colormap':
            self.fig.square(x=self.x,
                            y=self.y,
                            source=self.sourceCDS,
                            size=self.size,
                            color=self.color,
                            legend_group=self.legend_group,
                            fill_alpha=self.fill_alpha,
                            line_alpha=self.line_alpha,
                            name=self.glyph_name)
        else:
            self.fig.square(x=self.x,
                            y=self.y,
                            source=self.sourceCDS,
                            size=self.size,
                            color=self.color,
                            legend_label=self.legend_label,
                            fill_alpha=self.fill_alpha,
                            line_alpha=self.line_alpha,
                            name=self.glyph_name)

    def __build_triangle_glyph(self):

        if self.glyph_color is 'colormap':
            self.fig.triangle(x=self.x,
                              y=self.y,
                              source=self.sourceCDS,
                              size=self.size,
                              color=self.color,
                              angle=self.angle
                              legend_group=self.legend_group,
                              fill_alpha=self.fill_alpha,
                              line_alpha=self.line_alpha,
                              name=self.glyph_name)
        else:
            self.fig.triangle(x=self.x,
                              y=self.y,
                              source=self.sourceCDS,
                              size=self.size,
                              color=self.color,
                              angle=self.angle
                              legend_label=self.legend_label,
                              fill_alpha=self.fill_alpha,
                              line_alpha=self.line_alpha,
                              name=self.glyph_name)

    def build_vector_glyph(self, x_vbeg, y_vbeg, x_vend, y_vend, x_head, y_head,
                           x_scale_beg, y_scale_beg, x_scale_end, y_scale_end,
                           x_text, y_text, text, **vector_dict):
        """Generate a vector graphic.

            This method is specifically for generating a vector graphic.

            Parameters
            ----------
            x_vbeg, y_vbeg : float
                Coordinates of the starting points of the segments

            x_vend, y_vend : float
                Coordinates of the ending points of the segments

            x_head, y_head : float
                Coordinates for the center of the triangle markers

            x_scale_beg, y_scale_beg : float
                Coordinates of the starting point of the scale/reference segment

            x_scale_end, y_scale_end : float
                Coordinates of the ending point of the scale/reference segment

            x_text, y_text: str
                Coordinates of the text label for the scale/reference segment

            color : str, optional
                The color of the line segments and arrow head
                Default is 'black'

            line_width : float, optional
                The width of the line segments
                Default is 1.0 as set in HAPFigure instantiation

            angle : float, optional
                The color of the line segments and arrow head
                Default is 0.0 as set in HAPFigure instantiation

            marker_size : int, optional
                Size of the triangle glyph used as the arrow head
                Default is 10 as set in HAPFigure instantiation

            text_size : float, optional
                Size of the text label for the scale/reference segment
                Default is '10px' as set in HAPFigure instantiation

        """

        # Get/set the required attributes
        self.x_vbeg = x_vbeg
        self.y_vbeg = y_vbeg
        self.x_vend = x_vend
        self.y_vend = y_vend
        self.x_head = x_head
        self.y_head = y_head
        self.x_scale_beg = x_scale_beg
        self.y_scale_beg = y_scale_beg
        self.x_scale_end = x_scale_end
        self.y_scale_end = y_scale_end
        self.x_text = x_text
        self.y_text = y_text

        # Get/set optional attributes
        color = vector_dict.get('color', 'black')
        line_width = vector_dict.get('line_width', self.line_width)
        angle = vector_dict.get('angle', self.angle)
        size = vector_dict.get('marker_size', self.size)
        text_font_size = vector_dict.get('text_size', self.text_size)
        self.fig.segment(x0=self.x_vbeg,
                         y0=self.y_vbeg,
                         x1=self.x_vend,
                         y1=self.y_vend,
                         color=self.color,
                         line_width=self.line_width)
        self.fig.triangle(x=self.x_head,
                          y=self.y_head,
                          size=self.size,
                          color=self.color,
                          angle=self.angle)

        # Create the reference scale
        # Some attributes are hard-coded as they are never
        # intended to be altered
        self.fig.segment(x0=self.x_scale
                         y0=self.y_scale,
                         x1=self.xr_scale,
                         y1=self.yr_scale,
                         color='black',
                         line_width=self.line_width)
        self.fig.text(x=self.x_text,
                      y=self.y_text,
                      text=self.scale_text,
                      color='black',
                      text_font_size=self.text_size)

    def build_histogram(self, top, bottom, left, right, **data_dict):

        fill_color = data_dict.get('fill_color', 'gray')
        line_color = data_dict.get('line_color', 'black')
        alpha = data_dict.get('fill_transparency', 1.0)

        self.fig.quad(top=top,
                      bottom=bottom,
                      left=left,
                      right=right,
                      fill_color=fill_color,
                      line_color=line_color,
                      alpha=alpha)
