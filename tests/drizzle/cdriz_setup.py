import os
import numpy as np
from astropy import wcs
from drizzlepac import cdriz


def get_wcs(_grid, pscale=0.04):
    w = wcs.WCS()
    w.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    w.wcs.crpix = [_grid[0] / 2.0, _grid[1] / 2.0]  # define as center pixel in array
    w.wcs.crval = [10.000, 10.000]  # RA/Decl. in decimal degrees
    w.wcs.cdelt = [pscale / 3600.0, pscale / 3600.0]
    w.wcs.set()
    return w


class Get_Grid:
    """Sets up inputs to call the python wrapper of the c code: cdriz.tdriz.
    The size of the input and ouput grids can be specified in the arguments for
     the init(). For an input grid of 4 x 4 and output grid of 5 x 5 you would
     use the following call of the class.

        Get_Grid(inx=4,iny=4, outx=5, outy=5)

        mapping = cdriz.DefaultWCSMapping(
            input_wcs, output_wcs,
            input_wcs.pixel_shape[0], input_wcs.pixel_shape[1],
            stepsize
        )

    """

    def __init__(self, inx, iny, outx, outy):
        np.random.seed(0)  # keep same random across each instance
        self.in_grid = (iny, inx)  # use numpy indexing order (y,x)
        self.out_grid = (outy, outx)  # use numpy indexing order (y,x)
        self.insci = np.random.randn(self.in_grid[0], self.in_grid[1]).astype("float32")
        self.inwht = np.ones(self.in_grid, dtype=np.float32)
        self.dny = self.insci.shape[1]  # input y_grid
        self.outsci = np.zeros(self.out_grid, dtype=np.float32)
        self.outwht = np.zeros(self.out_grid, dtype=np.float32)
        self.outctx = np.zeros(self.out_grid, dtype=np.int32)
        self.w1 = get_wcs(self.in_grid)
        self.w2 = get_wcs(self.out_grid)
        self.mapping = cdriz.DefaultWCSMapping(
            self.w1, self.w2, self.in_grid[0], self.in_grid[1], 1
        )

    def zero_background(self):
        super().__init__()
        self.insci = np.zeros(self.in_grid, dtype=np.float32)


def cdriz_call(_set_kernel_pars, kernel):
    """
    parameters explained in c code (arrdrizmodule.c)

    _vers, nmiss, nskip = cdriz.tdriz(insci, inwht, outsci, outwht,
        outctx, uniqid, ystart, 1, 1, _dny,
        pix_ratio, 1.0, 1.0, 'center', pixfrac,
        kernel, in_units, expscale, wt_scl,
        fillval, nmiss, nskip, 1, mapping)

    """
    cdriz.tdriz(
        _set_kernel_pars.insci,
        _set_kernel_pars.inwht,
        _set_kernel_pars.outsci,
        _set_kernel_pars.outwht,
        _set_kernel_pars.outctx,
        1,  # uniqid
        0,  # ystart
        1,  # xmin pixel to start reading; fits so starts at 1
        1,  # ymin pixel to start reading; fits so starts at 1
        _set_kernel_pars.dny,  # _dny
        1.0,  # pix_ratio
        1.0,  # xscale; plate scale variations
        1.0,  # yscale; plate scale variations
        "corner",  # center of pixel is what 0 corresponds to
        1.0,  # pixfrac
        kernel,  # kernel
        "cps",  # units
        1.0,  # expscale
        1.0,  # wt_scale
        "INDEF",  # fillval
        0,  # nmiss
        0,  # nskip
        1,  # vflag; historical value only (not used)
        _set_kernel_pars.mapping,
    )


def get_output_fullpath(relative_path, output_filename):
    """Returns the final output path given a name and relative path. 

    Parameters
    ----------
    relative_path : str
        relative path from this file to desired path. 
    output_filename : str
        desired name of output file (png or csv). 

    Returns
    -------
    output_fullpath: str
        final output path
    """    
    local_path = os.path.dirname(__file__)
    output_name = os.path.join(relative_path, output_filename)
    output_fullpath = os.path.join(local_path, output_name)
    print(f"Output Path: {output_fullpath}")
    return output_fullpath


def save_array(_data, _name):
    np.savetxt(
        _name,
        X=_data,
        fmt="%1.8f",
        delimiter=",",
    )


def generate_png(_set_kernel_pars, _name):
    """Added function for displaying the test results as an image map.
    This should show the expected point kernel response.

    Parameters
    ----------
    _set_kernel_pars : Class instance
        an instance of the kernel_pars class.
    _name : str
        The name of the output png file.
    """
    # for generating truth files
    import matplotlib.pyplot as plt

    fig = plt.figure(figsize=(4, 2))
    ax1 = fig.add_subplot(111, projection=_set_kernel_pars.w1)
    ax1.imshow(_set_kernel_pars.outsci, origin="lower", cmap="Greys")
    ax1.set_ylabel(" ")
    ax1.set_xlabel(" ")
    fig.savefig(_name)


def error_message(_data, _name):
    """Saves new truth csv file on failure of test.

    Parameters
    ----------
    _data : np.array
        data to save to new truth file
    _name : str
        new name of truth file, should be slightly
        different than current truth file
    """
    save_array(_data, _name)
