import pytest
import os
import numpy as np
import cdriz_setup


@pytest.fixture
def kernel_pars(request):
    kernel = request.getfixturevalue("kernel")[0]
    # get cdriz.triz inputs from cdriz_setup.py

    # Only offset the output WCS for the "point" kernel,
    # so that we can avoid kernel falling on either side of an edge between
    # pixels.
    offset = 1.0e-5 * np.pi * np.ones(2, dtype=float) if kernel == "point" else None

    params = cdriz_setup.Get_Grid(
        inx=10, iny=10, outx=13, outy=13, offset=offset, background=0.0
    )
    return params

# "square", "point", "turbo", "gaussian", "lanczos3"
@pytest.mark.parametrize(
    "kernel",
    [
        ("square", 10000, 1e-7),
        ("point", 10000, 1e-7),
        ("turbo", 10000, 1e-7),
        ("gaussian", 10000, 1e-3),
        ("lanczos3", 9882.103, 1e-3)
    ]
)
def test_spt_kernels(kernel_pars, kernel, return_png=False):
    """Function tests different c code point kernels (inputs already created on instantiation).

    Parameters
    ----------
    kernel_pars : Class
        The Class inintialized in Get_Grid which includes all of the inputs need to run cdriz.tdriz.
    kernel: tuple [argument passed with parameterize]
        The tuple containing the name of the kernel, the expected sum, and the absolute tolerance.
    return_png:
        Whether to return a png of the outputs.
    """
    kernel_name, expected_sum, atol = kernel

    output_name = f"cdriz_{kernel_name}.png"
    relative_path = "truth_files"
    output_fullpath = cdriz_setup.get_output_fullpath(relative_path, output_name)

    # add bright pixel at center
    kernel_pars.insci[3, 3] = 1e4

    # resample:
    cdriz_setup.cdriz_call(kernel_pars, kernel_name)

    if return_png:
        # save truth file as figure
        cdriz_setup.generate_png(kernel_pars, output_fullpath)

    assert np.allclose(np.sum(kernel_pars.outsci), expected_sum, atol)
