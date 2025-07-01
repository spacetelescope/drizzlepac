# Drizzlepac

[![build](https://github.com/spacetelescope/drizzlepac/actions/workflows/build.yml/badge.svg?branch=main)](https://github.com/spacetelescope/drizzlepac/actions/workflows/build.yml)
[![codecov](https://codecov.io/gh/spacetelescope/drizzlepac/branch/master/graph/badge.svg)](https://codecov.io/gh/spacetelescope/drizzlepac)
[![docs](https://readthedocs.org/projects/drizzlepac/badge/?version=latest)](http://drizzlepac.readthedocs.io/en/latest/?badge=latest)
[![Powered by STScI Badge](https://img.shields.io/badge/powered%20by-STScI-blue.svg?colorA=707170&colorB=3e8ddd&style=flat)](http://www.stsci.edu)
[![Powered by Astropy Badge](http://img.shields.io/badge/powered%20by-AstroPy-orange.svg?style=flat)](http://www.astropy.org/)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3743274.svg)](https://doi.org/10.5281/zenodo.3743274)

Nightly regression test results are run using [github actions](https://github.com/spacetelescope/RegressionTests/actions/workflows/drizzlepac.yml).

The use of this software on HST data is described at:

    http://drizzlepac.stsci.edu/

A complete description of the documented interfaces in the code itself 
can be found at:

    http://drizzlepac.readthedocs.io


# Installation

## Conda (Recommended)

`Drizzlepac` is installed when you install the `stenv` conda environment (a replacement for `astroconda`). Select your desired release and follow the instructions on the [`stenv` installation page](https://stenv.readthedocs.io/en/latest/getting_started.html). 

## Install with pip

```bash
$ pip install git+https://github.com/spacetelescope/drizzlepac.git
```
The option `--no-use-pep517` MAY be required in order to correctly build 
the C extensions with `pip` versions up to 22.2, after commenting out 
the `build-backend` from the `pyproject.toml` config file.

Support for installing using `pip` is still evolving, so use of this 
command is provided on an experimental basis for now.

## From Source

### Clone this repository
```bash
$ git clone https://github.com/spacetelescope/drizzlepac
$ cd drizzlepac
```

### Build the documentation

*Note:* If you intend to use `drizzlepac`'s embedded help feature from within
an interactive `python` or `ipython` session, we recommend you do not skip
this step.

```bash
$ cd doc/
$ make html
```

### Install DrizzlePac

```bash
$ python setup.py install
```
