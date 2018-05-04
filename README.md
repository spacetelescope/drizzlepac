# Drizzlepac

[![Build Status](https://ssbjenkins.stsci.edu/job/STScI/job/drizzlepac/job/master/badge/icon)](https://ssbjenkins.stsci.edu/job/STScI/job/drizzlepac/job/master/) [![Documentation Status](https://readthedocs.org/projects/drizzlepac/badge/?version=latest)](http://drizzlepac.readthedocs.io/en/latest/?badge=latest)


The use of this software on HST data is described at:

    http://drizzlepac.stsci.edu/

A complete description of the documented interfaces in the code itself 
can be found at:

    http://drizzlepac.readthedocs.io


# Installation

## Conda (Recommended)

```bash
$ conda config --add channels http://ssb.stsci.edu/astroconda
$ conda create -n astroconda stsci
```

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
$ python setup.py build_sphinx
```

### Install drizzlepac

```bash
$ python setup.py install
```
