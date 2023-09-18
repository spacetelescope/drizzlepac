Installation
------------

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

##### SUPPORT for PIP Installation:
Installation tools are evolving to rely on commands such as `pip` 
to build and install software.  This package can now be installed 
using the following command:

```bash
$ pip install .
```
The option `--no-use-pep517` MAY be required in order to correctly build 
the C extensions with `pip` versions up to 22.2, after commenting out 
the `build-backend` from the `pyproject.toml` config file.

**Support for installing using `pip` is still evolving, so use of this 
command is provided on an experimental basis for now.**
