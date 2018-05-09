// Obtain files from source control system.
if (utils.scm_checkout()) return

// Globals
PIP_INST = "pip install"
CONDA_CHANNEL = "http://ssb.stsci.edu/astroconda"
CONDA_CREATE = "conda create -y -q -c ${CONDA_CHANNEL}"
CONDA_INST = "conda install -y -q -c ${CONDA_CHANNEL}"
PY_SETUP = "python setup.py"
PYTEST_ARGS = "tests --basetemp=tests_output --junitxml results.xml"
DEPS = "astropy fitsblender graphviz nictools numpy numpydoc \
        scipy spherical-geometry sphinx sphinx_rtd_theme \
        stsci_rtd_theme stsci.convolve stsci.image \
        stsci.imagemanip stsci.imagestats stsci.ndimage \
        stsci.skypac stsci.stimage stwcs pyregion setuptools"

matrix_python = ["2.7", "3.5", "3.6"]
matrix_astropy = ["2", "3"]
matrix = []


// RUN ONCE:
//    "sdist" is agnostic enough to work without any big dependencies
sdist = new BuildConfig()
sdist.nodetype = "linux-stable"
sdist.build_mode = "sdist"
sdist.build_cmds = ["${CONDA_CREATE} -n dist astropy numpy",
                    "with_env -n dist ${PY_SETUP} sdist"]
matrix += sdist


// RUN ONCE:
//    "build_sphinx" with default python
docs = new BuildConfig()
docs.nodetype = "linux-stable"
docs.build_mode = "docs"
docs.build_cmds = ["${CONDA_CREATE} -n docs ${DEPS}",
                   "with_env -n docs ${PY_SETUP} build build_ext --inplace -- build_sphinx -w"]
matrix += docs


// Generate installation compatibility matrix
for (python_ver in matrix_python) {
    for (astropy_ver in matrix_astropy) {
        // Astropy >=3.0 no longer supports Python 2.7
        if (python_ver == "2.7" && astropy_ver == "3") {
            continue
        }

        DEPS_INST = "python=${python_ver} astropy=${astropy_ver} " + DEPS

        install = new BuildConfig()
        install.nodetype = "linux-stable"
        install.build_mode = "install-${python_ver}"
        install.build_cmds = ["${CONDA_CREATE} -n ${python_ver} ${DEPS_INST}",
                              "with_env -n ${python_ver} ${PY_SETUP} egg_info",
                              "with_env -n ${python_ver} ${PY_SETUP} install"]
        matrix += install
    }
}

// Iterate over configurations that define the (distibuted) build matrix.
// Spawn a host of the given nodetype for each combination and run in parallel.
utils.run(matrix)
