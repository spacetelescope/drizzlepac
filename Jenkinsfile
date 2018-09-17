// @Library('utils@master') _ // For testing unreleased utils library

// Obtain files from source control system.
if (utils.scm_checkout()) return


// Globals
PIP_INST = "pip install --upgrade --upgrade-strategy 'only-if-needed'"
CONDA_CHANNEL = "http://ssb.stsci.edu/astroconda"
CONDA_ARGS = "-y -q -c ${CONDA_CHANNEL}"
CONDA_CREATE = "conda create ${CONDA_ARGS}"
CONDA_INST = "conda install ${CONDA_ARGS}"
PY_SETUP = "python setup.py"
PYTEST = "pytest --basetemp=tests_output --junitxml results.xml --bigdata --remote-data=any"

// The minimum modules required to execute setup.py at all
BASE_DEPS = "astropy numpy"
TEST_DEPS = "pytest pytest-remotedata crds"

// Conda needs explicit dependencies listed
DEPS = "fitsblender graphviz nictools numpydoc \
        pytest pytest-remotedata pyregion \
        scipy spherical-geometry sphinx sphinx_rtd_theme \
        stsci_rtd_theme stsci.convolve stsci.image \
        stsci.imagemanip stsci.imagestats stsci.ndimage \
        stsci.skypac stsci.stimage stwcs setuptools python"

matrix_python = ["3.6", "3.7"]
matrix_astropy = [">=3.0.5"]
matrix_numpy = [">=1.14", "==1.15.0rc1"]
matrix = []


// Configure artifactory ingest
data_config = new DataConfig()
data_config.server_id = 'bytesalad'
data_config.root = 'tests_output'
data_config.match_prefix = '(.*)_result' // .json is appended automatically


// RUN ONCE:
//    "sdist" is agnostic enough to work without any big dependencies
sdist = new BuildConfig()
sdist.nodetype = "linux"
sdist.name = "sdist"
sdist.build_cmds = ["${PIP_INST} ${BASE_DEPS}",
                    "${PY_SETUP} sdist"]
matrix += sdist


// RUN ONCE:
//    "build_sphinx" with default python
/*
docs = new BuildConfig()
docs.nodetype = "linux"
docs.name = "docs"
docs.build_cmds = ["${PIP_INST} ${BASE_DEPS}",
                   "${PY_SETUP} install build_sphinx"]
matrix += docs
*/

// Generate installation compatibility matrix
int matrix_id = 0
for (python_ver in matrix_python) {
for (astropy_ver in matrix_astropy) {
for (numpy_ver in matrix_numpy) {
    MATRIX_SUFFIX = "${matrix_id}_py${python_ver}_np${numpy_ver}_ap${astropy_ver}"
                    .replaceAll("[<>=\\!\\.]", "")
    MATRIX_TITLE = "mtx-${MATRIX_SUFFIX}"
    DEPS_PYTHON= "python=${python_ver} "
    DEPS_EX = ""
    DEPS_EX += "astropy${astropy_ver} "
    DEPS_EX += "numpy${numpy_ver} "

    // "with_env" is a `source activate [env]` wrapper for conda environments
    WRAPPER = "with_env -n ${python_ver}"

    install = new BuildConfig()
    install.nodetype = "linux"
    install.name = MATRIX_TITLE
    install.env_vars = ['BUILD_MATRIX_SUFFIX=${MATRIX_SUFFIX}',
                        'BUILD_MATRIX_ID=${matrix_id}',]
    install.build_cmds = [
        // Install python @ version
        "${CONDA_CREATE} -n ${python_ver} ${DEPS_PYTHON}",

        // Install custom required packages @ version
        "${WRAPPER} ${PIP_INST} ${DEPS_EX}",

        // Install self from source
        "${WRAPPER} ${PY_SETUP} install",
    ]

    install.test_cmds = ["${WRAPPER} pip install ${TEST_DEPS}",
                         "${WRAPPER} ${PYTEST}",]
    install.test_configs = [data_config]
    matrix += install
    matrix_id++
}}}

// Iterate over configurations that define the (distibuted) build matrix.
// Spawn a host of the given nodetype for each combination and run in parallel.
utils.run(matrix)
