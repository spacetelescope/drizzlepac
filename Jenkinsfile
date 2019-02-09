// Obtain files from source control system.
if (utils.scm_checkout()) return

// Generate installation compatibility matrix
matrix_python = ["3.6", "3.7"]
matrix_astropy = [">=3.1.0"]
matrix_numpy = [">=1.14", "==1.15.0"]
matrix = []

// Configure artifactory ingest
data_config = new DataConfig()
data_config.server_id = 'bytesalad'
data_config.root = 'tests_output'
data_config.match_prefix = '(.*)_result' // .json is appended automatically

int matrix_id = 0
for (python_ver in matrix_python) {
for (astropy_ver in matrix_astropy) {
for (numpy_ver in matrix_numpy) {

    MATRIX_SUFFIX = "${matrix_id}_py${python_ver}_np${numpy_ver}_ap${astropy_ver}"
                    .replaceAll("[<>=\\!\\.]", "")
    MATRIX_TITLE = "mtx-${MATRIX_SUFFIX}"

    bc = new BuildConfig()
    bc.nodetype = "linux"
    bc.name = MATRIX_TITLE
    bc.env_vars = ['BUILD_MATRIX_SUFFIX=' + MATRIX_SUFFIX,
                        'BUILD_MATRIX_ID=' + matrix_id]
    bc.conda_channels = ['http://ssb.stsci.edu/astroconda']
    bc.conda_packages = ['fitsblender',
                         'graphviz',
                         'nictools',
                         'numpydoc',
                         'crds',
                         'scipy',
                         'spherical-geometry',
                         'sphinx',
                         'sphinx_rtd_theme',
                         'stsci_rtd_theme',
                         'stsci.image',
                         'stsci.imagemanip',
                         'stsci.imagestats',
                         'stsci.skypac',
                         'stregion',
                         'stsci.stimage',
                         'setuptools',
                         // test dependencies
                         'pytest=3.8.2',
                         'pytest-remotedata',
                         'crds']
    // Dependencies that change based on build matrix entry.
    bc.conda_packages += ["astropy${astropy_ver}",
                          "numpy${numpy_ver}",
                          "python=${python_ver}"]
    bc.build_cmds = ["python setup.py install"]
    bc.test_cmds = ["pytest --basetemp=tests_output --junitxml results.xml --bigdata --remote-data=any"]
    bc.test_configs = [data_config]
    matrix += bc
    matrix_id++
}}}


// RUN ONCE:
//    "sdist" is agnostic enough to work without any big dependencies
sdist = new BuildConfig()
sdist.nodetype = "linux"
sdist.name = "sdist"
sdist.conda_packages = ['astropy',
                        'numpy']
sdist.build_cmds = ["python setup.py sdist"]
matrix += sdist


//    "build_sphinx" with default python
//docs = new BuildConfig()
//docs.nodetype = "linux"
//docs.name = "docs"
//sdist.conda_packages = ['astropy',
//                        'numpy']
//docs.build_cmds = ["python setup.py install build_sphinx"]
//matrix += docs


// Iterate over configurations that define the (distibuted) build matrix.
// Spawn a host of the given nodetype for each combination and run in parallel.
utils.run(matrix)
