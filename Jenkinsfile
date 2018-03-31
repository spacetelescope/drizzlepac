// Obtain files from source control system.
// if (utils.scm_checkout()) return
node {
    stage("Setup") {
        checkout(scm)
        stash includes: 'drizzlepac', name: 'source_tree', useDefaultExcludes: false
    }
}

// Globals
PIP_INST = "pip install"
CONDA_INST = "conda install -y -q -c http://ssb.stsci.edu/astroconda"
PY_SETUP = "python setup.py"
PYTEST_ARGS = "tests --basetemp=tests_output --junitxml results.xml"

matrix_python = ["2.7", "3.5", "3.6"]
matrix_astropy = ["2", "3"]
matrix = []


// RUN ONCE:
//    "sdist" is agnostic enough to work without any big dependencies
sdist = new BuildConfig()
sdist.nodetype = "linux-stable"
sdist.build_mode = "sdist"
sdist.build_cmds = ["${CONDA_INST} astropy numpy",
                    "${PY_SETUP} sdist"]
matrix += sdist

// Generate installation compatibility matrix
for (python_ver in matrix_python) {
    for (astropy_ver in matrix_astropy) {
        // Astropy >=3.0 no longer supports Python 2.7
        if (python_ver == "2.7" && astropy_ver == "3") {
            continue
        }

        DEPS = "python=${python_ver} astropy=${astropy_ver} \
                fitsblender graphviz nictools numpy numpydoc \
                scipy spherical-geometry sphinx sphinx_rtd_theme \
                stsci_rtd_theme stsci.convolve stsci.image \
                stsci.imagemanip stsci.imagestats stsci.ndimage \
                stsci.skypac stsci.stimage stwcs pyregion setuptools"

        install = new BuildConfig()
        install.nodetype = "linux-stable"
        install.build_mode = "install ${DEPS}"
        install.build_cmds = ["${CONDA_INST} ${DEPS}",
                              "${PY_SETUP} egg_info",
                              "${PY_SETUP} install"]
        matrix += install
    }
}

// Iterate over configurations that define the (distibuted) build matrix.
// Spawn a host of the given nodetype for each combination and run in parallel.
utils.run(matrix)
