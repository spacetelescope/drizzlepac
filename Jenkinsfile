// Obtain files from source control system.
if (utils.scm_checkout()) return

def test_bin(env_name, bin) {
    def result = "with_env -n ${env_name} ${bin}"
    return result
}

def test_import(env_name, module) {
    def result = "with_env -n ${env_name} python -c 'import ${module}'"
    return result
}

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
matrix_numpy = ["latest"]
matrix = []


// RUN ONCE:
//    "sdist" is agnostic enough to work without any big dependencies
sdist = new BuildConfig()
sdist.nodetype = "linux"
sdist.name = "sdist"
sdist.build_cmds = ["${CONDA_CREATE} -n dist astropy numpy",
                    "with_env -n dist ${PY_SETUP} sdist"]
matrix += sdist


// RUN ONCE:
//    "build_sphinx" with default python
docs = new BuildConfig()
docs.nodetype = "linux"
docs.name = "docs"
docs.build_cmds = ["${CONDA_CREATE} -n docs ${DEPS}",
                   "with_env -n docs ${PY_SETUP} build_sphinx"]
matrix += docs


// Generate installation compatibility matrix
for (python_ver in matrix_python) {
    for (astropy_ver in matrix_astropy) {
        for (numpy_ver in matrix_numpy) {
            // Astropy >=3.0 no longer supports Python 2.7
            if (python_ver == "2.7" && astropy_ver == "3") {
                continue
            }

            DEPS_INST = "python=${python_ver} "

            if (astropy_ver != "latest") {
                DEPS_INST += "astropy=${astropy_ver} "
            }

            if (numpy_ver != "latest") {
                DEPS_INST += "numpy=${numpy_ver} "
            }

            DEPS_INST += DEPS

            install = new BuildConfig()
            install.nodetype = "linux"
            install.name = "install-py=${python_ver},np=${numpy_ver},ap=${astropy_ver}"
            install.build_cmds = ["${CONDA_CREATE} -n ${python_ver} ${DEPS_INST}",
                                  "with_env -n ${python_ver} ${PY_SETUP} egg_info",
                                  "with_env -n ${python_ver} ${PY_SETUP} install",

                                  test_bin(python_ver, 'mdriz -h'),
                                  test_bin(python_ver, 'resetbits -h'),
                                  test_bin(python_ver, 'updatenpol -h'),
                                  test_bin(python_ver, 'runastrodriz -h'),

                                  test_import(python_ver, 'drizzlepac'),
                                  test_import(python_ver, 'drizzlepac.ablot'),
                                  test_import(python_ver, 'drizzlepac.acsData'),
                                  test_import(python_ver, 'drizzlepac.adrizzle'),
                                  test_import(python_ver, 'drizzlepac.astrodrizzle'),
                                  test_import(python_ver, 'drizzlepac.buildmask'),
                                  test_import(python_ver, 'drizzlepac.buildwcs'),
                                  test_import(python_ver, 'drizzlepac.catalogs'),
                                  test_import(python_ver, 'drizzlepac.cdriz'),
                                  test_import(python_ver, 'drizzlepac.createMedian'),
                                  test_import(python_ver, 'drizzlepac.drizCR'),
                                  test_import(python_ver, 'drizzlepac.findobj'),
                                  test_import(python_ver, 'drizzlepac.imageObject'),
                                  test_import(python_ver, 'drizzlepac.imagefindpars'),
                                  test_import(python_ver, 'drizzlepac.imgclasses'),
                                  test_import(python_ver, 'drizzlepac.irData'),
                                  test_import(python_ver, 'drizzlepac.linearfit'),
                                  test_import(python_ver, 'drizzlepac.mapreg'),
                                  test_import(python_ver, 'drizzlepac.mdriz'),
                                  test_import(python_ver, 'drizzlepac.mdzhandler'),
                                  test_import(python_ver, 'drizzlepac.minmed'),
                                  test_import(python_ver, 'drizzlepac.nicmosData'),
                                  test_import(python_ver, 'drizzlepac.outputimage'),
                                  test_import(python_ver, 'drizzlepac.photeq'),
                                  test_import(python_ver, 'drizzlepac.pixreplace'),
                                  test_import(python_ver, 'drizzlepac.pixtopix'),
                                  test_import(python_ver, 'drizzlepac.pixtosky'),
                                  test_import(python_ver, 'drizzlepac.processInput'),
                                  test_import(python_ver, 'drizzlepac.quickDeriv'),
                                  test_import(python_ver, 'drizzlepac.refimagefindpars'),
                                  test_import(python_ver, 'drizzlepac.regfilter'),
                                  test_import(python_ver, 'drizzlepac.resetbits'),
                                  test_import(python_ver, 'drizzlepac.runastrodriz'),
                                  test_import(python_ver, 'drizzlepac.sky'),
                                  test_import(python_ver, 'drizzlepac.skytopix'),
                                  test_import(python_ver, 'drizzlepac.staticMask'),
                                  test_import(python_ver, 'drizzlepac.stisData'),
                                  test_import(python_ver, 'drizzlepac.tweakback'),
                                  test_import(python_ver, 'drizzlepac.tweakreg'),
                                  test_import(python_ver, 'drizzlepac.updatehdr'),
                                  test_import(python_ver, 'drizzlepac.updatenpol'),
                                  test_import(python_ver, 'drizzlepac.util'),
                                  test_import(python_ver, 'drizzlepac.wcs_functions'),
                                  test_import(python_ver, 'drizzlepac.wfc3Data'),
                                  test_import(python_ver, 'drizzlepac.wfpc2Data'),]
            matrix += install
        }
    }
}

// Iterate over configurations that define the (distibuted) build matrix.
// Spawn a host of the given nodetype for each combination and run in parallel.
utils.run(matrix)
