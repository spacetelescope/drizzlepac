// Obtain files from source control system.
if (utils.scm_checkout()) return

withCredentials([string(
    credentialsId: 'drizzlepac-codecov',
    variable: 'codecov_token')]) {
// Generate installation compatibility matrix
matrix_python = ["3.6"]
matrix_astropy = [">=3.1,<3.2"]
matrix_numpy = [">=1.16,<1.17"]
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

    MATRIX_SUFFIX = utils.convert_specifiers("${matrix_id}_py${python_ver}_np${numpy_ver}_ap${astropy_ver}")
    MATRIX_TITLE = "mtx-${MATRIX_SUFFIX}"

    bc = new BuildConfig()
    bc.nodetype = "python${python_ver}"
    bc.name = MATRIX_TITLE
    bc.env_vars = ['BUILD_MATRIX_SUFFIX=' + MATRIX_SUFFIX,
                   'BUILD_MATRIX_ID=' + matrix_id,
                   'TEST_BIGDATA=https://bytesalad.stsci.edu/artifactory']
    bc.build_cmds = ["pip install codecov pytest-cov 'numpy${numpy_ver}' 'astropy${astropy_ver}'",
                     "pip install --upgrade -r requirements-dev.txt -e '.[test]'",
                     "pip freeze"]
    bc.test_cmds = ["pytest --cov=./ --basetemp=tests_output --junitxml results.xml --bigdata",
                    "codecov --token=${codecov_token}"]
    bc.test_configs = [data_config]
    matrix += bc
    matrix_id++
}}}


// RUN ONCE:
//    "sdist" is agnostic enough to work without any big, or targeted dependencies
sdist = new BuildConfig()
sdist.nodetype = "python${matrix_python[0]}"
sdist.name = "sdist"
sdist.build_cmds = ["pip install numpy astropy",
                    "python setup.py sdist"]
matrix += sdist


// Iterate over configurations that define the (distibuted) build matrix.
// Spawn a host of the given nodetype for each combination and run in parallel.
utils.run(matrix)
}
