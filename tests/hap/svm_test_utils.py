# Shared helpers for SVM regression tests.

from pathlib import Path

from astropy.io import ascii, fits

from ci_watson.artifactory_helpers import get_bigdata

TESTS_ROOT = Path(__file__).parent
DEFAULT_INPUTS_ROOT = "drizzlepac"
DEFAULT_ENV = "dev"
SVM_INPUT_SUBDIR = "svm_inputs"

def load_poller_filenames(poller_file):
    """Read a poller file and return the list of filenames it references."""
    path = TESTS_ROOT / poller_file
    table = ascii.read(path, format="no_header")
    filename_column = table.colnames[0]
    filenames = list(table[filename_column])
    return filenames


def retrieve_data_for_processing(
    filenames,
    pytestconfig=None,
    artifactory_subdir=SVM_INPUT_SUBDIR,
):
    inputs_root = DEFAULT_INPUTS_ROOT
    env = DEFAULT_ENV

    if pytestconfig:
        roots = pytestconfig.getini("inputs_root")
        first_root = roots if isinstance(roots, str) else next(iter(roots), inputs_root)
        inputs_root = first_root.strip() or inputs_root
        env = (pytestconfig.getoption("env") or env).strip() or env

    artifactory_path = (inputs_root, env, artifactory_subdir)

    expected = {Path(name).name for name in filenames}

    for name in expected:
        path = Path(name)
        if path.exists():
            continue
        get_bigdata(*artifactory_path, name)

    return expected


def build_manifest_name(reference_filename):
    """Derive the manifest filename based on an input exposure."""
    inst = fits.getval(reference_filename, "INSTRUME", ext=0).lower()
    root = fits.getval(reference_filename, "ROOTNAME", ext=0).lower()
    tokens = (inst, root[1:4], root[4:6], "manifest.txt")
    manifest_filename = "_".join(tokens)
    return manifest_filename


def read_manifest(manifest_filename):
    """Read a manifest and return every listed filename."""
    files = []
    with open(manifest_filename, "r", encoding="utf-8") as manifest:
        for line in manifest:
            files.append(line.rstrip("\n"))
    return files


def run_svm_pipeline(poller_file, runner):
    """Execute the SVM pipeline entry point with robust logging."""
    poller_path = TESTS_ROOT / poller_file

    try:
        runner(str(poller_path))
    except Exception as exc:
        raise Exception(f"\nsvm_setup. Exception Visit: {poller_path}\n") from exc