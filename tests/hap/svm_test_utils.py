# Shared helpers for SVM regression tests.

import datetime
import os
from pathlib import Path

import pytest
from astropy.io import ascii, fits

from ci_watson.artifactory_helpers import BigdataError, get_bigdata

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


def change_to_temp_working_dir(tmp_path_factory, test_filename):
    """Create a per-test temporary directory and change into it."""
    workdir = tmp_path_factory.mktemp(Path(test_filename).name)
    os.chdir(workdir)
    return workdir


def retrieve_data_for_processing(
    filenames,
    suffixes=("FLC", "FLT"),
    product_type="pipeline",
    *,
    pytestconfig=None,
    artifactory_subdir=SVM_INPUT_SUBDIR,
):
    """Materialize required observation files locally via Artifactory."""

    _ = suffixes, product_type  # Parameters retained for backward compatibility.

    inputs_root = DEFAULT_INPUTS_ROOT
    env = DEFAULT_ENV

    if pytestconfig is not None:
        try:
            configured_roots_raw = pytestconfig.getini("inputs_root")
        except (TypeError, ValueError):
            configured_roots_raw = []

        if isinstance(configured_roots_raw, str):
            configured_roots = [configured_roots_raw]
        else:
            configured_roots = list(configured_roots_raw)

        if configured_roots:
            candidate_root = configured_roots[0].strip()
            if candidate_root:
                inputs_root = candidate_root

        env_candidate = (pytestconfig.getoption("env") or "").strip()
        if env_candidate:
            env = env_candidate

    artifactory_path = (inputs_root, env, artifactory_subdir)

    expected = {Path(name).name for name in filenames}
    staged = set()
    missing = []

    for name in expected:
        if Path(name).exists():
            staged.add(name)
            continue

        try:
            local_path = get_bigdata(*artifactory_path, name)
        except BigdataError:
            missing.append(name)
            continue

        staged.add(Path(local_path).name)

    if missing:
        missing_list = ", ".join(sorted(missing))
        raise BigdataError(f"Missing SVM input data for: {missing_list}")

    files_to_remove = staged - expected
    for orphan in files_to_remove:
        try:
            os.remove(orphan)
        except FileNotFoundError:
            continue
        except Exception as exc:
            raise Exception(f"The file {orphan} could not be deleted from disk.") from exc

    return staged & expected


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

    current_dt = datetime.datetime.now()

    try:
        runner(str(poller_path))
    except Exception as exc:  # pragma: no cover - defensive logging
        raise Exception(f"\nsvm_setup. Exception Visit: {poller_path}\n") from exc

    final_dt = datetime.datetime.now()