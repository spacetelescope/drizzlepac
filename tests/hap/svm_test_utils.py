"""Shared helpers for SVM regression tests.

These utilities consolidate repeated logic across the ``test_svm_*`` modules
including poller-file parsing, data retrieval, manifest handling, and running
of pipeline steps. Keeping them in one place reduces duplication and makes it
simpler to update test mechanics in the future.
"""
from __future__ import annotations

import datetime
import os
from pathlib import Path
from typing import Iterable, Sequence, Set

import pytest
from astropy.io import ascii, fits

from drizzlepac.haputils import astroquery_utils as aqutils

TESTS_ROOT = Path(__file__).parent


def load_poller_filenames(poller_file: str) -> list[str]:
    """Read a poller file and return the list of filenames it references."""
    path = TESTS_ROOT / poller_file
    table = ascii.read(path, format="no_header")
    filename_column = table.colnames[0]
    filenames = list(table[filename_column])
    return filenames


def change_to_temp_working_dir(tmp_path_factory, test_filename: str) -> Path:
    """Create a per-test temporary directory and change into it."""
    workdir = tmp_path_factory.mktemp(Path(test_filename).name)
    os.chdir(workdir)
    return workdir


def _collect_wildcards(filenames: Iterable[str], suffixes: Sequence[str]) -> dict[str, str]:
    wildcards: dict[str, str] = {}
    for name in filenames:
        lower_name = name.lower()
        for suffix in suffixes:
            suffix_lower = suffix.lower()
            if lower_name.endswith(f"{suffix_lower}.fits") and suffix not in wildcards:
                wildcards[suffix] = name[:6] + "*"
    return wildcards


def retrieve_data_for_processing(
    filenames: Iterable[str],
    suffixes: Sequence[str] = ("FLC", "FLT"),
    product_type: str = "pipeline",
) -> Set[str]:
    """Download required observation files and cull unused extras."""
    suffixes = tuple(s.upper() for s in suffixes)
    wildcards = _collect_wildcards(filenames, suffixes)

    downloaded: list[str] = []
    for suffix in suffixes:
        flag = wildcards.get(suffix)
        if flag:
            files = aqutils.retrieve_observation(flag, suffix=[suffix], product_type=product_type)
            downloaded.extend(files)

    filenames_set = set(filenames)
    downloaded_set = set(downloaded)

    if not downloaded_set:
        existing_files = {fn for fn in filenames if Path(fn).exists()}
        if existing_files:
            return existing_files

    files_to_process = filenames_set & downloaded_set if downloaded_set else filenames_set

    files_to_remove = filenames_set.symmetric_difference(downloaded_set)
    for ftr in files_to_remove:
        try:
            os.remove(ftr)
        except FileNotFoundError:
            continue
        except Exception as exc:
            raise Exception(f"The file {ftr} could not be deleted from disk. ") from exc
    return files_to_process


def build_manifest_name(reference_filename: str) -> str:
    """Derive the manifest filename based on an input exposure."""
    inst = fits.getval(reference_filename, "INSTRUME", ext=0).lower()
    root = fits.getval(reference_filename, "ROOTNAME", ext=0).lower()
    tokens = (inst, root[1:4], root[4:6], "manifest.txt")
    manifest_filename = "_".join(tokens)
    return manifest_filename


def read_manifest(manifest_filename: str) -> list[str]:
    """Read a manifest and return every listed filename."""
    files: list[str] = []
    with open(manifest_filename, "r", encoding="utf-8") as manifest:
        for line in manifest:
            files.append(line.rstrip("\n"))
    return files


def run_svm_pipeline(poller_file: str, runner) -> None:
    """Execute the SVM pipeline entry point with robust logging."""
    poller_path = TESTS_ROOT / poller_file

    current_dt = datetime.datetime.now()

    try:
        runner(str(poller_path))
    except Exception as exc:  # pragma: no cover - defensive logging
        raise Exception(f"\nsvm_setup. Exception Visit: {poller_path}\n") from exc

    final_dt = datetime.datetime.now()