from pathlib import Path

import pytest

from drizzlepac.mapreg import map_region_files
from drizzlepac.util import count_sci_extensions

PROJECT_ROOT = Path(__file__).resolve().parents[1]


def _project_path(filename: str) -> Path:
    return PROJECT_ROOT / filename


def test_map_region_files_with_sci_string(tmp_path):
    input_reg = _project_path("mapreg_test2.reg")
    image = _project_path("jcdua3f4q_flc.fits")
    output_dir = tmp_path / "regions"
    output_dir.mkdir()

    # Disable filtering so partially off-chip regions still produce outputs.
    map_region_files(
        str(input_reg),
        images=str(image),
        img_wcs_ext="sci",
        outpath=str(output_dir),
        filter=None,
    )

    sci_indices = count_sci_extensions(str(image), return_ind=True)
    assert sci_indices, "Expected at least one SCI extension in test image"

    generated = {path.name for path in output_dir.glob("*.reg")}
    expected = {
        f"{image.stem}_extn{ext}_twreg.reg" for ext in sci_indices
    }

    assert generated == expected
    for filename in expected:
        assert (output_dir / filename).stat().st_size > 0


def test_map_region_files_with_explicit_extensions(tmp_path):
    input_reg = _project_path("mapreg_test2.reg")
    image = _project_path("jcdua3f4q_flc.fits")
    output_dir = tmp_path / "regions"
    output_dir.mkdir()

    sci_indices = count_sci_extensions(str(image), return_ind=True)
    assert sci_indices

    selected = sci_indices[:2] if len(sci_indices) > 1 else sci_indices

    # Disable filtering so partially off-chip regions still produce outputs.
    map_region_files(
        str(input_reg),
        images=str(image),
        img_wcs_ext=selected,
        outpath=str(output_dir),
        filter=None,
    )

    generated = {path.name for path in output_dir.glob("*.reg")}
    expected = {
        f"{image.stem}_extn{ext}_twreg.reg" for ext in selected
    }

    assert generated == expected
    for filename in expected:
        assert (output_dir / filename).stat().st_size > 0


def test_map_region_files_rejects_compound_string(tmp_path):
    input_reg = _project_path("mapreg_test2.reg")
    image = _project_path("jcdua3f4q_flc.fits")

    with pytest.raises(ValueError):
        map_region_files(
            str(input_reg),
            images=str(image),
            img_wcs_ext="SCI,1",
            outpath=str(tmp_path),
        )


def test_map_region_files_rejects_negative_extension(tmp_path):
    input_reg = _project_path("mapreg_test2.reg")
    image = _project_path("jcdua3f4q_flc.fits")

    with pytest.raises(ValueError):
        map_region_files(
            str(input_reg),
            images=str(image),
            img_wcs_ext=[-1],
            outpath=str(tmp_path),
        )
