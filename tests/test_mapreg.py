from pathlib import Path

import numpy as np
import pytest
from astropy.io import fits
from astropy.wcs import WCS

from drizzlepac.mapreg import map_region_files
from drizzlepac.util import count_sci_extensions


def _create_mock_fits(base_dir: Path, name: str = "mock_image.fits", n_sci: int = 3) -> Path:
    base_dir.mkdir(parents=True, exist_ok=True)

    wcs = WCS(naxis=2)
    wcs.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    wcs.wcs.crval = [150.0, 2.0]
    wcs.wcs.crpix = [50.0, 50.0]
    wcs.wcs.cd = np.array([
        [-2.7777777e-4, 0.0],
        [0.0, 2.7777777e-4],
    ])

    data = np.ones((100, 100), dtype=np.float32)
    hdus = [fits.PrimaryHDU()]

    for extver in range(1, n_sci + 1):
        header = wcs.to_header()
        header["EXTNAME"] = "SCI"
        header["EXTVER"] = extver
        hdus.append(fits.ImageHDU(data=data.copy(), header=header))

    image_path = base_dir / name
    fits.HDUList(hdus).writeto(image_path, overwrite=True)
    return image_path


def _create_mock_region(base_dir: Path, name: str = "mock_regions.reg") -> Path:
    base_dir.mkdir(parents=True, exist_ok=True)

    region_content = (
        "# Region file format: DS9 version 4.1\n"
        "fk5\n"
        "circle(150d, 2d, 30\")\n"
        "circle(150.0003d, 2.0003d, 20\")\n"
    )

    region_path = base_dir / name
    region_path.write_text(region_content, encoding="utf-8")
    return region_path


def test_map_region_files_with_sci_string(tmp_path):
    data_dir = tmp_path / "inputs"
    input_reg = _create_mock_region(data_dir)
    image = _create_mock_fits(data_dir)
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
    data_dir = tmp_path / "inputs"
    input_reg = _create_mock_region(data_dir)
    image = _create_mock_fits(data_dir)
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
    data_dir = tmp_path / "inputs"
    input_reg = _create_mock_region(data_dir)
    image = _create_mock_fits(data_dir)

    with pytest.raises(ValueError):
        map_region_files(
            str(input_reg),
            images=str(image),
            img_wcs_ext="SCI,1",
            outpath=str(tmp_path),
        )
