import numpy as np
import pytest
from astropy.io import fits
from astropy.wcs import WCS

from drizzlepac.mapreg import map_region_files
from drizzlepac.util import count_sci_extensions


def _create_mock_fits(base_dir, name="mock_image.fits", n_sci=3):
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


def _create_mock_region(base_dir, name="mock_regions.reg"):
    base_dir.mkdir(parents=True, exist_ok=True)

    region_content = (
        "# Region file format: DS9 version 4.1\n"
        "fk5\n"
        "circle(150d, 2d, 30\")\n"
        "circle(177.4219814,22.2723991,1.405\")\n"
        "box(177.4191036,22.2725963,6.188\",3.119\",302.78378)\n"
        "polygon(177.4241323,22.2771491,177.4222819,22.2767891,177.4241703,22.2753066,177.4247268,22.2761007,177.4246434,22.2768395)\n"

    )

    region_path = base_dir / name
    region_path.write_text(region_content, encoding="utf-8")
    return region_path

# test helpers -----------------------------------------------------------------
def _single_region_input(data_dir):
    return str(_create_mock_region(data_dir))


def _list_region_input(data_dir):
    return [
        str(_create_mock_region(data_dir, name="mock_regions1.reg")),
        str(_create_mock_region(data_dir, name="mock_regions2.reg")),
    ]


def _glob_region_input(data_dir):
    _create_mock_region(data_dir, name="mock_regions1.reg")
    _create_mock_region(data_dir, name="mock_regions2.reg")
    return str(data_dir / "mock_regions*.reg")


# test image inputs ------------------------------------------------------------
@pytest.mark.parametrize("image_specs, expect_list", [pytest.param([("mock_image.fits", 3)], False, id="single-image"),
        pytest.param(
            [("mock_image1.fits", 3), ("mock_image2.fits", 2)],
            True,
            id="multiple-images",
        ),
    ],
)
def test_map_region_files_with_image_variants(tmp_path, image_specs, expect_list):
    data_dir = tmp_path / "inputs"
    input_reg = _create_mock_region(data_dir)
    output_dir = tmp_path / "regions"
    output_dir.mkdir()

    images = [
        _create_mock_fits(data_dir, name=name, n_sci=n_sci)
        for name, n_sci in image_specs
    ]
    images_arg = [str(path) for path in images]
    if not expect_list:
        images_arg = images_arg[0]

    # Disable filtering so partially off-chip regions still produce outputs.
    map_region_files(
        str(input_reg),
        images=images_arg,
        outpath=str(output_dir),
    )

    expected_files = {
        f"{image.stem}_extn{ext}_twreg.reg"
        for image, (_, n_sci) in zip(images, image_specs)
        for ext in range(1, n_sci + 1)
    }

    generated = {path.name for path in output_dir.glob("*.reg")}
    assert generated == expected_files
    for filename in expected_files:
        assert (output_dir / filename).stat().st_size > 0


# test region file inputs -----------------------------------------------------
@pytest.mark.parametrize(
    "region_factory",
    [
        pytest.param(_single_region_input, id="single"),
        pytest.param(_list_region_input, id="list"),
        pytest.param(_glob_region_input, id="glob"),
    ],
)
def test_map_region_files_with_region_variants(tmp_path, region_factory):
    data_dir = tmp_path / "inputs"
    image = _create_mock_fits(data_dir)
    output_dir = tmp_path / "regions"
    output_dir.mkdir()

    region_input = region_factory(data_dir)

    # Disable filtering so partially off-chip regions still produce outputs.
    map_region_files(
        region_input,
        images=str(image),
        outpath=str(output_dir),
    )

    expected_files = {
        f"{image.stem}_extn1_twreg.reg",
        f"{image.stem}_extn2_twreg.reg",
        f"{image.stem}_extn3_twreg.reg",
    }

    generated = {path.name for path in output_dir.glob("*.reg")}
    assert generated == expected_files
    for filename in expected_files:
        assert (output_dir / filename).stat().st_size > 0


# test extension inputs -------------------------------------------------------
@pytest.mark.parametrize(
    "ext_variant",
    [
        pytest.param("sci-string", id="sci-string"),
        pytest.param("subset-list", id="subset-list"),
    ],
)
def test_map_region_files_with_extension_variants(tmp_path, ext_variant):
    data_dir = tmp_path / "inputs"
    input_reg = _create_mock_region(data_dir)
    image = _create_mock_fits(data_dir)
    output_dir = tmp_path / "regions"
    output_dir.mkdir()

    sci_indices = count_sci_extensions(str(image), return_ind=True)
    assert sci_indices, "Expected at least one SCI extension in test image"

    if ext_variant == "sci-string":
        ext_arg = "sci"
        expected_indices = sci_indices
    else:
        ext_arg = sci_indices[:2] if len(sci_indices) > 1 else sci_indices
        expected_indices = ext_arg

    # Disable filtering so partially off-chip regions still produce outputs.
    map_region_files(
        str(input_reg),
        images=str(image),
        img_wcs_ext=ext_arg,
        outpath=str(output_dir),
    )

    generated = {path.name for path in output_dir.glob("*.reg")}
    expected = {
        f"{image.stem}_extn{ext}_twreg.reg" for ext in expected_indices
    }

    assert generated == expected
    for filename in expected:
        assert (output_dir / filename).stat().st_size > 0


