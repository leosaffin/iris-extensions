import pytest
import numpy as np
from iris.tests import stock

from irise import interpolate


def test_interpolate():
    pass


@pytest.mark.parametrize(
    "cube, result",
    [
        (stock.realistic_4d(), 0),
    ]
)
def test_cross_section(cube, result):
    result = interpolate.cross_section(cube[0], 359.59, 359.60, -0.12, -0.11, 100)

    assert result.ndim == 2
    assert result.dim_coords[0].name() == "grid_latitude"
    assert result.dim_coords[1].name() == "model_level_number"
    np.testing.assert_almost_equal(result.data[0, 0], 287.97137, decimal=5)


def test_to_level():
    cube = stock.realistic_4d()
    result = interpolate.to_level(cube[0], altitude=[5000])

    assert result.ndim == 3
    assert result.dim_coords[0].name() == "altitude"
    assert result.dim_coords[1].name() == "grid_latitude"
    assert result.dim_coords[2].name() == "grid_longitude"
    np.testing.assert_almost_equal(result.data[0, 0, 0], 312.4549, decimal=4)


