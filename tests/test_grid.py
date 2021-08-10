import pytest
from iris.tests import stock

import irise


@pytest.mark.parametrize("cube", [stock.realistic_4d()])
def test_polar_coords(cube):
    theta, phi, rho = irise.grid.polar_coords(cube)


@pytest.mark.parametrize("cube", [stock.realistic_3d(), stock.realistic_4d()])
def test_true_coords(cube):
    lon, lat = irise.grid.true_coords(cube)


@pytest.mark.parametrize("cube", [stock.realistic_4d()])
def test_volume(cube):
    volume = irise.grid.volume(cube)


@pytest.mark.parametrize("cube", [stock.realistic_4d()])
def test_make_coord(cube):
    irise.grid.make_coord(cube)


@pytest.mark.parametrize("cube", [stock.realistic_4d()])
def test_make_cube(cube):
    irise.grid.make_cube(cube, "altitude")


def test_add_cube_as_coord():
    pass


def test_add_hybrid_height():
    pass


def test_add_hybrid_pressure():
    pass


@pytest.mark.parametrize("cube", [stock.realistic_3d(), stock.realistic_4d()])
def test_get_datetime(cube):
    irise.grid.get_datetime(cube)
