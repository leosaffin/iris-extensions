import pytest
from iris.tests import stock

from irise import grid


@pytest.mark.parametrize("cube", [stock.realistic_4d()])
def test_polar_coords(cube):
    theta, phi, rho = grid.polar_coords(cube)


@pytest.mark.parametrize("cube", [stock.realistic_3d(), stock.realistic_4d()])
def test_true_coords(cube):
    lon, lat = grid.true_coords(cube)


@pytest.mark.parametrize("cube", [stock.realistic_4d()])
def test_volume(cube):
    volume = grid.volume(cube)


@pytest.mark.parametrize("cube", [stock.realistic_4d()])
def test_make_coord(cube):
    grid.make_coord(cube)


@pytest.mark.parametrize("cube", [stock.realistic_4d()])
def test_make_cube(cube):
    grid.make_cube(cube, "altitude")


def test_add_cube_as_coord():
    pass


def test_add_hybrid_height():
    pass


def add_hybrid_pressure():
    pass


def test_extract_dim_coord():
    pass


@pytest.mark.parametrize("cube", [stock.realistic_3d(), stock.realistic_4d()])
def test_get_datetime(cube):
    grid.get_datetime(cube)
