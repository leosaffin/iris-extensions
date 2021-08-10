import pytest
from iris.tests import stock

import irise

@pytest.mark.parametrize(
    "cube",
    [
        stock.realistic_4d(),
    ]
)
def test_multidim(cube):
    irise.calculus.multidim(cube, "altitude", "z")


@pytest.mark.parametrize(
    "cube, axes",
    [
        (stock.realistic_3d(), ["t", "y", "x"]),
        (stock.realistic_4d(), ["t", "z", "y", "x"])
    ]
)
def test_diff_by_axis(cube, axes):
    for axis in axes:
        irise.calculus.diff_by_axis(cube, axis)


@pytest.mark.parametrize(
    "cube",
    [
        stock.realistic_4d(),
    ]
)
def test_polar_horizontal(cube):
    irise.calculus.polar_horizontal(cube, "x")
    irise.calculus.polar_horizontal(cube, "y")


@pytest.mark.parametrize(
    "cube",
    [
        stock.realistic_4d(),
    ]
)
def test_div(cube):
    irise.calculus.div(cube, cube, cube)


@pytest.mark.parametrize(
    "cube",
    [
        stock.realistic_4d(),
    ]
)
def test_grad(cube):
    irise.calculus.grad(cube)


@pytest.mark.parametrize(
    "cube",
    [
        stock.realistic_4d(),
    ]
)
def test_curl(cube):
    irise.calculus.curl(cube, cube, cube)


@pytest.mark.parametrize(
    "cube",
    [
        stock.realistic_4d(),
    ]
)
def test_laplacian(cube):
    irise.calculus.laplacian(cube)
