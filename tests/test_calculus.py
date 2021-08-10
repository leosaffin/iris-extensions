import pytest
from iris.tests import stock

from irise import calculus

@pytest.mark.parametrize(
    "cube",
    [
        stock.realistic_4d(),
    ]
)
def test_multidim(cube):
    calculus.multidim(cube, "altitude", "z")


@pytest.mark.parametrize(
    "cube, axes",
    [
        (stock.realistic_3d(), ["t", "y", "x"]),
        (stock.realistic_4d(), ["t", "z", "y", "x"])
    ]
)
def test_diff_by_axis(cube, axes):
    for axis in axes:
        calculus.diff_by_axis(cube, axis)


@pytest.mark.parametrize(
    "cube",
    [
        stock.realistic_4d(),
    ]
)
def test_polar_horizontal(cube):
    calculus.polar_horizontal(cube, "x")
    calculus.polar_horizontal(cube, "y")


@pytest.mark.parametrize(
    "cube",
    [
        stock.realistic_4d(),
    ]
)
def test_div(cube):
    calculus.div(cube, cube, cube)


@pytest.mark.parametrize(
    "cube",
    [
        stock.realistic_4d(),
    ]
)
def test_grad(cube):
    calculus.grad(cube)


@pytest.mark.parametrize(
    "cube",
    [
        stock.realistic_4d(),
    ]
)
def test_curl(cube):
    calculus.curl(cube, cube, cube)


@pytest.mark.parametrize(
    "cube",
    [
        stock.realistic_4d(),
    ]
)
def test_laplacian(cube):
    calculus.laplacian(cube)
