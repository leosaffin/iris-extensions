import iris.cube
import pytest

from irise import convert

def test_recursion():
    cubes = iris.cube.CubeList()
    pytest.raises(ValueError, convert.calc, "air_pressure", cubes)