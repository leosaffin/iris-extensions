import pytest
import numpy as np
import iris

from irise import variable



@pytest.mark.parametrize(
    "function, inputs, result",
    [
        (variable.pressure, [0.5], 8838.834764831841),
        (variable.exner, [85000], 0.9546275831395891),
        (variable.density, [1e5, 300], 1.161248029720029),
        (variable.wind_speed, [10, 5, 0.1], 11.180787092150535),
        (variable.theta_e, [300, 0.01, 290], 326.8880672467619),
        (variable.r_vs, [1e5, 1e3], 0.006282393030885892),
        (variable.vapour_pressure, [1e5, 0.01], 1607.828426421411),
        (variable.tetens, [290], 1934.1100957308577)
    ]
)
def test_scalars(function, inputs, result):
    inputs = [iris.cube.Cube(data=value) for value in inputs]

    x = function(*inputs)

    np.testing.assert_almost_equal(x.data, result)


# category_map
# column_integral
# height
# surface_height
# geopotential_height
# coriolis_parameter
# potential_vorticity
# isentropic_circulation
# mslp
# brunt_vaisala
# brunt_vaisala_squared
# Eady
# isentropic_density
# vorticity

