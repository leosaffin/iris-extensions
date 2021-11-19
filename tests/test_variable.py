import pytest
import numpy as np
import iris

import irise



@pytest.mark.parametrize(
    "function, inputs, result",
    [
        (irise.variable.pressure, [0.5], 8838.834764831841),
        (irise.variable.exner, [85000], 0.9546275831395891),
        (irise.variable.density, [1e5, 300], 1.161248029720029),
        (irise.variable.wind_speed, [10, 5, 0.1], 11.180787092150535),
        (irise.variable.theta_e, [300, 0.01, 290], 326.8880672467619),
        (irise.variable.r_vs, [1e5, 1e3], 0.006282393030885892),
        (irise.variable.vapour_pressure, [1e5, 0.01], 1607.828426421411),
        (irise.variable.tetens, [290], 1934.1100957308577)
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
# cloud_top_height
# cloud_thickness
# coriolis_parameter
# potential_vorticity
# isentropic_circulation
# mslp
# brunt_vaisala
# brunt_vaisala_squared
# Eady
# isentropic_density
# vorticity

