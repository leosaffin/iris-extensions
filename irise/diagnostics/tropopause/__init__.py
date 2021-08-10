import numpy as np
import iris.cube

from irise import grid
from irise.fortran import diagnostic as fdiagnostic


def dynamical(pv, q, pvtrop=2.0, qmax=0.001):
    """Locates the tropopause based on the PV definition

    Finds the tropopause using a downwards search for a PV values and returns
    the height followed by subsequent downward searches for the height of any
    tropopause folds. Specific humidity is used to distinguish tropopause folds
    from diabatic PV anomalies.

    Args:
        pv (iris.cube.Cube): Potential Vorticity

        q (iris.cube.Cube): Specific humidity

        pvtrop (float, optional): The value of PV at the tropopause. Default is
            2.0

        qmax (float, optional): The maximum value of specific humidity at the
            tropopause. Default is 0.001

    Returns:
        tuple(iris.cube.Cube, iris.cube.Cube, iris.cube.Cube):
            The height of the dynamical tropopause and any folds (top/bottom).
    """
    z = pv.coord('altitude').points

    # Find the tropopause height
    trop, fold_top, fold_bottom = fdiagnostic.dynamical_tropopause(
        np.abs(pv.data), q.data, z, pvtrop, qmax)

    # Put the result in a cube
    xcoord = pv.coord(axis='X', dim_coords=True)
    ycoord = pv.coord(axis='Y', dim_coords=True)

    trop = iris.cube.Cube(
        trop, long_name='dynamical_tropopause_altitude', units='m',
        dim_coords_and_dims=[(ycoord, 0), (xcoord, 1)])

    fold_top = iris.cube.Cube(
        fold_top, long_name='dynamical_tropopause_fold_top_altitude',
        units='m', dim_coords_and_dims=[(ycoord, 0), (xcoord, 1)])

    fold_bottom = iris.cube.Cube(
        fold_bottom, long_name='dynamical_tropopause_fold_bottom_altitude',
        units='m', dim_coords_and_dims=[(ycoord, 0), (xcoord, 1)])

    return trop, fold_top, fold_bottom


def thermal(temperature, zmin=4500, threshold=-0.002, dz=2000):
    r"""Calculate the height of the thermal tropopause

    The thermal tropopause is defined as the lowest level where the lapse rate

    :math:`\frac{\partial T}{\partial z} < 0.002 K m^{-1}`,

    and the average lapse rate

    :math:`\bar{\frac{\partial T}{\partial z}} < 0.002 K m^{-1}`,

    for 2 km above.

    Args:
        temperature (iris.cube.Cube): Temperature in Kelvin. Must include a
            vertical coordinate

        zmin (float, optional): Height to start upwards search for tropopause.
            This avoids diagnosing a boundary layer inversion as the
            tropopause. Default is 2000m.

        threshold (float, optional): The temperature lapse rate threshold to
            define the thermal tropopause. Default is 0.002.

        dz (float, optional): Distance above the tropopause which must have an
            average stratospheric lapse rate. Default is 2000m

    Returns:
        iris.cube.Cube: A 2d map of the thermal tropopause height
    """
    z = temperature.coord('altitude').points
    # Find the tropopause height
    ztrop = fdiagnostic.thermal_tropopause(
        temperature.data, z, zmin, threshold, dz)

    # Put the result in a cube
    try:
        xcoord = temperature.coord(axis='X', dim_coords=True)
        ycoord = temperature.coord(axis='Y', dim_coords=True)
        dim_coords_and_dims = [(ycoord, 0), (xcoord, 1)]
    except NameError:
        dim_coords_and_dims = []

    ztrop = iris.cube.Cube(
        ztrop, long_name='lapse_rate_tropopause_altitude', units='m',
        dim_coords_and_dims=dim_coords_and_dims)

    return ztrop
