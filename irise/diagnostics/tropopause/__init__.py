import numpy as np
import iris.cube
from numba import jit

from irise.interpolate import linear, search_downwards


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
    trop, fold_top, fold_bottom = dynamical_tropopause(
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
    ztrop = thermal_tropopause(temperature.data, z, zmin, threshold, dz)

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


@jit(nopython=True)
def dynamical_tropopause(pv, q, z, pvtrop, qtrop):
    """Find the height of the dynamical tropopause for each grid point using a PV based definition
    """
    nz, ny, nx = pv.shape

    ztrop = np.zeros((ny, nx))
    fold_top = np.zeros((ny, nx))
    fold_bottom = np.zeros((ny, nx))

    # Loop over grid columns
    for j in range(ny):
        for i in range(nx):
            k = nz - 1
            ztrop[j, i] = dynamical_tropopause_sc(pv[:,j,i], q[:,j,i], z[:,j,i], pvtrop, qtrop, k)
            if k > 1:
                dynamical_tropopause_sc(pv[:,j,i], q[:,j,i], z[:,j,i], pvtrop, qtrop, fold_top(j,i), k)
            if k > 1:
                dynamical_tropopause_sc(pv[:,j,i], q[:,j,i], z[:,j,i], pvtrop, qtrop, fold_bottom(j,i), k)

    return  ztrop, fold_top, fold_bottom


@jit(nopython=True)
def dynamical_tropopause_sc(pv, q, z, pvtrop, qtrop, k):
    """Find the height of the dynamical tropopause"""
    found = True

    # Search downward for where PV crosses threshold
    while k >= 0 and found:
        k = k-1
        # Check whether the pv crosses the threshold
        if (pv[k+1] > pvtrop) != (pv[k] > pvtrop):
            # Check the humidity threshold
            if q[k] < qtrop:
                # Interpolate to the exact height where the threshold is crossed
                return linear(z, pv, pvtrop, k)

    # Return zero if nothing found
    return 0.0


@jit(nopython=True)
def thermal_tropopause(T, z, zmin, lapse_rate_trop, dz):
    """
    Find the height of the thermal tropopause for each grid point using the WMO
    temperature lapse rate definition
    """
    nz, ny, nx = T.shape
    ztrop = np.zeros((ny, nx), dtype=float)

    # Loop over grid columns
    for j in range(ny):
        for i in range(nx):
            k = 1
            while z[k, j, i] < zmin:
                k = k+1

            ztrop[j, i] = thermal_tropopause_sc(T[k:,j,i], z[k:,j,i], lapse_rate_trop, dz)

    return ztrop


@jit(nopython=True)
def thermal_tropopause_sc(T, z, lapse_rate_trop, dz):
    """
    Find the height of the thermal tropopause using the WMO temperature lapse rate
    definition
    """
    nz = len(T)
    z_lr = np.zeros(2, dtype=float)
    lapse_rate = np.zeros(2, dtype=float)

    # Search upward for where lapse rate drops below threshold
    k = 0
    while k < nz - 1:
        k += 1

        # Store the previous lapse rate and altitude
        lapse_rate[0] = lapse_rate[1]
        z_lr[0] = z_lr[1]

        # Calculate the lapse rate between two grid boxes
        z_lr[1] = 0.5 * (z[k+1] + z[k])
        lapse_rate[1] = (T[k+1] - T[k]) / (z[k+1] - z[k])

        # Check whether the lapse crosses the threshold
        if lapse_rate[0] < lapse_rate_trop < lapse_rate[1]:
            # Interpolate to the exact height where the threshold is crossed
            ztrop = linear(z_lr, lapse_rate, lapse_rate_trop, 1)

            # Interpolate to temperature at the tropopause
            T_trop = linear(T, z, ztrop, k, nz)

            # Check that the average lapse rate is stratospheric above the tropopause
            # (next 2km for WMO)
            z_above = ztrop + dz

            # If the top of the domain is within 2km use the average lapse rate to the
            # top of the domain
            if z[nz] <= z_above:
                z_above = z[nz]
                T_above = T[nz]
            else:
                k_above = search_downwards(z_above, z)
                # Interpolate to temperature above the tropopause
                T_above = linear(T, z, z_above, k_above)

            # Calculate the average lapse rate over the height
            av_lapse_rate = (T_above - T_trop) / (z_above - ztrop)

            if av_lapse_rate > lapse_rate_trop:
                # If tropopause is found then end while loop
                return ztrop

    # Return zero if nothing found
    return 0.0
