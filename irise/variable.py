"""This module contains a set of diagnostics to be calculated from a standard
set of prognostic fields.
"""

import numpy as np
import iris.cube
import iris.analysis
from iris.analysis import maths, cartography, MAX

from irise import calculus, constants, grid, interpolate
from irise.fortran import variable as fvariable


def pressure(Pi):
    r"""Calculate the pressure from the Exner function

    :math:`P = P_0 \Pi^{\frac{c_p}{R}}`

    Args:
        Pi (iris.cube.Cube or numpy.ndarray): Exner function of pressure

    Returns:
        Air Pressure (Same type as input)
    """
    return Pi ** (constants.Cp_d / constants.Rd).data * constants.P0


def exner(P):
    r"""Calculate the Exner function of pressure

    :math:`\Pi = (\frac{P}{P_0})^{\frac{R}{c_p}}`

    Args:
        P (iris.cube.Cube or numpy.ndarray): Air pressure

    Returns:
        Exner Pressure (Same type as input)
    """
    return (P / constants.P0) ** (constants.Rd / constants.Cp_d).data


def density(P, T):
    r"""Calculate density using the ideal gas equation

    :math:`\rho = \frac{P}{R T}`

    Args:
        P (iris.cube.Cube or numpy.ndarray): Air pressure

        T (iris.cube.Cube or numpy.ndarray): Temperature

    Returns:
        Air Density (Same type as inputs)
    """
    return P / (T * constants.Rd)


def wind_speed(u, v, w):
    """Calculate the magnitude of the wind speed

    Args:
        u,v,w (iris.cube.Cube or numpy.ndarray): Components of the wind speed

    Returns:
        Wind Speed (Same type as inputs)
    """
    return (u**2 + v**2 + w**2)**0.5


def theta_e(theta, rvs, T):
    """Calculate the equivalent potential temperature

    Args:
        theta (iris.cube.Cube): Potential temperature

        rvs (iris.cube.Cube): Vapour mixing ratio

        T (iris.cube.Cube): Temperature

    Returns:
        theta_e (iris.cube.Cube)
    """
    return theta * maths.exp((rvs * constants.Lv) / (T * constants.Cp_d))


def r_vs(P, e_s):
    """Calculate the relative vapour mixing ratio

    Args:
        P (iris.cube.Cube or numpy.ndarray): Pressure

        e_s (iris.cube.Cube or numpy.ndarray): Saturated vapour mixing ratio

    Returns:
        r_vs (Same type as inputs)
    """
    return (e_s / (P - e_s)) * (constants.water_molecular_weight / constants.dry_air_molecular_weight).data


def vapour_pressure(P, q):
    """Calculate the partial pressure of water vapour

    Args:
        P (iris.cube.Cube or numpy.ndarray): Pressure

        q (iris.cube.Cube or numpy.ndarray): Specific humidity

    Returns:
        Vapour Pressure (Same type as inputs)
    """
    return P * q * (constants.Rv / constants.Rd)


def tetens(T):
    """Calculate the saturated vapour pressure using Teten's approximation

    Args:
        T (iris.cube.Cube): Temperature

    Returns:
        e_s (iris.cube.Cube)
    """
    x = iris.cube.Cube(610.94, units='Pa')
    return maths.exp(17.625 * (T - 273) / (T - 30)) * x


def category_map(*args):
    """Create a category map from multiple binary cubes

    Args:
        *args: Any number of :py:class:`iris.cube.Cube`'s with a 0/1 mapping

    Returns:
        iris.cube.Cube: A single cube with data 0-N, where N is the number of
        cubes input. The data then represents the cube number which is True at
        that location. A `names` attribute, which is a dictionary mapping of
        numbers to cube names, is also added to the cube.
    """
    names = {}

    # Copy the shape and coordinates from one of the cubes
    mapping = args[0].copy()
    mapping.data[:] = 0

    # For each cube set a number for it's mapping
    for n, cube in enumerate(args):
        names[n] = cube.name()
        mapping.data += n * cube.data

    # Save the key mapping numbers to names
    mapping.attributes['names'] = names

    return mapping


def column_integral(cube, **kwargs):
    """Calculate the sum of a cube along the z-axis

    Args:
        cube (iris.cube.Cube): The 3d cube to be collapsed

        **kwargs: Additional keywords to pass to :py:class:`iris.analysis.SUM`

    Returns:
        iris.cube.Cube: A 2d collapsed cube
    """

    z = cube.coord(axis='z', dim_coords=True)

    return cube.collapsed(z.name(), iris.analysis.SUM, **kwargs)


def height(cube):
    """Extract height coordinate as a cube

    Args:
        cube (iris.cube.Cube):

    Returns:
        iris.cube.Cube:
    """
    z = grid.make_cube(cube, 'altitude')

    if z.shape != cube.shape:
        z = grid.broadcast_to_cube(z, cube)

    return z


def surface_height(cube):
    """Extract height coordinate as a cube

    Args:
        cube (iris.cube.Cube)

    Returns:
        iris.cube.Cube:
    """
    z_s = grid.make_cube(cube[0], 'surface_altitude')
    for aux_factory in z_s.aux_factories:
        z_s.remove_aux_factory(aux_factory)

    for coord in z_s.aux_coords:
        z_s.remove_coord(coord)

    return z_s


def geopotential_height(P, levels):
    """Calculate height above sea level on pressure levels

    Args:
        P (iris.cube.Cube):

        levels (list or numpy.ndarray): The values of pressure to output the
            geopotential height on

    Returns:
        iris.cube.Cube: Height on pressure levels
    """
    # Create a height cube
    z = grid.make_cube(P, 'altitude')

    # Add pressure as a coordinate to the height
    pressure_coord = grid.make_coord(P)
    z.add_aux_coord(pressure_coord, range(z.ndim))

    # Interpolate the height on to pressure levels
    geopotential = interpolate.to_level(z, **{P.name(): levels})
    geopotential.rename('Geopotential_height')

    return geopotential


def cloud_top_height(cloud, altitude):
    """The altitude of the highest gridbox containing cloud in each column

    Args:
        cloud (iris.cube.Cube): Cloud amount

        altitude (iris.cube.Cube):

    Returns:
        iris.cube.Cube: cloud top height

    """
    # Mask the altitude so we only see where there is cloud
    z_masked = altitude.copy(data=np.ma.masked_where(cloud.data <= 0, altitude.data))

    # Find the highest point in each column
    z_coord = z_masked.coord(axis="z", dim_coords=True)
    z_cld_top = z_masked.collapsed([z_coord], MAX)

    return z_cld_top


def cloud_thickness(cloud, dz):
    """Total thickness of cloud in the grid column

    Sum of gridbox depths for each gridbox in the column containing cloud

    Args:
        cloud (iris.cube.Cube): Cloud amount

        dz (iris.cube.Cube): Gridbox vertical thickness

    Returns:
        iris.cube.Cube: Total cloud thickness

    """
    cloud.data = (cloud.data > 0).astype(int) * dz.data

    cloud_depth = cloud.collapsed(
        [cloud.coord(axis="z", dim_coords=True)], iris.analysis.SUM
    )

    cloud_depth.units = dz.units

    return cloud_depth


def coriolis_parameter(cube):
    r"""Calculate the Coriolis parameter on the xy grid of the given cube

    :math:`f = 2 \Omega sin(\phi)`

    Args:
        cube (iris.cube.Cube): Any cube with lon/lat coordinates

    Returns:
        iris.cube.Cube: Coriolis parameter as function of lon/lat
    """
    # Calculate the Coriolis parameter
    lat = grid.true_coords(cube)[1]
    f = 2 * constants.omega.data * np.sin(np.deg2rad(lat))

    # Put the output into a cube
    f = iris.cube.Cube(
        f, long_name='coriolis_parameter', units='s-1',
        attributes=cube.attributes,
        dim_coords_and_dims=[(cube.coord(axis='y', dim_coords=True), 0),
                             (cube.coord(axis='x', dim_coords=True), 1)])

    return f


def potential_vorticity(u, v, w, theta, rho):
    r"""Calculate PV

    .. math::

        q = \frac{1}{\rho} (\nabla \times \mathbf{u} + 2 \boldsymbol{\Omega})
            \cdot \nabla \theta

    Args:
        u (iris.cube.Cube): Zonal velocity

        v (iris.cube.Cube): Meridional velocity

        w (iris.cube.Cube): Vertical velocity

        theta (iris.cube.Cube): Potential temperature

        rho (iris.cube.Cube): Density

    Returns:
        iris.cube.Cube: PV in PVU
    """

    # Relative Vorticity
    xterm, yterm, zterm = vorticity(u, v, w)

    # Absolute vorticity
    lat = grid.true_coords(theta)[1]
    f = 2 * constants.omega.data * np.sin(lat * np.pi / 180)
    zterm.data = zterm.data + f

    # Grad(theta)
    dtheta_dx = calculus.polar_horizontal(theta, 'x')
    dtheta_dx = interpolate.remap_3d(dtheta_dx, theta)

    dtheta_dy = calculus.polar_horizontal(theta, 'y')
    dtheta_dy = interpolate.remap_3d(dtheta_dy, theta)

    z = grid.make_cube(theta, 'altitude')
    dtheta_dz = calculus.multidim(theta, z, 'z')
    dtheta_dz = interpolate.remap_3d(dtheta_dz, theta)

    # PV components
    PV_x = xterm * dtheta_dx
    PV_y = yterm * dtheta_dy
    PV_z = zterm * dtheta_dz

    # Full PV
    rho_theta = interpolate.remap_3d(rho, theta)
    epv = (PV_x + PV_y + PV_z) / rho_theta

    epv.rename('ertel_potential_vorticity')
    epv.convert_units('PVU')

    return epv


def isentropic_circulation(pv, pressure, mask=None):
    r"""Calculate the circulation associated with a given PV anomaly

    :math:`C_{\theta} = \iint_R \sigma Q dx dy`

    args:
        pv (iris.cube.Cube or iris.cube.CubeList): Potential vorticity on
            isentropic levels

        pressure (iris.cube.Cube): Pressure at gridpoints

        mask (numpy.ndarray): Mask for area to integrate over

    returns:
        iris.cube.CubeList: The circulation from each PV in q
    """
    # Make sure to iterate over list or cubelist
    if type(pv) == iris.cube.Cube:
        pv = [pv]

    # Calculate pressure halfway between theta levels
    theta_levs = pv[0].coord('air_potential_temperature').points
    dtheta = theta_levs[1] - theta_levs[0]
    theta_levs = list(theta_levs - dtheta / 2)
    theta_levs.append(theta_levs[-1] + dtheta)
    p_theta = interpolate.to_level(
        pressure, air_potential_temperature=theta_levs)

    # Calculate isentropic theta gradient at theta levels
    dp_dtheta = calculus.differentiate(p_theta, 'air_potential_temperature')

    # Calculate isentropic density
    sigma = -1 * dp_dtheta / constants.g

    # Apply the mask to sigma as it is used in each calculation
    if mask is not None:
        sigma.data = np.ma.masked_where(mask, sigma.data)

    # Extract the xy coordinate names to collapse the cube over
    coords = [sigma.coord(axis=axis).name()
              for axis in ['x', 'y']]

    # Calculate the weights to perform area integration
    # weights = iris.analysis.cartography.area_weights(sigma)
    weights = grid.volume(pressure)[-1].data * np.ones_like(pv[0].data)

    circulation = iris.cube.CubeList()
    for pv_i in pv:
        # Calculate the circulation
        c_i = (sigma * pv_i).collapsed(
            coords, iris.analysis.MEAN, weights=weights)

        # Give the cube a meaningful name
        c_i.rename('circulation_due_to_' + pv_i.name())
        circulation.append(c_i)

    return circulation


def mslp(theta, Pi, w, lapse_rate=0.0065, npmsl_height=500.0):
    """Calculate mean sea-level pressure

    Args:
        theta (iris.cube.Cube): Potential temperature at theta-points

        Pi (iris.cube.Cube): Exner pressure at rho-points

        w (iris.cube.Cube): Vertical velocity at theta-points, including
            surface level

        lapse_rate (float): (default is 0.0065)

        npmsl_height (float): (default is 500.0)

    Returns:
        iris.cube.Cube:
    """
    # Extract height coordinate of each grid position
    z_theta = constants.earth_avg_radius.data + w.coord('altitude').points
    z_rho = constants.earth_avg_radius.data + Pi.coord('altitude').points[:-1]
    level_height = theta.coord('level_height').points

    # Calculate grid variables
    grid_lons, grid_lats = cartography.get_xy_grids(theta)
    cos_theta_lat = np.cos(grid_lats)
    delta_lambda = grid_lons[0, 1] - grid_lons[0, 0]
    delta_phi = grid_lats[1, 0] - grid_lats[0, 0]

    # Calculate pressure and exner on theta levels
    Pi = interpolate.interpolate(
        Pi, level_height=w.coord('level_height').points)
    P = pressure(Pi)

    # Extract surface pressure and remove from other variables
    P_star = P[0]
    Pi = Pi[1:]
    P = P[1:]

    # Calculate mean sea-level pressure using hacked MetUM routine
    P_msl = fvariable.calc_pmsl(
        theta.data, Pi.data, P.data, P_star.data, z_theta, z_rho, cos_theta_lat,
        level_height, constants.Cp_d.data, constants.g.data, constants.Rd.data,
        lapse_rate, constants.earth_avg_radius.data, delta_lambda, delta_phi,
        npmsl_height)

    # Copy data into a new cube
    P_msl = P[0].copy(data=P_msl)
    P_msl.rename('air_pressure_at_mean_sea_level')

    return P_msl


def brunt_vaisala(theta, zcoord='altitude'):
    r"""Calculate the Brunt-Vaisala frequency

    :math:`N = \sqrt{\frac{g}{\theta} \frac{d\theta}{dz}}`


    Args:
        theta (iris.cube.Cube): Potential temperature at theta-points

        zcoord (str): The name of the vertical (z) coordinate

    Returns:
        iris.cube.Cube: Brunt-Vaisala frequency at theta points
    """
    N_sq = brunt_vaisala_squared(theta, zcoord=zcoord)
    N = N_sq ** 0.5

    return N


def brunt_vaisala_squared(theta, zcoord='altitude'):
    r"""Calculate the squared Brunt-Vaisala frequency

    :math:`N^2 = \frac{g}{\theta} \frac{d\theta}{dz}`

    Uses a Centred difference to calculate the vertical gradient of potential
    temperature

    Args:
        theta (iris.cube.Cube): Potential temperature

        zcoord (str): The name of the vertical (z) coordinate

    Returns:
        iris.cube.Cube: Squared Brunt-Vaisala frequency at theta points
    """
    # Differentiate theta with respect to altitude
    z = grid.make_cube(theta, zcoord)
    dtheta_dz = calculus.multidim(theta, z, 'z')

    # Interpolate derivatives back to original levels
    zdim = theta.coord(axis='z')
    dtheta_dz = interpolate.interpolate(dtheta_dz, **{zdim.name(): zdim.points})

    # Calculate Brunt-Vaisala frequency squared
    N_sq = (dtheta_dz / theta) * constants.g

    # Correct units because iris is shit sometimes
    N_sq.units = 's-2'

    return N_sq


def Eady(theta, u, P, pressure_levels=(85000, 40000),
         theta_0=iris.cube.Cube(300, units='K')):
    r"""Calculate the Eady growth rate

    :math:`\sigma = 0.31 \frac{f}{N} |\frac{dU}{dz}|`

    Where :math:`f = 2 \Omega sin(\phi)`

    and :math:`N^2 = \frac{g}{\theta_0} \frac{d \theta}{dz}`

    Args:
        theta (iris.cube.Cube): Air potential temperature

        u (iris.cube.Cube): Zonal wind speed

        P (iris.cube.Cube): Air pressure

        pressure_levels:

        theta_0:

    Returns:
        iris.cube.Cube: Eady growth rate
    """
    # Calculate Coriolis parameter
    f = coriolis_parameter(theta)

    # Extract altitude as a cube
    z = grid.make_cube(P, 'altitude')

    # Add pressure as a coordinate to the cubes
    P = grid.make_coord(P)
    for cube in (theta, u, z):
        cube.add_aux_coord(P, [0, 1, 2])

    # Interpolate variables to pressure levels
    theta = interpolate.to_level(theta, air_pressure=pressure_levels)
    u = interpolate.to_level(u, air_pressure=pressure_levels)
    z = interpolate.to_level(z, air_pressure=pressure_levels)

    # Calculate the Brunt-Vaisala frequency
    dtheta_dz = calculus.multidim(theta, z, 'z')
    N_sq = dtheta_dz * (constants.g / theta_0)
    N_sq.units = 's-2'
    N = N_sq ** 0.5

    # Calclate the wind shear
    du_dz = calculus.multidim(u, z, P.name())
    du_dz.data = np.abs(du_dz.data)

    # Calculate the Eady index
    sigma = 0.31 * (du_dz / N) * f

    return sigma


def isentropic_density(p, theta):
    r"""Calculate isentropic density at rho-points

    The isentropic density is defined as
    :math:`\sigma = -\frac{1}{g} \frac{\partial p}{\partial \theta}`

    Since pressure and potential temperature are stored at theta points and g
    is a constant this diagnostic is easily calculated at rho points using
    finite differencing
    """
    # Differentiate pressure with respect to theta
    dp_dtheta = calculus.multidim(p, theta, 'z')

    # Interpolate back to theta levels
    z = theta.coord(axis='z')
    dp_dtheta = interpolate.interpolate(dp_dtheta,  **{z.name(): z.points})

    # Calculate isentropic density
    sigma = -1 * dp_dtheta / constants.g

    return sigma


def vorticity(u, v, w):
    r"""Calculate vorticity at theta-points

    :math:`\xi = \nabla \times \mathbf{u}`

    Args:
        u (iris.cube.Cube): Zonal velocity at u-points on rho-levels

        v (iris.cube.Cube): Meridional velocity at v-points on rho-levels

        w (iris.cube.Cube): Vertical velocity at theta-points

    Returns:
        tuple (iris.cube.Cube, iris.cube.Cube, iris.cube.Cube):
            Three cubes of the different vorticity components
    """
    xi_i, xi_j, xi_k = calculus.curl(u, v, w)

    # Rename components of vorticity
    xi_i.rename('zonal_vorticity')
    xi_j.rename('meridional_vorticity')
    xi_k.rename('vertical_vorticity')

    return xi_i, xi_j, xi_k
