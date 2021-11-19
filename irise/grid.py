import numpy as np
import iris.aux_factory
import iris.coords
import iris.util
from iris.analysis import cartography

from irise import constants


def polar_coords(cube):
    """Get spherical polar coordinates from a cube

    Args:
        cube (iris.cube.Cube): An iris cube with the relevant coordinate
            information. XY grids and the altitude coordinate are required.

    Returns:
        tuple (numpy.ndarray, numpy.ndarray, numpy.ndarray):
            The zonal and meridional co-ordinates in radians and the vertical
            co-ordinate in metres
    """
    rho = constants.earth_avg_radius.data + cube.coord('altitude').points
    x, y = cartography.get_xy_grids(cube)

    theta = x * (np.pi / 180)
    phi = (90 - y) * (np.pi / 180)

    return theta, phi, rho


def true_coords(cube):
    """Extract latitude and longitude from a cube with rotated coordinates

    Args:
        cube (iris.cube.Cube):

    Returns:
        tuple (numpy.ndarray, numpy.ndarray):
            N-dimensional arrays of longitude and latitude
    """
    lon, lat = cartography.get_xy_grids(cube)
    cs = cube.coord_system()
    if cs.grid_mapping_name == 'rotated_latitude_longitude':
        lon, lat = cartography.unrotate_pole(
            lon, lat, cs.grid_north_pole_longitude, cs.grid_north_pole_latitude)

    return lon, lat


def area(cube):
    # Calling area_weights on a cube with a rotated pole will cause an error if the cube
    # longitude/latitude and grid_longitude/grid_latitude coordinates. Just remove the
    # longitude/latitude coordinates temporarily in this case
    try:
        a = iris.analysis.cartography.area_weights(cube)
    except ValueError:
        cube = cube.copy()
        cube.remove_coord("longitude")
        cube.remove_coord("latitude")
        a = iris.analysis.cartography.area_weights(cube)

    # The area weights calculation uses the Earth radius as part of the calculation
    # Rescale by a factor to account for the altitude above earth
    #z = ...
    #a = (z**2 / iris.analysis.cartography.DEFAULT_SPHERICAL_EARTH_RADIUS**2) * area

    return a


def thickness(cube, z_name="altitude"):
    """Calculate the vertical thickness of gridboxes

    Args:
        cube (iris.cube.Cube):

        z_name (str):

    Returns:

    """
    z_bounds = cube.coord(z_name).bounds
    dz = z_bounds[..., 1] - z_bounds[..., 0]
    dz = iris.util.broadcast_to_shape(dz, cube.shape, cube.coord_dims(z_name))

    # Convert the output to a cube
    output = cube.copy(data=dz)
    output.units = "m"

    return output


def volume(cube, z_name="altitude"):
    """Calculate the volume of grid boxes

    Args:
        cube (iris.cube.Cube): An iris cube with the relevant coordinate
            information. x/y grids and the altitude coordinate are required.

        z_name (str):

    Returns:
        iris.cube.Cube: A cube copied from the original with the volume
        information
    """
    vol = area(cube) * thickness(cube, z_name=z_name)

    # Convert the output to a cube
    output = cube.copy(data=vol)
    output.rename('volume')
    output.units = 'm^3'

    return output


def make_coord(cube, bounds=None):
    """Makes an auxiliary coordinate from the cube

    Args:
        cube (iris.cube.Cube):

        bounds (iris.cube.Cube): Optional bounds for the output coordinate. Must
            have a size of N+1, where N is the size of the primary cube used to
            make the coordinate.

    Returns:
        iris.coords.AuxCoord:
    """
    if bounds is not None:
        # The bounds attribute of an iris coord has shape (N, 2)
        bounds = np.array([bounds.data[:-1], bounds.data[1:]]).transpose()

    coord = iris.coords.AuxCoord(
        cube.data, long_name=cube.name(), units=cube.units, bounds=bounds)

    return coord


def make_cube(cube, coord_name):
    """Make a cube from a coordinate

    Args:
        cube (iris.cube.Cube): Cube with coordinate attributes

        coord_name (str): The name of the coordinate to extract from the cube

    Returns:
        iris.cube.Cube:
    """
    coord = cube.coord(coord_name)

    # Subselect the cube by the dimensions of the requested coordinate
    axes = cube.coord_dims(coord_name)

    sl = [slice(None)] * cube.ndim

    for n in range(cube.ndim):
        if n not in axes:
            sl[n] = slice(0, 1)

    newcube = iris.util.squeeze(cube[tuple(sl)]).copy(data=coord.points)
    newcube.units = coord.units
    newcube.rename(coord.name())

    return newcube


def add_cube_as_coord(cube, coord_cube):
    """Add a coordinate to the cube along the common dimensions

    Args:
        cube (iris.cube.Cube): The cube to be modified

        coord_cube (iris.cube.Cube): The cube to be added as a coordinate. All
            dimensions in coord_cube must be in cube
    """
    coord = make_coord(coord_cube)

    common_dimensions = None

    cube.add_aux_coord(coord, common_dimensions)

    return


def add_hybrid_height(cube):
    """Add the hybrid height factory to the given cube
    """
    delta = cube.coord('level_height')
    sigma = cube.coord('sigma')
    orography = cube.coord('surface_altitude')

    factory = iris.aux_factory.HybridHeightFactory(
        delta=delta, sigma=sigma, orography=orography)

    cube.add_aux_factory(factory)

    return


def add_hybrid_pressure(cube):
    """Add the hybrid pressure factory to the given cube
    """
    sigma = cube.coord('sigma')
    surface_pressure = cube.coord('surface_pressure')

    factory = iris.aux_factory.HybridPressureFactory(
        sigma=sigma, surface_air_pressure=surface_pressure)

    cube.add_aux_factory(factory)

    return


def get_datetime(cube):
    """Extract the time coordinate from a cube in datetime format
    """
    tcoord = cube.coord('time')
    tcoord_as_datetime = tcoord.units.num2date(tcoord.points)

    return tcoord_as_datetime


def broadcast_to_cube(cube, target):
    """
    Extend `iris.util.broadcast_to_shape` to broadcast a cube to a larger shape

    Args:
        cube:
        target:

    Returns:

    """
    c = target.dim_coords

    common_dimensions = [c.index(coord) for coord in cube.dim_coords]

    broadcast_data = iris.util.broadcast_to_shape(
        cube.data, target.shape, common_dimensions,
    )

    broadcast_cube = iris.cube.Cube(
        data=broadcast_data,
        standard_name=cube.standard_name,
        long_name=cube.long_name,
        var_name=cube.var_name,
        units=cube.units,
        attributes=cube.attributes,
        cell_methods=target.cell_methods,
        dim_coords_and_dims=[
            (coord, target.coord_dims(coord)) for coord in target.dim_coords
        ],
        aux_coords_and_dims=[
            (coord, target.coord_dims(coord)) for coord in target.aux_coords
        ],
        aux_factories=target.aux_factories,
        cell_measures_and_dims=[
            (cell_measure, target.cell_measure_dims(cell_measure))
            for cell_measure in target.cell_measures()
        ],
        ancillary_variables_and_dims=[
            (ancil_var, target.ancillary_variable_dims(ancil_var))
            for ancil_var in target.ancillary_variables()
        ],
    )

    return broadcast_cube
