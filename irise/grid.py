import numpy as np
import iris.aux_factory
import iris.coords
import iris.util
from iris.analysis import cartography

from irise import constants
from irise.fortran import grid as fgrid



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


def volume(cube):
    """Calculate the volume of grid boxes

    Args:
        cube (iris.cube.Cube): An iris cube with the relevant coordinate
            information. x/y grids and the altitude coordinate are required.

    Returns:
        iris.cube.Cube: A cube copied from the original with the volume
        information
    """
    # Get the coordinates
    bounds = constants.earth_avg_radius.data + cube.coord('altitude').bounds
    theta, phi, rho = polar_coords(cube)

    # Calculate the volume in each gridbox as a spherical integral
    data = fgrid.volume(rho, bounds, theta[0, :], phi[:, 0])

    # Convert the output to a cube
    output = cube.copy(data=data)
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
    newcube = cube.copy(data=coord.points)
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


def extract_dim_coord(cube, axis):
    """Extracts a coordinate from a cube based on the axis it represents

    Uses :py:func:`iris.util.guess_coord_axis` to match a dimensional coordinate
    with the requested axis

    Args:
        cube (iris.cube.Cube):

        axis (str): The co-ordinate axis to take from the cube. Must be one of
            X,Y,Z,T.

    Returns:
        iris.coords.DimCoord:
            The dimensional coordinate matching the requested axis on the given
            cube.

    Raises:
        ValueError: If axis is not one of X,Y,Z,T

        KeyError: If the cube does not contain a coord with the requested axis
    """
    axis = axis.upper()
    # If the axis supplied is not correct raise an error
    if axis not in ['X', 'Y', 'Z', 'T']:
        raise ValueError('Must specify a co-ordinate axis')

    # Loop over dimensional coords in the cube
    for coord in cube.dim_coords:
        # Return the coordinate if it matches the axis
        if axis == iris.util.guess_coord_axis(coord):
            return coord

    # If no coordinate matches raise an error
    raise KeyError('Cube has no coordinate for axis ' + axis)


def get_datetime(cube):
    """Extract the time coordinate from a cube in datetime format
    """
    tcoord = cube.coord('time')
    tcoord_as_datetime = tcoord.units.num2date(tcoord.points)

    return tcoord_as_datetime
