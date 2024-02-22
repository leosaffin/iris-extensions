import numpy as np
import iris
import iris.analysis
import iris.cube
import iris.coords
import iris.util
import numba
from numba import jit


def interpolate(cube, interpolator=iris.analysis.Linear, extrapolation_mode='linear',
         **kwargs):
    """Interpolates to the given coordinates

    Works as a wrapper function to :py:func:`iris.cube.Cube.interpolate`

    Example:

        >>> newcube = interpolate(cube, longitude=0, latitude=45)

    Args:
        cube (iris.cube.Cube):

        interpolator (optional): Instance of the iris interpolator to use.
            Default is :py:class:`iris.analysis.Linear`

        extrapolation_mode (optional): The extrapolation mode for iris to use
            in the case of data outside the co-ordinate bounds. Default is
            linear

        **kwargs: Provides the coordinate value pairs to be interpolated to.

    Returns:
        iris.cube.Cube:
    """
    # Extract the specified output co-ordinates
    points = [(coord, kwargs[coord]) for coord in kwargs]

    # Call the cube's built in interpolation method
    newcube = cube.interpolate(points, interpolator(
                               extrapolation_mode=extrapolation_mode))
    return newcube


def cross_section(cube, xs, xf, ys, yf, npoints,
                  interpolator=iris.analysis.Linear,
                  extrapolation_mode='linear'):
    """ Interpolate a cross section between (xs, ys) and (xf, yf)

    Args:
        cube (iris.cube.Cube): Must be 3D (z,y,x)

        xs (float): Start point of x-coordinate

        xf (float): End point of x-coordinate

        ys (float): Start point of y-coordinate

        yf (float): End point of y-coordinate

        npoints (int, optional): Number of points to interpolate to along the
            cross-section

        interpolator (optional): Instance of the iris interpolator to use.
            Default is :py:class:`iris.analysis.Linear`

        extrapolation_mode (optional): The extrapolation mode for iris to use
            in the case of data outside the co-ordinate bounds. Default is
            linear

    Returns:
        iris.cube.Cube: A 2D cube of the vertical cross-section along the given
        points
    """
    # Create arrays for the set of points along the cross-section
    xpoints = np.linspace(xs, xf, npoints)
    ypoints = np.linspace(ys, yf, npoints)

    # Extract the names of the x and y co-ordinates
    xcoord = cube.coord(axis='x', dim_coords=True).name()
    ycoord = cube.coord(axis='y', dim_coords=True).name()

    # Call the cube's built in interpolation method
    newcube = cube.interpolate(
        [(xcoord, xpoints), (ycoord, ypoints)],
        interpolator(extrapolation_mode=extrapolation_mode)
    )

    # The interpolation returns a box with all corresponding xpoints and
    # ypoints (i.e. a 3d cube). We need to extract the diagonal line along this
    # box to return a 2d cube
    # Demote the y co-ordinate to an auxiliary coord prior to reducing the
    # shape
    iris.util.demote_dim_coord_to_aux_coord(newcube, ycoord)

    # Take the diagonal line along the cube. Use the cubelist functionalityto reduce
    # this to a single cube.
    # TODO change this so we aren't assuming the horizontal coordinates are the last
    # two dimensions
    newcubelist = iris.cube.CubeList()
    for i in range(npoints):
        newcubelist.append(newcube[..., i, i])

    # Reduce to a single cube. This currently has the side effect of always
    # making the y coordinate the 1st dimension and the z coordinate the 2nd
    # dimension. Still plots OK though
    newcube = newcubelist.merge()[0]

    return newcube


def to_level(cube, order=0, **kwargs):
    """ Interpolates to the vertical co-ordinate level

    Args:
        cube (iris.cube.Cube):

        order (int): Order of interpolation. Currently only supports linear (1).

        **kwargs: Provides the coordinate value pair to be interpolated to. Must
            be specified as a list.

            e.g. to interpolate a cube onto a vertical surface of 1000m

            >>> to_level(cube, altitude=[1000])

    Returns:
        iris.cube.Cube: A cube interpolated onto the new vertical co-ordinate.
        Has the same properties as the input cube but with new vertical
        co-ordinates
    """
    if len(kwargs) > 1:
        raise Exception('Can only specify a single vertical co-ordinate')

    # Extract the specified output co-ordinate information
    coord_name = list(kwargs)[0]
    coord_in = cube.coord(coord_name)
    coord_out = kwargs[coord_name]

    # Broadcast coordinate arrays to 3d
    # Input array should match cube shape
    coord_in_points = iris.util.broadcast_to_shape(coord_in.points, cube.shape, coord_in.cube_dims(cube))

    # Output coordinate should be 3d with nz specified by the input
    dims = np.ndim(coord_out)
    if dims == 1:
        ny, nx = cube.shape[1:]
        coord_out_3d = coord_out * np.ones([nx, ny, len(coord_out)])
        coord_out_3d = coord_out_3d.transpose()
    elif dims == 2:
        coord_out_3d = np.expand_dims(coord_out, axis=0)
    elif dims == 3:
        coord_out_3d = coord_out

    else:
        raise Exception('Coordinate must be 2d surface(s) or a list of levels')

    # Select the interpolation flag based on the coordinate
    if 'pressure' in coord_name:
        # Air pressure is interpolated logarithmically
        interp_flag = 1
    else:
        # Otherwise interpolation is linear
        interp_flag = 0

    # Interpolate data
    newdata, mask = down_to_level(cube.data, coord_in_points, coord_out_3d, interp_flag)
    newdata = np.ma.masked_where(mask, newdata)

    # Create a new cube with the new number of vertical levels
    newcube = iris.cube.Cube(
        newdata, long_name=cube.name(), units=cube.units,
        attributes=cube.attributes,
        dim_coords_and_dims=[(cube.coord(axis='y', dim_coords=True), 1),
                             (cube.coord(axis='x', dim_coords=True), 2)])

    if dims == 2:
        newcube = iris.util.squeeze(newcube)

    # Add the new co-ordinate to the output cube
    newcoord = iris.coords.AuxCoord(
        coord_out, long_name=coord_name, units=coord_in.units)
    newcube.add_aux_coord(newcoord, range(newcoord.ndim))

    # Promote single dimensional coordinates to dimensional coordinates
    if dims == 1:
        iris.util.promote_aux_coord_to_dim_coord(newcube, coord_name)

    # If the input was a stack of 2d surfaces then add a dummy vertical coordinate
    # If there is only a single 2d surface then a dummy coordinate is not needed
    elif dims == 3:
        dummy_coord = iris.coords.DimCoord(range(len(coord_out)),
                               long_name='level_number')
        newcube.add_dim_coord(dummy_coord, 0)

    # Add single value coordinates back to the newcube
    add_scalar_coords(cube, newcube)

    return newcube


def add_scalar_coords(cube, newcube):
    for coord in cube.aux_coords:
        if len(coord.points) == 1:
            newcube.add_aux_coord(coord)


def remap_3d(cube, target, vert_coord=None):
    """Remap one cube on to the target mapping

    Args:
        cube (iris.cube.Cube): The cube to be re-mapped

        target (iris.cube.Cube): The cube to re-map to

        vert_coord (str, optional): The name of the coordinate for the vertical
            re-mapping to be done on. Default is None and will use the DimCoord
            for the z-axis

    Returns:
        iris.cube.Cube:
    """
    # Regrid in the horizontal
    cube = cube.regrid(target, iris.analysis.Linear())

    # Interpolate in the vertical
    if vert_coord is None:
        z = target.coord(axis='z', dim_coords=True)
    else:
        z = target.coord(vert_coord)
    cube = cube.interpolate([(z.name(), z.points)], iris.analysis.Linear())

    # Match coordinate information
    newcube = target.copy(data=cube.data)
    newcube.rename(cube.name())
    newcube.units = cube.units

    # Put back correct time information
    for coord in newcube.aux_coords:
        if iris.util.guess_coord_axis(coord) == 'T':
            newcube.remove_coord(coord)

    for coord in cube.aux_coords:
        if iris.util.guess_coord_axis(coord) == 'T':
            newcube.add_aux_coord(coord)

    return newcube


@jit(nopython=True)
def down_to_level(variable, coord_in, coord_out, type):
    # Interpolate to the output point
    if type == 0:
        interp_function = linear
    elif type == 1:
        interp_function = log_linear

    # Loop over zcoord_out points
    npr, ny, nx = coord_out.shape
    output = np.zeros((npr, ny, nx), dtype=float)
    mask = np.zeros((npr, ny, nx), dtype=numba.boolean)
    for l in range(npr):
        for j in range(ny):
            for i in range(nx):
                z = coord_out[l, j, i]
                zcoord = coord_in[:, j, i]
                if z < zcoord.min() or z > zcoord.max():
                    # Set the mask to True if the requested value is not within the
                    # bounds of zcoord
                    mask[l, j, i] = True
                else:
                    # Set the mask to False and do the interpolation if the requested
                    # value is within the bounds of zcoord
                    mask[l, j, i] = False

                    # Find the index of the first zcoord point that has crossed the
                    # output coordinate by searching downwards
                    k = search_downwards(z, zcoord)

                    # Interpolate to the output point
                    output[l, j, i] = interp_function(variable[:, j, i], zcoord, z, k)

    return output, mask

@jit(nopython=True)
def search_downwards(z, zcoord):
    """
    Search downwards to find the index where zcoord first crosses the value z.
    Before using this function it must be checked that z is within the bounds of zcoord
    """
    sign_at_top = np.sign(z - zcoord[-1])
    k = -2
    while sign_at_top == np.sign(z - zcoord[k]):
        k -= 1

    return k


@jit(nopython=True)
def linear(variable, coord, z, k):
    """
    Calculate the interpolation weight assuming linear variation across two points
    """
    alpha = (z - coord[k]) / (coord[k+1] - coord[k])

    return alpha * variable[k+1] + (1 - alpha) * variable[k]


@jit(nopython=True)
def log_linear(variable, coord, z, k):
    """
    Calculate the interpolation weight assuming linear logarithmic variation across two
    points
    """
    alpha = np.log(z / coord[k]) / np.log(coord[k+1] / coord[k])

    return alpha * variable[k+1] + (1 - alpha) * variable[k]
