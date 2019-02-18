import numpy as np
import iris
from iris.coords import DimCoord, AuxCoord
from iris.analysis import MEAN

from irise import interpolate, grid


def averaged_over(x, bins, bin_variable, weights, mask=None):
    """Calculates the average of x in each of the bins

    Args:
        x (iris.cube.Cube or iris.cube.CubeList): The variable(s) to average
            between the bins

        bins: Sequence of scalars defining the bin edges

        bin_variable (iris.cube.Cube): The variable that decides where to
            average the input variable. Must be the same shape.

        weights (iris.cube.Cube or None): The weights for calculating a weighted
            average (e.g. mass). If weights= :py:class:`None`, an even weighting
            is assumed.

        mask (np.ndarray of bool, optional): An array to subset the domain
            averaged over. Elements set as :py:class:`True` will be ignored.

    Returns:
        cubes (iris.cube.CubeList): A cubelist containing a cube of the
        weighted mean of each variable across the bins and a cube of the
        integrated weight in each bin
    """
    # Check that inputs are correct
    if type(x) == iris.cube.Cube:
        x = [x]
    elif type(x) == iris.cube.CubeList or type(x) == list:
        pass
    else:
        raise TypeError('Must specify cube or cubelist')

    if weights is None:
        weights = bin_variable.copy()
        weights.data[...] = 1.0

    # Mask input arrays
    bin_array = np.ma.masked_where(mask, bin_variable.data).compressed()
    weight_array = np.ma.masked_where(mask, weights.data).compressed()
    x_array = [np.ma.masked_where(mask, cube.data).compressed() for cube in x]

    # Calculate the integral of the weights in each bin
    weight_hist = np.histogram(bin_array, bins=bins, weights=weight_array)[0]

    # Calculate the weighted integral of each variable in the bins
    x_hist = [np.histogram(bin_array, bins=bins, weights=y * weight_array)[0] /
              weight_hist for y in x_array]

    # Put the output in to a cubelist
    cubes = iris.cube.CubeList()
    # Create a coordinate based on the bins
    bounds = np.zeros([len(bins)-1, 2])
    bounds[:, 0] = bins[:-1]
    bounds[:, 1] = bins[1:]
    points = [(x1 + x2) / 2 for x1, x2 in bounds]
    coord = DimCoord(points, long_name=bin_variable.name(),
                     units=bin_variable.units, bounds=bounds)

    # Make a cube for each input variable
    for n, data in enumerate(x_hist):
        cube = iris.cube.Cube(data, long_name=x[n].name(),
                              units=x[n].units,
                              dim_coords_and_dims=[(coord, 0)])
        cubes.append(cube)

    # Make a cube for the weights
    cube = iris.cube.Cube(weight_hist, long_name=weights.name(),
                          units=weights.units,
                          dim_coords_and_dims=[(coord, 0)])
    cubes.append(cube)

    return cubes


def depth_average(cube, coord, minimum, maximum, axis, aggregator=MEAN,
                  **kwargs):
    """Calculate the average of the cube across a set range of a coordinate

    Example:
        I want to calculate the average low-level PV from my model data. I have
        PV and pressure as 3d cubes output from my model. If I define low
        levels as being between 85000 Pa and 70000 Pa I would calculate the
        average PV in that depth as

        >>> low_level_pv = depth_average(PV, pressure, 70000, 85000, 'z')

    Args:
        cube (iris.cube.Cube): The cube to average

        coord (str, iris.coords.Coord or iris.cube.Cube): The coordinate or name
            of coordinate on the cube to perform the averaging of the cube
            along

        minimum (scalar): The minimum of the specified coordinate

        maximum (scalar): The maximum of the specified coordinate

        axis (str): The axis to collapse the cube along (x,y,z,t)

        aggregator (iris.analysis.Aggregator): The aggregator used to collapse
            the cube on each level. Default is :py:class:`iris.analysis.MEAN`.

        **kwargs: Additional keywords to be supplied to the
            :py:func:`iris.cube.Cube.collapsed` method (e.g. weights).

    Returns:
        iris.cube.Cube: The integral of the input cube along one dimension.
        The result will have a shape reduced by one dimension compared to the
        input cube.
    """
    # Don't mess with the input cube
    cube = cube.copy()

    # Mask the cube outside the coordinate bounds
    if type(coord) == str:
        coord = cube.coord(coord)
    elif type(coord) == iris.cube.Cube:
        coord = grid.make_coord(coord)

    mask = np.logical_or(coord.points < minimum, coord.points > maximum)
    cube.data = np.ma.masked_where(mask, cube.data)

    # Calculate the mean of the cube across the coordinate
    dim_coord = grid.extract_dim_coord(cube, axis)
    integral = cube.collapsed(dim_coord, aggregator, **kwargs)

    return integral


def profile(x, surface, dz, mask=None, aggregator=MEAN, **kwargs):
    """Calculate an averaged vertical profile relative to a specified surface

    Args:
        x (iris.cube.Cube or iris.cube.CubeList): The variable(s) to produce
            profiles of.

        surface (iris.cube.Cube): 2d cube specifying the altitude of the
            surface to interpolate relative to.

        dz: Iterable container of floats prescribing the distances from the
            given surface to interpolate to.

        mask (array, optional): A binary array with the same shape as the cubes
            data. Default is None.

        aggregator (iris.analysis.Aggregator): The aggregator used to collapse
            the cube on each level. Default is MEAN.

        **kwargs: Additional keywords to be supplied to the
            :py:func:`iris.cube.Cube.collapsed` method (e.g. weights).

    Returns:
        iris.cube.CubeList: A vertical profile for each variable in x
    """
    # Check that inputs are correct
    if type(x) == iris.cube.Cube:
        x = [x.copy()]
    elif type(x) == iris.cube.CubeList or type(x) == list:
        x = [cube.copy() for cube in x]
    else:
        raise TypeError('Must specify cube or cubelist')

    # Loop over each cube
    cubes = iris.cube.CubeList()
    for cube in x:
        # Create a new coordinate for the cube
        newcoord = cube.coord('altitude').points - surface.data
        coord_name = 'distance_from_' + surface.name()
        newcoord = AuxCoord(newcoord, long_name=coord_name, units='m')
        cube.add_aux_coord(newcoord, [0, 1, 2])

        # Interpolate the cube relative to the surface
        newcube = interpolate.to_level(cube, **{coord_name: dz})

        # Apply the mask to the cube
        if mask is not None:
            newcube.data = np.ma.masked_where(
                mask * np.ones_like(newcube.data), newcube.data)

        # Collapse the cube along the horizontal dimensions
        xcoord = grid.extract_dim_coord(cube, 'X')
        ycoord = grid.extract_dim_coord(cube, 'Y')
        cubes.append(newcube.collapsed([xcoord, ycoord], aggregator, **kwargs))

    return cubes
