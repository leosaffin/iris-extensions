import iris.util
import numpy as np
from iris.analysis.calculus import differentiate
from irise import constants, grid, interpolate


def multidim(y, x, axis):
    """Calculate a derivative against a multi dimensional coordinate

    The function :py:func:`iris.analysis.calculus.differentiate` is a useful and
    efficient function for calculating a centred finite difference derivative of
    a multidimensional array with respect to a coordinate spanning a single
    dimension.

    This function is a way of extending that function to a coordinate spanning
    more than one dimension. An example is if the cube is in terms of model
    levels but we want the vertical derivative with respect to pressure. Simply
    by calculating the derivative of the cube and pressure with respect to
    model levels then applying the chain rule we get the derivative of the cube
    with respect to pressure at the midpoints in the vertical for all
    grid-points.

    args:
        y, x (iris.cube.Cube): The two multidimensional coordinates to
            calculate the derivative. Must be on the same grid and have a
            one dimensional coordinate for the specified axis. Alternatively x
            can be a string corresponding to a 3d coordinate in y

        axis (str): A letter marking the axis to differentiate with respect to
            (x,y,z,t) or the name of a coordinate on one of these axes.

    returns:
        iris.cube.Cube: The input cube y differentiated with respect to x along
        the chosen axis.
    """
    # Extract x from y as cube if it is specified as a string
    if type(x) == str:
        x = grid.make_cube(y, x)

    # Get the coordinate for the axis to be differentiated over
    try:
        k = y.coord(axis=axis, dim_coords=True)
    except ValueError:
        # Allow referencing the coordinate by name
        k = y.coord(axis)

    # Calculate derivatives with respect to single dimensional coordinate
    dy_dk = differentiate(y, k)
    dx_dk = differentiate(x, k)

    # Apply the chain rule
    dy_dx = dy_dk / dx_dk

    return dy_dx


def diff_by_axis(cube, axis):
    """Differentiate a cube with respect to a single dimension

    args:
        cube (iris.cube.Cube):

        axis (str): A letter marking the axis to differentiate with respect to
            (x,y,z,t).

    returns:
        iris.cube.Cube: The derivative of the cube with respect to its
        dimensional coordinate on the chosen axis
    """
    # Extract the dimensional coordinate corresponding to the axis
    coord = cube.coord(axis=axis, dim_coords=True)

    # Calculate the derivative
    diff = differentiate(cube, coord)

    return diff


def polar_horizontal(cube, axis):
    r"""Differentiate a cube with respect to lon/lat and convert to x/y


    :math:`dx = r cos(\phi) d\lambda`

    :math:`dy = r d\phi`

    args:
        cube (iris.cube.Cube):

        axis (str): A letter marking the axis to differentiate with respect to
            (x,y).

    returns:
        iris.cube.Cube: The derivative of the cube along the chosen polar axis
    """
    # Calculate derivative with respect to polar axis in radians
    diff = diff_by_axis(cube, axis) / constants.radians_per_degree

    # Convert the differential to x/y by considering the distance in lon/lat
    # Coordinates
    # Calculate radius relative to Earth centre
    radius = grid.make_cube(diff, 'altitude')
    radius.data += constants.earth_avg_radius.data

    if radius.ndim == 1:
        radius = iris.util.broadcast_to_shape(radius.data, diff.shape, [0])
        radius = diff.copy(data=radius)
        radius.rename("altitude")
        radius.units = "m"

    if axis.lower() == 'x':
        lat = (diff.coord(axis='y').points *
               constants.radians_per_degree.data)
        lat = np.outer(lat, np.ones(diff.shape[-1]))
        metres_per_radian = radius * np.cos(lat)

    elif axis.lower() == 'y':
        metres_per_radian = radius
    else:
        raise ValueError('Can only specify x or y axis')

    metres_per_radian.units = 'm radian-1'
    diff = diff / metres_per_radian
    return diff


def div(u, v, w):
    r"""Calculate divergence of the input vector

    :math:`\nabla \cdot \mathbf{x}`

    Args:
        u,v,w (iris.cube.Cube): i,j,k components of the input vector

    Returns:
        iris.cube.Cube: The divergence
    """
    # Calculate individual components
    du_dx = polar_horizontal(u, 'x')
    du_dx = interpolate.remap_3d(du_dx, w)

    dv_dy = polar_horizontal(v, 'y')
    dv_dy = interpolate.remap_3d(dv_dy, w)

    dw_dz = multidim(w, 'altitude', 'z')
    dw_dz = interpolate.remap_3d(dw_dz, w)

    # Sum to divergence
    divergence = du_dx.data + dv_dy.data + dw_dz.data

    # Copy to a cube
    divergence = dw_dz.copy(data=divergence)

    return divergence


def grad(cube):
    r"""Calculate gradient of the input field

    :math:`\nabla x`

    Args:
        cube (iris.cube.Cube): Input field

    Returns:
        tuple (iris.cube.Cube, iris.cube.Cube, iris.cube.Cube):
            Three cubes of the different components of the vector gradient
    """
    d_dx = polar_horizontal(cube, 'x')
    d_dx = interpolate.remap_3d(d_dx, cube)

    d_dy = polar_horizontal(cube, 'y')
    d_dy = interpolate.remap_3d(d_dy, cube)

    d_dz = multidim(cube, 'altitude', 'z')
    d_dz = interpolate.remap_3d(d_dz, cube)

    return d_dx, d_dy, d_dz


def curl(u, v, w):
    r"""Calculate curl of the input vector

    :math:`\nabla \times \mathbf{x}`

    Args:
        u,v,w (iris.cube.Cube): i,j,k components of the input vector

    Returns:
        tuple (iris.cube.Cube, iris.cube.Cube, iris.cube.Cube):
            Three cubes of the different components of the vector curl
    """
    # Calculate individual gradients
    dw_dx = polar_horizontal(w, 'x')
    dw_dx = interpolate.remap_3d(dw_dx, w)

    dw_dy = polar_horizontal(w, 'y')
    dw_dy = interpolate.remap_3d(dw_dy, w)

    du_dz = multidim(u, 'altitude', 'z')
    du_dz = interpolate.remap_3d(du_dz, w)

    dv_dz = multidim(v, 'altitude', 'z')
    dv_dz = interpolate.remap_3d(dv_dz, w)

    du_dy = polar_horizontal(u, 'y')
    du_dy = interpolate.remap_3d(du_dy, w)

    dv_dx = polar_horizontal(v, 'x')
    dv_dx = interpolate.remap_3d(dv_dx, w)

    # Calculate the components of vorticity
    xi_i = dw_dy.data - dv_dz.data
    xi_i = dw_dx.copy(data=xi_i)

    xi_j = du_dz.data - dw_dx.data
    xi_j = dw_dx.copy(data=xi_j)

    xi_k = dv_dx.data - du_dy.data
    xi_k = dw_dx.copy(data=xi_k)

    return xi_i, xi_j, xi_k


def laplacian(cube):
    r"""Calculate the Laplacian of the input field

    :math:`\nabla^2 x = \nabla \cdot \nabla x`

    Args:
        cube (iris.cube.Cube): Input field

    Returns:
        iris.cube.Cube: Laplacian of the input field
    """
    return div(*grad(cube))
