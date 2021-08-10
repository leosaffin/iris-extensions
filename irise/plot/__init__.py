"""Generic plotting functions
"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
import cartopy.crs as ccrs
import iris
import iris.plot as iplt
import iris.quickplot as qplt

from irise import grid, interpolate


def colored_line_plot(x, y, color, vmin=None, vmax=None, cmap='gray'):
    """Add a multicolored line to an existing plot

    Args:
        x (np.array): The x points of the plot

        y (np.array): The y points of the plot

        color (np.array): The color of the line at the xy points

        vmin (scalar, optional): The minimum of the colorscale. Defaults to the
            minimum of the color array.

        vmax (scalar, optional): The maximum of the colorscale. Defaults to the
            maximum of the color array.

        cmap (str, optional): Colormap to plot. Default is grey.

    returns:
        matplotlib.collections.LineCollection:
            The plotted LineCollection. Required as argument to
            :py:func:`matplotlib.pyplot.colorbar`
    """
    # Set the color scalings
    if vmin is None:
        vmin = color.min()
    if vmax is None:
        vmax = color.max()

    # Break the xy points up in to line segments
    points = np.array([x, y]).T.reshape(-1, 1, 2)
    segments = np.concatenate([points[:-1], points[1:]], axis=1)

    # Collect the line segments
    lc = LineCollection(segments, linewidth=1, cmap=plt.get_cmap(cmap),
                        norm=plt.Normalize(vmin, vmax))

    # Set the line color to the specified array
    lc.set_array(color)

    # Add the colored line to the existing plot
    plt.gca().add_collection(lc)

    return lc


def contour(cube, *args, **kwargs):
    """
    Plot a variable on a single height level automatically adding title,
    coastlines and gridlines

    Parameters:
        cube (iris.cube.Cube): The variable to be plotted. Must be 2D

        *args: Variable length arguments passed to
            :py:func:`iris.quickplot.contour`.

    Keyword Args:
        pv (iris.cube.Cube):
        projection (cartopy.crs.Projection or str):
        **kwargs: Remaining keyword arguments are passed to
            :py:func:`iris.quickplot.contour`.
    """
    return _plot_3d(cube, qplt.contour, *args, **kwargs)


def contourf(cube, *args, **kwargs):
    """
    Plot a variable on a single height level automatically adding colorbar,
    title, coastlines and gridlines

    Args:
        cube (iris.cube.Cube): The variable to be plotted. Must be 2D
    """
    return _plot_3d(cube, qplt.contourf, *args, **kwargs)


def pcolormesh(cube, *args, **kwargs):
    """
    Plot a variable on a single height level automatically adding colorbar,
    title, coastlines and gridlines

    Args:
        cube (iris.cube.Cube): The variable to be plotted. Must be 2D
    """
    result = _plot_3d(cube, qplt.pcolormesh, *args, **kwargs)

    return result


def errorbar(*args, **kwargs):
    """
    Args:
        *args: Supplies the cubes and coords to make the basic plot in the same
            way as iris.plot.plot

    Keyword Args:
        xerr (iris.cube.Cube or iris.coords.Coord): Can replace usual argument
            type to :py:func:`matplotlib.pyplot.errorbar`. The scalar/array will
            be extracted from the cube or coord.
        yerr (iris.cube.Cube or iris.coords.Coord): Can replace usual argument
            type to :py:func:`matplotlib.pyplot.errorbar`. The scalar/array will
            be extracted from the cube or coord.
        **kwargs: Remaining keyword arguments to
            :py:func:`matplotlib.pyplot.errorbar` can be specified as normal
    """
    if 'coords' in kwargs:
        raise TypeError('"coords" is not a valid plot keyword. Coordinates '
                        'and cubes may be passed as arguments for '
                        'full control of the plot axes.')

    # Allow errorbars to be specified as cubes or coordinates
    for err in ['xerr', 'yerr']:
        if err in kwargs:
            try:
                kwargs[err] = iplt._data_from_coord_or_cube(kwargs[err])
            except TypeError:
                pass

    result = iplt._draw_1d_from_points('errorbar', None, *args, **kwargs)
    qplt._label_1d_plot(*args)

    return result


def multiline(cubelist, xcoord=None, ycoord=None, legend=True,
              with_units=False, **kwargs):
    """Plot multiple cubes on the same axis

    Args:
        cubelist (iris.cube.CubeList): A list of 1 dimensional cube

        xcoord:

        ycoord:

        legend:

        with_units:
    """
    for cube in cubelist:
        # Add the cube label to the arguments for each plot
        kwargs['label'] = qplt._title(cube, with_units)

        # Determine which way to plot coord vs cube
        if xcoord is None:
            if ycoord is None:
                iplt.plot(cube, **kwargs)
            else:
                iplt.plot(cube, ycoord, **kwargs)
        else:
            iplt.plot(cube, xcoord, **kwargs)

    # Create a second figure containing the legend
    ax = plt.gca()
    if legend is True:
        legend = plt.figure()
        plt.figlegend(*ax.get_legend_handles_labels(), loc='upper left')

    return ax, legend


def overlay_winds(u, v, nx, ny, **kwargs):
    """Overlay a quiver plot on an existing iris plot

    Args:
        u (iris.cube.Cube): The x-component of the vector to be plotted

        v (iris.cube.Cube): The y-component of the vector to be plotted

        nx (int): The coarse grained resolution on the x-axis

        ny (int): The coarse grained resolution on the y-axis

    """
    # Extract a coarse representation of the x and y coordinates
    xcoord = u.coord(axis='x', dim_coords=True)
    x_fac = int(len(xcoord.points) / nx)
    ycoord = u.coord(axis='y', dim_coords=True)
    y_fac = int(len(ycoord.points) / ny)
    x_coarse = xcoord.points[::x_fac]
    y_coarse = ycoord.points[::y_fac]

    # Interpolate the vector field to the coarse resolution
    interp_kwargs = {xcoord.name(): x_coarse, ycoord.name(): y_coarse}
    u_coarse = interpolate.interpolate(u, **interp_kwargs)
    v_coarse = interpolate.interpolate(v, **interp_kwargs)

    # Plot the coarse data
    return plt.quiver(
        x_coarse, y_coarse, u_coarse.data, v_coarse.data,**kwargs)


def _plot_3d(cube, func, *args, **kwargs):
    # Create a figure with the requested projection
    kwargs = _setup_projection(**kwargs)

    # Add a contour of PV=2 if requested
    kwargs = _add_pv2(**kwargs)

    # Plot the variable
    result = func(cube, *args, **kwargs)

    # Add coastlines and gridlines
    _add_map()

    return result


def _setup_projection(**kwargs):
    # Remove projection from kwargs
    if 'projection' in kwargs:
        projection = kwargs.pop('projection')

        # Convert name of projection to projection instance
        if type(projection) == str:
            projection = eval('ccrs.' + projection + '()')

        plt.subplot(projection=projection)

    return kwargs


def _add_map():
    # Adds coastlines and gridlines if plot has cartopy projection
    try:
        # Add coastlines
        plt.gca().coastlines()
        # Add gridlines
        plt.gca().gridlines()
    except AttributeError:
        pass


def _add_pv2(**kwargs):
    # Plots the 2 PVU tropopause if contained in kwargs
    # Remove pv from kwargs
    if 'pv' in kwargs:
        pv = kwargs.pop('pv').copy()

        pv.convert_units('PVU')
        # Allow PV=2 to be plotted for both hemispheres
        pv.data = np.abs(pv.data)

        iplt.contour(pv, [2], colors='k', linewidths=2)

    return kwargs
