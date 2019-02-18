from string import ascii_lowercase
import numpy as np
import matplotlib.pyplot as plt
import iris.plot as iplt


class PlotParameter(object):
    """Acts as a struct type object for holding specific parameters for plotting

    Args:
            color (str): See :py:func:`matplotlib.pyplot.plot`
            linestyle (str): See :py:func:`matplotlib.pyplot.plot`
            label (str): See :py:func:`matplotlib.pyplot.plot`
            idx (int): Index for determining ordering of separate variables.
                Such as ordering in a legend using :py:func:`legend`

    Example:
        I want to plot the evolution of temperature (T) and
        specific humidity (q) with time at some point. I always want to use the
        same colours for these variables to make it easier when looking between
        different plots. Colours that make sense, such as red for temperature
        and blue for humidity. So I create PlotParameter objects for both
        variables that can be imported in my plotting scripts.

        >>> temperature = PlotParameter(color='r', label=r'$T$', idx=1)
        >>> specific_humidity = PlotParameter(color='b', label=r'$q$', idx=1)

        I could also put these objects in a dictionary so that each
        PlotParameter object could be accessed by using the variable name in an
        :py:class:`iris.cube.Cube`.
    """

    def __init__(self, color='k', linestyle='-', label=r'$x$', idx=0):
        self.color = color
        self.linestyle = linestyle
        self.label = label
        self.idx = idx

        return

    def __repr__(self):
        string = 'Color: {}\n' \
                 'Linestyle: {}\n' \
                 'Label: {}\n' \
                 'Index: {}' \
                 ''.format(self.color, self.linestyle, self.label, self.idx)
        return string

    def _add_to_kwargs(self, kwargs):
        # Add saved parameters to keyword arguments
        for parameter in ['color', 'linestyle', 'label']:
            if parameter not in kwargs:
                kwargs[parameter] = getattr(self, parameter)
        return kwargs

    def plot(self, *args, **kwargs):
        """Call to :py:func:`iris.plot.plot` with parameters automatically added

        Same arguments as :py:func:`iris.plot.plot`. Any keyword arguments
        passed that are named the same as parameters stored in this object will
        overwrite those parameters for this plot.
        """
        self._add_to_kwargs(kwargs)
        iplt.plot(*args, **kwargs)

        return


def legend(ax=None, key=None, **kwargs):
    """Place a sorted legend without errorbar markers
    """
    if ax is None:
        ax = plt.gca()

    # Get handles
    handles, labels = ax.get_legend_handles_labels()

    # Remove errorbars
    try:
        handles = [h[0] for h in handles]
    except TypeError:
        # No errorbars just continue
        pass

    # Put the labels in a specified order
    if key is not None:
        labels, handles = zip(*sorted(zip(labels, handles), key=key))

    # Place legend
    ax.legend(handles, labels, **kwargs)

    return


def multilabel(axis, n, factor=0.03, **kwargs):
    """Mark a plot in a multipanel figure with a letter

    Args:
        axis (matplotlib.axes.Axes): Axis within a multi-panel figure

        n (int): The n'th plot

        factor (scalar): Place the label slightly above the axis
    """
    axis.text(0, 1 + factor, '(' + ascii_lowercase[n] + ')',
              transform=axis.transAxes, **kwargs)

    return


def even_cscale(value, levels=17):
    """Creates an even colourscale that skips over zero

    Args:
        value (int): The +/- value at each extreme of the scale
        levels (int): The number of levels in the colourscale. Must be odd.
            Default is 17

    Returns:
        list:
            Numbers from -value to value with even spacing but skipping zero
    """
    if (levels % 2) == 0:
        raise ValueError('Number of levels must be odd')
    cscale = list(np.linspace(-value, value, levels))
    cscale.remove(0)
    return cscale
