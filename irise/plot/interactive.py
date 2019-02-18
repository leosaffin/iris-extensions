"""
"""
import string
import matplotlib.pyplot as plt
from irise import interpolate


class CrossSection(object):
    """Object to faciliate interactive cross-section plotting

    Passed to matplotlib.widgets.RectangleSelector as the onselect function
    """

    def __init__(self, plotters, npoints=100):
        self.plotters = plotters
        self.npoints = npoints
        self.letters = string.ascii_uppercase
        self._letter_idx = -1

    @property
    def letter(self):
        self._letter_idx = (self._letter_idx + 1) % len(self.letters)
        return self.letters[self._letter_idx]

    def __call__(self, click, release):
        # Get the xy points at each end of the line
        xs = click.xdata
        xf = release.xdata
        ys = click.ydata
        yf = release.ydata

        # Plot the cross section on the initial plot
        plt.plot([xs, xf], [ys, yf], '-x', color='k', linewidth=2)
        # Add the coordinate information to the plot
        plt.text(xs, ys, self.letter)
        plt.text(xf, yf, self.letter)

        # Plot the cross section
        plt.figure()
        for plotter in self.plotters:
            plotter.plot(xs, xf, ys, yf, self.npoints)

        plt.show()


class CrossSectionPlotter(object):

    def __init__(self, newfig, plotfunc, cube, *args, **kwargs):
        self.newfig = newfig
        self.plotfunc = plotfunc
        self.cube = cube
        self.args = args
        self.kwargs = kwargs

    def plot(self, xs, xf, ys, yf, npoints):
        if self.newfig:
            plt.figure()
        # Interpolate to each point along the line
        newcube = interpolate.cross_section(self.cube, xs, xf, ys, yf, npoints)

        # Apply the plotting function to the interpolated cube
        self.plotfunc(newcube, *self.args, **self.kwargs)
