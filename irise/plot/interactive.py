"""
"""
import string
import matplotlib.pyplot as plt
from irise import interpolate


class CrossSection(object):
    """Object to faciliate interactive cross-section plotting

    Draw a cross section on the plot by clicking on it to draw a line. Line starts
    drawing on click and is finalised when the click it released. A line will be drawn
    on the plot showing where the cross section was made and a new set of figures
    defined by the "plotters" will be drawn for the cross section
    """
    def __init__(self, fig, ax, plotters=[], npoints=100):
        self.fig = fig
        self.ax = ax
        self.plotters = plotters
        self.npoints = npoints
        self.letters = string.ascii_uppercase
        self._letter_idx = -1

        self.cid = fig.canvas.mpl_connect('button_press_event', self.press)
        self.cid = fig.canvas.mpl_connect('motion_notify_event', self.move)
        self.cid = fig.canvas.mpl_connect('button_release_event', self.release)

        self.line = None
        self.drawing = False
        self.start_x = 0.0
        self.start_y = 0.0

    @property
    def letter(self):
        self._letter_idx = (self._letter_idx + 1) % len(self.letters)
        return self.letters[self._letter_idx]

    def press(self, event):
        self.drawing = True
        self.start_x = event.xdata
        self.start_y = event.ydata

        self.ax.plot(self.start_x, self.start_y, "kx")
        self.line, = self.ax.plot(self.start_x, self.start_y)

    def move(self, event):
        if self.drawing:
            self.line.remove()
            self.line, = self.ax.plot([self.start_x, event.xdata], [self.start_y, event.ydata], "--k")

            self.fig.canvas.draw()

    def release(self, event):
        self.drawing = False

        self.ax.plot(event.xdata, event.ydata, "kx")
        self.line.remove()
        self.ax.plot([self.start_x, event.xdata], [self.start_y, event.ydata], "--k")
        plt.text(self.start_x, self.start_y, self.letter)
        plt.text(event.xdata, event.ydata, self.letter)

        self.fig.canvas.draw()

        # Plot the cross section
        print(self.start_x, event.xdata, self.start_y, event.ydata)
        for plotter in self.plotters:
            plotter.plot(
                self.start_x,
                event.xdata,
                self.start_y,
                event.ydata,
                self.npoints
            )

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
        # TODO account for rotated grids
        newcube = interpolate.cross_section(self.cube, xs, xf, ys, yf, npoints)

        # Apply the plotting function to the interpolated cube
        self.plotfunc(newcube, *self.args, **self.kwargs)
