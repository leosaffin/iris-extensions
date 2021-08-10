"""Methods for handling large amounts of forecast data in multiple files
"""
import copy
import datetime
import numpy as np
import matplotlib.pyplot as plt
import iris

import irise


class Forecast(object):
    """A collection of data from a single forecast

    Args:
        start_time (datetime.datetime): The time of the initialisation of
            the forecast
        mapping (dict): A mapping of `datetime.datetime` to filenames

    Attributes:
        start_time (datetime.datetime): Same as input
        current_time (datetime.datetime): The currently loaded lead time
        cubelist (iris.cube.CubeList): The currently loaded lead time data
    """

    def __init__(self, start_time, mapping):
        self.start_time = start_time
        self.current_time = None
        self.cubelist = None
        self._loader = _CubeLoader(mapping)

    def set_lead_time(self, *args, **kwargs):
        """Loads forecast data for the given lead time

        Args:
            *args: Can pass a single ``datetime.timedelta`` object as the
                requested lead time

            **kwargs: Alternatively pass the argument to the
                ``datetime.timedelta`` object directly to this function

        Returns:
            iris.cube.CubeList: The cubelist at the specified lead time
        """
        # Extract the lead time from the input
        if len(args) == 0:
            lead_time = datetime.timedelta(**kwargs)
        elif len(args) == 1:
            lead_time = args[0]
        else:
            raise ValueError('Can only pass one lead time')

        newtime = self.start_time + lead_time
        return self.set_time(newtime)

    def set_time(self, newtime):
        """Loads forecast data for the given time

        Args:
            newtime (datetime.datetime):

        Returns:
            iris.cube.CubeList: The cubelist at the specified time
        """
        # Only change the lead time if it's not already loaded
        if newtime != self.current_time:
            self.cubelist = self._loader.extract(newtime)
            self.current_time = newtime

        return self.cubelist

    @property
    def times(self):
        return sorted(list(self._loader.files))

    def animate(self, func, filename, *args, **kwargs):
        """Produce the same plot for each lead time in the forecast indexed by
        time

        Args:
            func (callable): A function that produces a matplotlib figure
                given an input cubelist

            filename (str): The path of the output figure

            *args: Additional parameters to pass to the plotting function

            **kwargs: Additional parameters to pass to the plotting function
        """
        for n, cubes in enumerate(self):
            # Run the function theta produces the plot
            func(cubes, *args, **kwargs)

            # Put the timestamp on the figure
            plt.title(str(self.current_time))

            # Save the produced figure with the lead time index
            plt.savefig(filename + str(n).zfill(3) + '.png')

    def produce_diagnostics(self, func, filename, *args, **kwargs):
        """Create a netcdf file with derived variables for each lead time

        Args:
            func: A function that produces an iris.cube.CubeList of
                variables derived from an input cubelist

            filename (str): The path of the output netcdf file

            *args: Additional parameters to pass to the function

            **kwargs: Additional parameters to pass to the function
        """
        for n, cubes in enumerate(self):
            # Call the function
            output = func(cubes, *args, **kwargs)

            # Save the derived output
            iris.save(output, filename + '_' + str(n).zfill(3) + '.nc')

        # Merge the output across all timesteps
        cubes = iris.load(filename + '_*.nc')
        iris.save(cubes, filename + '.nc')

    def copy(self):
        return copy.deepcopy(self)

    def __len__(self):
        return len(self._loader.files)

    def __iter__(self):
        for time in self.times:
            yield self.set_time(time)

    def __reversed__(self):
        for time in reversed(self.times):
            yield self.set_time(time)


class _CubeLoader(object):
    """Loads data from files so that a limited number are loaded at once

    Attributes:
        _loaded (dict): Currently loaded lead times paired to cubelists
        files (dict): Pairings of lead time to files
        limit (int): Limit on number of cubelists to be loaded at once
    """

    def __init__(self, mapping, limit=2):
        """

        Args:
            mapping (dict): A dictionary mapping times to filenames
            limit (int, optional): The maximum number of times allowed to be
                loaded at one time, default is 2
        """
        self.limit = limit
        self.files = mapping
        self._loaded = {}

    def extract(self, time):
        """Gets the data corresponding to the given time

        Args:
            time(datetime.datetime): The time of data to extract
        """
        # Return the cubelist if it is already loaded
        if time in self._loaded:
            return self._loaded[time]

        # Load a file if one corresponds to the lead time
        elif time in self.files:
            self._load_new_time(time)
            return self._loaded[time]

        # If the lead time does not match the available files raise an error
        else:
            raise KeyError('Requested time not available in files')

    def _load_new_time(self, time):
        """ Loads a new cubelist and removes others if necessary

        Args:
            time (datetime.datetime): The new time to be loaded
        """
        # Clear space for the new files
        self._make_space(time)

        # Load data from files with that lead time
        cubes = irise.load(self.files[time])

        # Add the data to the loaded files
        self._loaded[time] = cubes

    def _make_space(self, time):
        """ Removes a cubelist if the limit is about to be exceeded

        Args:
            time (datetime.datetime): The time of the new file to be loaded
        """
        # Work out if the limit is going to be exceeded
        while (len(self._loaded) - self.limit) >= 0:
            self._remove_time(time)

    def _remove_time(self, time):
        """ Removes data furthest from the given time
        """
        # Get the index of the loaded files furthest from the requested time
        index = abs(np.array(list(self._loaded)) - time).argmax()
        time = list(self._loaded)[index]
        self._loaded.pop(time)
