"""Creates files used for testing when loading is required
"""

import datetime
import pickle
import iris
from iris.tests import stock


def forecast():
    """Creates a set of files with different lead times
    """
    start_time = datetime.datetime(2013, 12, 1)
    dt = datetime.timedelta(hours=1)
    cube = stock.simple_4d_with_hybrid_height()
    nt = len(cube.coord('time').points)
    mapping = {}
    for n in range(nt):
        filename = 'forecast_' + str(n) + '.nc'
        subcube = cube[n]
        iris.save(subcube, filename)
        mapping[start_time + n * dt] = filename

    with open('forecast_times.pkl', 'wb') as output:
        pickle.dump((start_time, mapping), output)


if __name__ == '__main__':
    forecast()
