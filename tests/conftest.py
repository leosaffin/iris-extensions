import datetime

import pytest
import iris
from iris.tests import stock

from irise import grid
from irise.forecast import Forecast


@pytest.fixture
def forecast():
    """Creates a set of files with different lead times
    """
    cube = stock.simple_4d_with_hybrid_height()
    times = grid.get_datetime(cube)

    mapping = {}
    for n, time in enumerate(times):
        filename = 'forecast_' + str(n) + '.nc'
        iris.save(cube[n], filename)
        mapping[
            datetime.datetime(
                time.year,
                time.month,
                time.day,
                time.hour,
                time.minute,
            )
        ] = filename

    yield Forecast(sorted(mapping.keys())[0], mapping)


@pytest.fixture(scope="session", params=[
    stock.realistic_3d(),
    stock.realistic_4d(),
])
def cube(request):
    return request.param


@pytest.fixture
def realistic_3d():
    yield stock.realistic_3d()


@pytest.fixture
def realistic_4d():
    yield stock.realistic_4d()
