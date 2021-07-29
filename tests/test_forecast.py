import datetime

import pytest


def test_deletes(forecast):
    print(forecast._loader.files)
    # Load three files, when the third is loaded it should delete the first
    forecast.set_lead_time(datetime.timedelta(hours=0))
    forecast.set_lead_time(datetime.timedelta(hours=1))
    forecast.set_lead_time(datetime.timedelta(hours=2))

    # Check the number of files are equal to the default limit of 2
    nfiles = len(forecast._loader._loaded)
    assert nfiles == 2

    # Check the loaded files are for lead time 1 and 2
    assert forecast.start_time + datetime.timedelta(hours=1) \
           in forecast._loader._loaded
    assert forecast.start_time + datetime.timedelta(hours=2) \
           in forecast._loader._loaded


def test_out_of_bounds(forecast):
    with pytest.raises(KeyError):
        forecast.set_lead_time(datetime.timedelta(hours=10))

    with pytest.raises(KeyError):
        forecast.set_lead_time(datetime.timedelta(hours=-1))
