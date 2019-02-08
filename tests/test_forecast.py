"""
"""

import pickle
import unittest
import datetime
from irise import forecast


class TestForecast(unittest.TestCase):

    def setUp(self):
        with open('forecast_times.pkl', 'rb') as infile:
            start_time, mapping = pickle.load(infile)
        self.forecast = forecast.Forecast(start_time, mapping)

    def tearDown(self):
        del self.forecast

    def test_deletes(self):
        # Load three files, when the third is loaded it should delete the first
        self.forecast.set_lead_time(datetime.timedelta(hours=0))
        self.forecast.set_lead_time(datetime.timedelta(hours=1))
        self.forecast.set_lead_time(datetime.timedelta(hours=2))

        # Check the number of files are equal to the default limit of 2
        nfiles = len(self.forecast._loader._loaded)
        self.assertEqual(nfiles, 2)

        # Check the loaded files are for lead time 1 and 2
        self.assertTrue(self.forecast.start_time + datetime.timedelta(hours=1)
                        in self.forecast._loader._loaded)
        self.assertTrue(self.forecast.start_time + datetime.timedelta(hours=2)
                        in self.forecast._loader._loaded)

    def test_out_of_bounds(self):
        with self.assertRaises(KeyError):
            self.forecast.set_lead_time(datetime.timedelta(hours=10))

        with self.assertRaises(KeyError):
            self.forecast.set_lead_time(datetime.timedelta(hours=-1))
