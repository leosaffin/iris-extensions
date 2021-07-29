"""
Useful scientific constants defined as iris cubes

Uses constants defined by metpy.constansts
"""

from scipy import constants
from metpy import constants as mpconst
from iris.cube import Cube


# Ratio of radians to degrees
radians_per_degree = Cube(constants.degree, units='radians degrees-1')


# Import all constants from metpy and turn them into cubes
for constant in mpconst.__all__:
    x = mpconst.__dict__[constant]
    units = str(x.units).replace(" ", "").replace("dimensionless", "")
    exec(constant + " = Cube("
         "data=x.magnitude,"
         "long_name=constant,"
         "units=units)"
         )

    # Some units are in grams while other are in kilograms
    if "gram" in units and "kilogram" not in units:
        fixed_units = units.replace("gram", "kilogram")

    # Use SI units for pressure variables
    elif "millibar" in units:
        fixed_units = units.replace("millibar", "pascal")

    else:
        fixed_units = None

    if fixed_units is not None:
        exec("{}.convert_units('{}')".format(constant, fixed_units))
