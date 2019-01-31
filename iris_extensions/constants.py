"""
Useful scientific constants defined as cubes
"""

from scipy import constants
from iris.cube import Cube


# Ideal gas constant
R_star=Cube(constants.R, units='J mol-1 K-1')

# Specific Gas Constant for  dry air
R = Cube(287.058, units='J kg-1 K-1')

# Specific Gas Constant for water vapour
R_v = Cube(461.5, units='J kg-1 K-1')

# Effective molar mass of dry air
mu_d = Cube(28.97, units='g mol-1')

# Molar mass of water vapour
mu_v = Cube(18.015, units='g mol-1')

# Heat capacity of air at constant volume
c_v = Cube(718.0, units='J kg-1 K-1')

# Heat capacity of air at constant pressure
c_p = Cube(1006.0, units='J kg-1 K-1')

# Enthalpy of vaporisation of water
L = Cube(2.501e6, units='J kg-1')

# Standard atmospheric pressure
P_0 = Cube(constants.bar, units='Pa')

# Standard pressure at Earth's surface
P_atm = Cube(constants.atm, units='Pa')

# Earth's rotation rates
Omega = Cube(7.292116e-5, units='s-1')

# Radius of Earth
earth_radius = Cube(6.371e6, units='m')

# Gravitational acceleration
g = Cube(constants.g, units='m s-2')

# Ratio of radians to degrees
radians_per_degree = Cube(constants.degree, units='radians degrees-1')
