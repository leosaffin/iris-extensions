# -*- coding: utf-8 -*-

"""Top-level package for iris-extensions."""

__author__ = """Leo Saffin"""
__email__ = 'string_buster@hotmail.com'
__version__ = '0.2.0'

import iris
from iris.fileformats.pp import STASH
from iris.util import squeeze

from . import calculus, constants, convert, forecast, grid, interpolate, variable


def load(uris, constraints=None, callback=None):
    cubes = iris.load(uris, constraints=constraints, callback=callback)
    fix_by_stash(cubes)
    cubes.sort(key=get_stash)

    return iris.cube.CubeList([squeeze(cube) for cube in cubes])


def fix_by_stash(cubes):
    for cube in cubes:
        stash = get_stash(cube)
        if stash in stash_map:
            cube.rename(stash_map[stash]["name"])
            cube.units = stash_map[stash]["units"]


def get_stash(cube):
    try:
        return cube.attributes["STASH"]
    except KeyError:
        return STASH(1, 0, 1)


stash_map = {
    STASH(model=1, section=0, item=389): dict(name="air_density", units="kg m-3"),

    # PV tracers
    STASH(model=1, section=0, item=593):
        dict(name="advection_only_pv", units="m2 K s-1 kg-1"),
    STASH(model=1, section=0, item=577):
        dict(name="radiation_pv", units="m2 K s-1 kg-1"),
    STASH(model=1, section=0, item=578):
        dict(name="short_wave_radiation_pv", units="m2 K s-1 kg-1"),
    STASH(model=1, section=0, item=579):
        dict(name="long_wave_radiation_pv", units="m2 K s-1 kg-1"),
    STASH(model=1, section=0, item=580):
        dict(name="microphysics_pv", units="m2 K s-1 kg-1"),
    STASH(model=1, section=0, item=581):
        dict(name="gravity_wave_drag_pv", units="m2 K s-1 kg-1"),
    STASH(model=1, section=0, item=582):
        dict(name="slow_physics_pv", units="m2 K s-1 kg-1"),
    STASH(model=1, section=0, item=583):
        dict(name="convection_pv", units="m2 K s-1 kg-1"),
    STASH(model=1, section=0, item=584):
        dict(name="boundary_layer_pv", units="m2 K s-1 kg-1"),
    STASH(model=1, section=0, item=585):
        dict(name="stochastic_physics_pv", units="m2 K s-1 kg-1"),
    STASH(model=1, section=0, item=586):
        dict(name="cloud_pv", units="m2 K s-1 kg-1"),
    STASH(model=1, section=0, item=587):
        dict(name="analyis_increment_pv", units="m2 K s-1 kg-1"),
    STASH(model=1, section=0, item=588):
        dict(name="nudging_pv", units="m2 K s-1 kg-1"),
    STASH(model=1, section=0, item=589):
        dict(name="total_pv", units="m2 K s-1 kg-1"),
    STASH(model=1, section=0, item=590):
        dict(name="dynamics_tracer_inconsistency", units="m2 K s-1 kg-1"),
    STASH(model=1, section=0, item=591):
        dict(name="solver_pv", units="m2 K s-1 kg-1"),
    STASH(model=1, section=0, item=592):
        dict(name="mass_update_pv", units="m2 K s-1 kg-1"),

    # Potential temperature tracers
    STASH(model=1, section=0, item=600):
        dict(name="advection_only_theta", units="K"),
    STASH(model=1, section=0, item=601):
        dict(name="boundary_layer_theta", units="K"),
    STASH(model=1, section=0, item=602):
        dict(name="boundary_layer_mixing_theta", units="K"),
    STASH(model=1, section=0, item=603):
        dict(name="boundary_layer_heating_theta", units="K"),
    STASH(model=1, section=0, item=604):
        dict(name="convection_theta", units="K"),
    STASH(model=1, section=0, item=605):
        dict(name="microphysics_theta", units="K"),
    STASH(model=1, section=0, item=606):
        dict(name="radiation_theta", units="K"),
    STASH(model=1, section=0, item=607):
        dict(name="short_wave_radiation_theta", units="K"),
    STASH(model=1, section=0, item=608):
        dict(name="long_wave_radiation_theta", units="K"),
    STASH(model=1, section=0, item=609):
        dict(name="slow_physics_theta", units="K"),
    STASH(model=1, section=0, item=610):
        dict(name="cloud_theta", units="K"),
    STASH(model=1, section=0, item=611):
        dict(name="pc2_checks_theta", units="K"),
    # Moisture tracers
    STASH(model=1, section=0, item=701):
        dict(name="advection_only_q", units="1"),
    STASH(model=1, section=0, item=702):
        dict(name="short_wave_radiation_q", units="1"),
    STASH(model=1, section=0, item=703):
        dict(name="long_wave_radiation_q", units="1"),
    STASH(model=1, section=0, item=704):
        dict(name="microphysics_q", units="1"),
    STASH(model=1, section=0, item=705):
        dict(name="convection_q", units="1"),
    STASH(model=1, section=0, item=706):
        dict(name="boundary_layer_q", units="1"),
    STASH(model=1, section=0, item=707):
        dict(name="microphysics_cloud_q", units="1"),
    STASH(model=1, section=0, item=708):
        dict(name="PC2_checks_q", units="1"),
    STASH(model=1, section=0, item=709):
        dict(name="boundary_layer_cloud_q", units="1"),
    STASH(model=1, section=0, item=710):
        dict(name="cloud_q", units="1"),
    STASH(model=1, section=0, item=711):
        dict(name="methane_oxidation_q", units="1"),
    STASH(model=1, section=0, item=712):
        dict(name="leonard_terms_q", units="1"),
    STASH(model=1, section=0, item=713):
        dict(name="perturbations_q", units="1"),
    STASH(model=1, section=0, item=714):
        dict(name="rain_evaporation_q", units="1"),
}
