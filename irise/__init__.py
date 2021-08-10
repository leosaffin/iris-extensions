# -*- coding: utf-8 -*-

"""Top-level package for iris-extensions."""

__author__ = """Leo Saffin"""
__email__ = 'string_buster@hotmail.com'
__version__ = '0.2.0'

import iris
from iris.fileformats.pp import STASH
from iris.util import squeeze


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
    STASH(model=1, section=0, item=635):
        dict(name="advection_only_q", units="1"),
    STASH(model=1, section=0, item=636):
        dict(name="advection_only_qcl", units="1"),
    STASH(model=1, section=0, item=637):
        dict(name="advection_only_qcf", units="1"),
    STASH(model=1, section=0, item=638):
        dict(name="short_wave_radiation_q", units="1"),
    STASH(model=1, section=0, item=639):
        dict(name="short_wave_radiation_qcl", units="1"),
    STASH(model=1, section=0, item=640):
        dict(name="long_wave_radiation_q", units="1"),
    STASH(model=1, section=0, item=641):
        dict(name="long_wave_radiation_qcl", units="1"),
    STASH(model=1, section=0, item=642):
        dict(name="microphysics_q", units="1"),
    STASH(model=1, section=0, item=643):
        dict(name="microphysics_qcl", units="1"),
    STASH(model=1, section=0, item=644):
        dict(name="microphysics_qcf", units="1"),
    STASH(model=1, section=0, item=645):
        dict(name="pc2_checks_q", units="1"),
    STASH(model=1, section=0, item=646):
        dict(name="pc2_checks_qcl", units="1"),
    STASH(model=1, section=0, item=647):
        dict(name="pc2_checks_qcf", units="1"),
    STASH(model=1, section=0, item=648):
        dict(name="convection_q", units="1"),
    STASH(model=1, section=0, item=649):
        dict(name="convection_qcl", units="1"),
    STASH(model=1, section=0, item=650):
        dict(name="convection_qcf", units="1"),
    STASH(model=1, section=0, item=651):
        dict(name="boundary_layer_q", units="1"),
    STASH(model=1, section=0, item=652):
        dict(name="boundary_layer_qcl", units="1"),
    STASH(model=1, section=0, item=653):
        dict(name="boundary_layer_qcf", units="1"),
    STASH(model=1, section=0, item=654):
        dict(name="cloud_q", units="1"),
    STASH(model=1, section=0, item=655):
        dict(name="cloud_qcl", units="1"),
    STASH(model=1, section=0, item=656):
        dict(name="cloud_qcf", units="1"),
    STASH(model=1, section=0, item=657):
        dict(name="pc2_erosion_q", units="1"),
    STASH(model=1, section=0, item=658):
        dict(name="pc2_erosion_qcl", units="1"),
    STASH(model=1, section=0, item=659):
        dict(name="slow_physics_q", units="1"),
    STASH(model=1, section=0, item=660):
        dict(name="slow_physics_qcl", units="1"),
    STASH(model=1, section=0, item=661):
        dict(name="slow_physics_qcf", units="1"),
    STASH(model=1, section=0, item=662):
        dict(name="advection_correction_q", units="1"),
    STASH(model=1, section=0, item=663):
        dict(name="advection_correction_qcl", units="1"),
    STASH(model=1, section=0, item=664):
        dict(name="advection_correction_qcf", units="1"),
    STASH(model=1, section=0, item=665):
        dict(name="solver_q", units="1"),
    STASH(model=1, section=0, item=666):
        dict(name="solver_qcl", units="1"),
    STASH(model=1, section=0, item=667):
        dict(name="solver_qcf", units="1"),
    STASH(model=1, section=0, item=668):
        dict(name="methane_oxidation_q", units="1"),
    STASH(model=1, section=0, item=669):
        dict(name="microphysics_settling_q", units="1"),
    STASH(model=1, section=0, item=670):
        dict(name="microphysics_fixes_q", units="1"),
    STASH(model=1, section=0, item=671):
        dict(name="microphysics_nucleation_q", units="1"),
    STASH(model=1, section=0, item=672):
        dict(name="microphysics_deposition_q", units="1"),
    STASH(model=1, section=0, item=673):
        dict(name="microphysics_evaporation_q", units="1"),
    STASH(model=1, section=0, item=674):
        dict(name="microphysics_melting_q", units="1"),
    STASH(model=1, section=0, item=675):
        dict(name="boundary_layer_entrainment_q", units="1"),
    STASH(model=1, section=0, item=676):
        dict(name="boundary_layer_surface_fluxes_q", units="1"),
    STASH(model=1, section=0, item=677):
        dict(name="boundary_layer_other_q", units="1"),
    STASH(model=1, section=0, item=678):
        dict(name="fast_physics_q", units="1"),
    STASH(model=1, section=0, item=679):
        dict(name="theta_perturbations_q", units="1"),
    STASH(model=1, section=0, item=680):
        dict(name="leonard_terms_q", units="1"),
}
