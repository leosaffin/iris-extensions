"""Calculating meteorological variables
"""
import iris.cube
import iris.exceptions

from irise import grid, interpolate, variable


def calc(names, cubelist, levels=None):
    """Wrapper function to ``calculate`` with the option to interpolate to a
    level

    Args:
        names (str or list[str]): The CF standard name(s) of the
            variable(s) to be calculated.

        cubelist (iris.cube.CubeList): Contains either the requested variable or
            the variables required to calculate the requested variable.

        levels (tuple or None): The name and levels to interpolate the
            requested variable to. Default is None.

    Returns:
        iris.cube.Cube or iris.cube.CubeList: The variable(s) requested.
        Optionally interpolated to the given level
    """
    if type(names) == str:
        names = [names]
    # Call calculate to get the cube
    cubes = [calculate(name, cubelist) for name in names]

    # Return the cube if interpolation is not requested
    if levels is None:
        if len(names) == 1:
            return cubes[0]
        else:
            return iris.cube.CubeList(cubes)

    else:
        output = iris.cube.CubeList()

        # Extract the coordinate requested
        coord_name, values = levels

        for cube in cubes:
            # Extract the coordinate from the cubes
            try:
                coord = cube.coord(coord_name)

            # Alternatively use a cube from the cubelist
            except iris.exceptions.CoordinateNotFoundError:
                coord = calculate(coord_name, cubelist)
                coord = grid.make_coord(coord)
                cube.add_aux_coord(coord, range(cube.ndim))

            # Interpolate to the requested coordinate levels
            if coord.points.ndim == 1:
                result = cube.interpolate(
                    [(coord_name, values)],
                    iris.analysis.Linear()
                )
            else:
                result = interpolate.to_level(cube, **{coord_name: values})
            output.append(result)

        if len(names) == 1:
            return output[0]
        else:
            return output


def calculate(name, cubelist):
    """Calculates any variable available from the cubelist

    Args:
        name (string): The CF standard name of the variable to be calculated

        cubelist (iris.cube.CubeList): A cubelist containing either the
            requested variable or the variables required to calculate the
            requested variable

    Returns:
        iris.cube.Cube: The variable requested

    Raises:
        ValueError: If the requested variable is not available or there are
            multiple matching cubes in the input `cubelist`.
    """
    # If the cube is in the cubelist simply extract and return it
    newcubelist = cubelist.extract(name)
    if len(newcubelist) == 1:
        return newcubelist[0].copy()

    elif len(newcubelist) > 1:
        raise ValueError('Multiple cubes found matching ' + name +
                         ' not sure which to use')

    # If an equation is present for the requested variable the try to derive the
    # variable from existing cube in the input cubelist
    elif name in available:
        # Calculate all required variables
        args = [calc(var, cubelist) for var in available[name]['required']]
        # Call the function to calculate the requested variable
        cube = available[name]['function'](*args)
        cube.rename(name)
        return cube
    else:
        raise ValueError('Can not get ' + name + ' from cubelist.')


# Define simple maths functions
def _nothing(x):
    return x


def _subtract(x, y):
    return x - y


def _multiply(x, y):
    return x * y


def _divide(x, y):
    return x / y


def _summation(*args):
    return sum(args)


available = {
    'air_potential_temperature': {
        'function': _divide,
        'required': ['air_temperature', 'dimensionless_exner_function']},

    'air_pressure': {
        'function': variable.pressure,
        'required': ['dimensionless_exner_function']},

    'air_temperature': {
        'function': _multiply,
        'required': ['dimensionless_exner_function',
                     'air_potential_temperature']},

    'altitude': {
        'function': variable.height,
        'required': ['air_potential_temperature']},

    'atmosphere_boundary_layer_height': {
        'function': _summation,
        'required': ['surface_altitude',
                     'atmosphere_boundary_layer_thickness']},

    'boundary_layer_type': {
        'function': variable.category_map,
        'required': ['Stable boundary layer indicator',
                     'Stratocumulus over stable boundary layer indicator',
                     'Well-mixed boundary layer indicator',
                     'Decoupled stratocumulus not over cumulus indicator',
                     'Decoupled stratocumulus over cumulus indicator',
                     'Cumulus capped boundary layer indicator',
                     'Shear driven boundary layer indicator']},

    'boundary_layer_height': {
        'function': _summation,
        'required': ['surface_altitude',
                     'atmosphere_boundary_layer_thickness']},

    'brunt_vaisala_frequency': {
        'function': variable.brunt_vaisala,
        'required': ['air_potential_temperature']},

    'brunt_vaisala_frequency_squared': {
        'function': variable.brunt_vaisala_squared,
        'required': ['air_potential_temperature']},

    'air_density': {
        'function': variable.density,
        'required': ['air_pressure', 'air_temperature']},

    'derived_pv': {
        'function': variable.potential_vorticity,
        'required': ['x_wind', 'y_wind', 'upward_air_velocity',
                     'air_potential_temperature', 'air_density']},

    'dimensionless_exner_function': {
        'function': variable.exner,
        'required': ['air_pressure']},

    'eady_growth_rate': {
        'function': variable.Eady,
        'required': ['air_potential_temperature', 'x_wind', 'air_pressure']},

    'equivalent_potential_temperature': {
        'function': variable.theta_e,
        'required': ['air_potential_temperature',
                     'saturated_humidity_mixing_ratio',
                     'air_temperature']},

    'ertel_potential_vorticity': {
        'function': _nothing,
        'required': ['derived_pv']},

    'mass': {
        'function': _multiply,
        'required': ['air_density', 'volume']},

    'mass_fraction_of_cloud': {
        'function': _summation,
        'required': ['mass_fraction_of_cloud_liquid_water_in_air',
                     'mass_fraction_of_cloud_ice_in_air']},

    'mass_fraction_of_water': {
        'function': _summation,
        'required': ['mass_fraction_of_cloud',
                     'specific_humidity']},

    'relative_humidity': {
        'function': _divide,
        'required': ['vapour_pressure', 'saturated_vapour_pressure']},

    'saturated_humidity_mixing_ratio': {
        'function': variable.r_vs,
        'required': ['air_pressure', 'saturated_vapour_pressure']},

    'saturated_vapour_pressure': {
        'function': variable.tetens,
        'required': ['air_temperature']},

    'specific_total_water_content': {
        'function': _summation,
        'required': ['specific_humidity',
                     'mass_fraction_of_cloud_liquid_water_in_air',
                     'mass_fraction_of_cloud_ice_in_air']},

    'surface_altitude': {
        'function': variable.surface_height,
        'required': ['air_potential_temperature']},

    'total_column_water': {
        'function': variable.column_integral,
        'required': ['total_water_content']},

    'total_column_liquid_water': {
        'function': variable.column_integral,
        'required': ['liquid_water_content']},

    'total_column_ice': {
        'function': variable.column_integral,
        'required': ['ice_content']},

    'total_water_content': {
        'function': _multiply,
        'required': ['thickness', 'specific_total_water_content']},

    'liquid_water_content': {
        'function': _multiply,
        'required': ['thickness', 'mass_fraction_of_cloud_liquid_water_in_air']},

    'ice_content': {
        'function': _multiply,
        'required': ['thickness', 'mass_fraction_of_cloud_ice_in_air']},

    'vapour_pressure': {
        'function': variable.vapour_pressure,
        'required': ['air_pressure', 'specific_humidity']},

    'volume': {
        'function': grid.volume,
        'required': ['air_potential_temperature']},

    'wind_speed': {
        'function': variable.wind_speed,
        'required': ['x_wind', 'y_wind', 'upward_air_velocity']},


    # PV tracer functions
    'total_minus_advection_only_pv': {
        'function': _subtract,
        'required': ['ertel_potential_vorticity', 'advection_only_pv']},

    'total_minus_advection_only_theta': {
        'function': _subtract,
        'required': ['air_potential_temperature',
                     'advection_only_theta']},

    'sum_of_physics_pv_tracers': {
        'function': _summation,
        'required': ['short_wave_radiation_pv', 'long_wave_radiation_pv',
                     'microphysics_pv', 'gravity_wave_drag_pv', 'convection_pv',
                     'boundary_layer_pv', 'cloud_rebalancing_pv']},

    'sum_of_physics_theta_tracers': {
        'function': _summation,
        'required': ['short_wave_radiation_theta', 'long_wave_radiation_theta',
                     'microphysics_theta', 'convection_theta',
                     'boundary_layer_theta', 'cloud_rebalancing_theta']},

    'epsilon': {
        'function': _subtract,
        'required': ['total_minus_advection_only_pv',
                     'sum_of_physics_pv_tracers']},

    'residual_pv': {
        'function': _subtract,
        'required': ['epsilon', 'dynamics_tracer_inconsistency']},

    'residual_theta': {
        'function': _subtract,
        'required': ['total_minus_advection_only_theta',
                     'sum_of_physics_theta_tracers']},

    'diabatic_pv': {
        'function': variable.potential_vorticity,
        'required': ['x_wind', 'y_wind', 'upward_air_velocity',
                     'total_minus_advection_only_theta', 'air_density']},

    # Moisture traces
    'total_minus_advection_only_q': {
        'function': _subtract,
        'required': ['specific_humidity', 'advection_only_q']},

    'sum_of_physics_q_tracers': {
        'function': _summation,
        'required': ['short_wave_radiation_q', 'long_wave_radiation_q',
                     'microphysics_q', 'convection_q',
                     'boundary_layer_q']},

    'residual_q': {
        'function': _subtract,
        'required': ['total_minus_advection_only_q',
                     'sum_of_physics_q_tracers']},
}
