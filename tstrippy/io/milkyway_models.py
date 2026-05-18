import yaml
import astropy.constants as const
from pathlib import Path

# Get the directory of the current script
current_script_directory = Path(__file__).parent

# Construct the path to the data file relative to the current script
path_to_data = current_script_directory / ".." / "data" 

# Resolve the path to make it absolute (and normalize it)
path_to_unit_basis = path_to_data / "unit_basis.yaml"
absolute_path_to_unit_basis = path_to_unit_basis.resolve()


#### MAKE THE G constant a global variable 
with open(absolute_path_to_unit_basis, 'r') as basis:
    try:
        unitbasis = yaml.safe_load(basis)
    except yaml.YAMLError as exc:
        print(exc)
G=const.G.to(unitbasis['G']).value

            
def pouliasis2017pii():
    path_to_potential = path_to_data / "pouliasis2017pii.yaml"
    absolute_path_to_potential = path_to_potential.resolve()
    with open(absolute_path_to_potential, 'r') as potential:
        try:
            potential_parameters = yaml.safe_load(potential)
        except yaml.YAMLError as exc:
            print(exc)

    components = []
    for component in potential_parameters['components']:
        parameters = component.get('parameters', [])
        if isinstance(parameters, dict):
            # Convert ordered mapping to the runtime-friendly list format.
            parameter_names = component.get('parameter_names')
            if parameter_names is None:
                parameter_names = list(parameters.keys())
            parameters = [parameters[name] for name in parameter_names]
        components.append({
            'name': component['name'],
            'parameters': parameters,
            'parameter_names': component.get('parameter_names', []),
            'parameter_units': component.get('parameter_units', []),
        })

    return components


def pouliasis2017pii_flat():
    """Return the historical flat parameter list used by older call sites."""
    components = pouliasis2017pii()
    return [value for component in components for value in component['parameters']]

    
