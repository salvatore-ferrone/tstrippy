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
            
    myparams= []
    myparams.append(G)
    myparams.append(potential_parameters['components'][0]["parameters"]['M'])
    myparams.append(potential_parameters['components'][0]["parameters"]['a'])
    myparams.append(potential_parameters['components'][0]["parameters"]['exp'])
    myparams.append(potential_parameters['components'][0]["parameters"]['cutoffradius'])
    myparams.append(potential_parameters['components'][1]["parameters"]['M'])
    myparams.append(potential_parameters['components'][1]["parameters"]['a'])
    myparams.append(potential_parameters['components'][1]["parameters"]['b'])
    myparams.append(potential_parameters['components'][2]["parameters"]['M'])
    myparams.append(potential_parameters['components'][2]["parameters"]['a'])
    myparams.append(potential_parameters['components'][2]["parameters"]['b'])
    return myparams

    
