from pathlib import Path
import yaml
from astropy import coordinates
from astropy import units as u

# Get the directory of the current script
current_script_directory = Path(__file__).parent

# Construct the path to the data file relative to the current script
path_to_data = current_script_directory / ".." / "data" 

# Resolve the path to make it absolute (and normalize it)
path_to_unit_basis = path_to_data / "unit_basis.yaml"
absolute_path_to_unit_basis = path_to_unit_basis.resolve()


def MWrefframeFerrone2023():
    path_to_frame= path_to_data / "MWrefframeFerrone2023.yaml"
    absolute_path_to_frame = path_to_frame.resolve()
    with open(absolute_path_to_frame, 'r') as frame:
        try:
            frame_data = yaml.safe_load(frame)
        except yaml.YAMLError as exc:
            print(exc)
        galcen_distance=frame_data['value']['galcen_distance']*u.Unit(frame_data['unit']['galcen_distance'])
        z_sun   =   frame_data['value']['z_sun']    *   u.Unit(frame_data['unit']['z_sun'])
        vSun    =   frame_data['value']['vSun']     *   u.Unit(frame_data['unit']['vSun'])
        vLSR    =   frame_data['value']['vLSR']     *   u.Unit(frame_data['unit']['vLSR'])
    return coordinates.Galactocentric(galcen_distance = galcen_distance, galcen_v_sun=vLSR+vSun, z_sun=z_sun)
