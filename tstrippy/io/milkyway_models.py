import yaml
import astropy.constants as const
import logging
from pathlib import Path

# Get the directory of the current script
current_script_directory = Path(__file__).parent

# Construct the path to the data file relative to the current script
path_to_data = current_script_directory / ".." / "data" 
path_to_MWmodels = path_to_data / "MWmodels"
# Resolve the path to make it absolute (and normalize it)
absolute_path_to_MWmodels = path_to_MWmodels.resolve()

# configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

# Supported file handlers
FILE_HANDLERS = {}
def register_handler(extension):
    """Decorator to register a handler for a specific file extension."""
    def decorator(func):
        FILE_HANDLERS[extension] = func
        return func
    return decorator

@register_handler(".yaml")
@register_handler(".yml")
def load_yaml(file_path):
    """Load a YAML file."""
    with open(file_path, 'r') as fp:
        return yaml.safe_load(fp)

def get_model(modelname):
    """Load a Milky Way model based on its name."""
    # Get all files in the MWmodels directory
    models = sorted(path_to_MWmodels.resolve().iterdir())
    matches = sorted(Path(path_to_MWmodels).rglob(f"*{modelname}*"))

    if not matches:
        logger.error(f"Model '{modelname}' not found.")
        logger.info("Available models:")
        for model in models:
            logger.info(f"- {model.name}")
        raise FileNotFoundError(f"Model '{modelname}' not found.")

    if len(matches) > 1:
        logger.warning(f"Multiple matches found for '{modelname}'. Be more specific.")
        for match in matches:
            logger.warning(f"- {match.name}")

    # Get the first match
    model_file = matches[0]
    extension = model_file.suffix

    # Check if the format is supported
    if extension not in FILE_HANDLERS:
        logger.error(f"Unsupported file format: {extension}")
        raise ValueError(f"Unsupported file format: {extension}")

    # Load the model using the appropriate handler
    handler = FILE_HANDLERS[extension]
    return handler(model_file)




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

    
