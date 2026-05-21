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


def milkyway_models(modelname):
    """Return a Milky Way model payload by name."""
    return get_model(modelname)


__all__ = ["milkyway_models", "get_model", "register_handler", "FILE_HANDLERS"]


