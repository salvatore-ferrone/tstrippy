from pathlib import Path
import yaml

# Get the directory of the current script
current_script_directory = Path(__file__).parent

# Construct the path to the data file relative to the current script
path_to_data = current_script_directory / ".." / "data" / "MWreferenceframes"


def load_yaml(file_path):
    """Load a YAML file."""
    with open(file_path, "r", encoding="utf-8") as fp:
        return yaml.safe_load(fp)


def available_reference_frames():
    """Return available reference frame YAML file stems."""
    names = []
    for candidate in sorted(path_to_data.resolve().glob("*.y*ml")):
        names.append(candidate.stem)
    return names


def _resolve_reference_frame_path(referencename):
    if not referencename:
        raise ValueError("referencename must be a non-empty string")

    matches = sorted(path_to_data.resolve().rglob(f"*{referencename}*.y*ml"))
    if not matches:
        available = available_reference_frames()
        available_msg = ", ".join(available) if available else "<none>"
        raise FileNotFoundError(
            f"Reference frame '{referencename}' not found. "
            f"Available reference frames: {available_msg}"
        )

    if len(matches) > 1:
        found = ", ".join(m.name for m in matches)
        raise ValueError(
            f"Multiple matches found for '{referencename}': {found}. "
            "Please be more specific."
        )

    return matches[0]


def _validate_galactocentric_schema(frame_data, source="<unknown>"):
    if not isinstance(frame_data, dict):
        raise ValueError(f"Reference frame '{source}' must be a mapping")

    if "value" not in frame_data or "unit" not in frame_data:
        raise ValueError(
            f"Reference frame '{source}' must contain top-level 'value' and 'unit' maps"
        )

    values = frame_data["value"]
    units = frame_data["unit"]
    if not isinstance(values, dict) or not isinstance(units, dict):
        raise ValueError(
            f"Reference frame '{source}' entries 'value' and 'unit' must be mappings"
        )

    required = ["ra", "dec", "vLSR", "vSun", "galcen_distance", "z_sun"]
    for key in required:
        if key not in values:
            raise ValueError(f"Reference frame '{source}' is missing value['{key}']")
        if key not in units:
            raise ValueError(f"Reference frame '{source}' is missing unit['{key}']")


def _build_galactocentric_from_yaml(frame_data):
    """Build an astropy Galactocentric frame from a standardized YAML payload."""
    try:
        from astropy import coordinates
        from astropy import units as u
    except Exception as exc:
        raise ImportError(
            "Astropy is required to build Galactocentric reference frames. "
            "Install a NumPy/Astropy compatible pair in this environment."
        ) from exc

    values = frame_data["value"]
    unit_map = frame_data["unit"]

    galcen_distance = values["galcen_distance"] * u.Unit(unit_map["galcen_distance"])
    ra = values["ra"] * u.Unit(unit_map["ra"])
    dec = values["dec"] * u.Unit(unit_map["dec"])
    z_sun = values["z_sun"] * u.Unit(unit_map["z_sun"])
    v_sun = values["vSun"] * u.Unit(unit_map["vSun"])
    v_lsr = values["vLSR"] * u.Unit(unit_map["vLSR"])
    galcen_coord = coordinates.ICRS(ra=ra.to(u.deg), dec=dec.to(u.deg))

    return coordinates.Galactocentric(
        galcen_coord=galcen_coord,
        galcen_distance=galcen_distance,
        galcen_v_sun=v_lsr + v_sun,
        z_sun=z_sun,
    )


def get_reference_frame(referencename):
    """Load a Galactocentric reference frame by file name or partial name."""
    frame_path = _resolve_reference_frame_path(referencename)
    frame_data = load_yaml(frame_path)
    _validate_galactocentric_schema(frame_data, source=frame_path.name)
    return _build_galactocentric_from_yaml(frame_data)


def reference_frames(referencename):
    """Return a reference frame payload by name."""
    return get_reference_frame(referencename)


def MWrefframeFerrone2023():
    """Backward-compatible convenience accessor for Ferrone 2023 frame."""
    return get_reference_frame("MWrefframeFerrone2023")


__all__ = [
    "reference_frames",
    "get_reference_frame",
    "available_reference_frames",
    "load_yaml",
    "MWrefframeFerrone2023",
]
