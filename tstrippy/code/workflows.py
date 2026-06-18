import numpy as np 


MAX_ORBIT_STEPS = int(1e6) # sorry, beyond this is excessive 
# the recipees
def compute_backward_and_forward_host_orbit(simulator, schemeName, dt, NstepsBack, NstepsForward, initialconditions, mwmodel, t0=0, trim_orbits=1 ):
    _require_methods(
        simulator,
        (
            "clear",
            "setinitialconditions",
            "setscheme",
            "add_component",
            "setbackwardorbit",
            "finalize",
            "run",
        ),
        "simulator",
    )
    if not isinstance(initialconditions, (list, tuple)) or len(initialconditions) != 6:
        raise ValueError("initialconditions must be a sequence of length 6")
    _validate_components(mwmodel)

    # compute the orbit backward 
    # BUILD THE BACK SCHEME
    SCHEME_BACKWARD = schemeName, [t0,dt,NstepsBack]
    SCHEME_FORWARD = schemeName, [t0, dt, NstepsForward]
    simulator.clear()
    simulator.setinitialconditions(*initialconditions)
    simulator.setscheme(*SCHEME_BACKWARD)
    for comp in mwmodel["components"]:
        simulator.add_component(comp["name"], comp["parameters"])
    simulator.setbackwardorbit()
    simulator.trim_orbits(trim_orbits)
    simulator.finalize()
    simulator.run()
    orbits=simulator.orbits.copy()
    timestamps_backward=simulator.orbits_timestamps.copy()
    orbits_backward = orbits[:,:,0] # only one particle. go down a dimension
    # now do it forward
    simulator.clear()
    simulator.setinitialconditions(*initialconditions)
    simulator.setscheme(*SCHEME_FORWARD)
    for comp in mwmodel["components"]:
        simulator.add_component(comp["name"], comp["parameters"])
    simulator.trim_orbits(trim_orbits)
    simulator.finalize()
    simulator.run()
    orbits=simulator.orbits.copy()
    timestamps_forward=simulator.orbits_timestamps.copy()
    orbits_forward = orbits[:,:,0] # only one particle. go down a dimension
    simulator.clear()
    # concatenate the two
    timestamps = np.concatenate((timestamps_backward[::-1], timestamps_forward[1:]))
    orbit = np.concatenate( (orbits_backward[::-1, :], orbits_forward[1:, :]), axis=0)
    return timestamps, orbit.T




def compute_host_orbit(simulator, initialconditions, mwmodel, scheme, trim_orbits=100):
    _require_methods(
        simulator,
        (
            "clear",
            "setinitialconditions",
            "setscheme",
            "add_component",
            "setbackwardorbit",
            "finalize",
            "run",
        ),
        "simulator",
    )
    if not isinstance(initialconditions, (list, tuple)) or len(initialconditions) != 6:
        raise ValueError("initialconditions must be a sequence of length 6")
    if not isinstance(scheme, (list, tuple)) or len(scheme) != 2:
        raise ValueError("scheme must be a 2-item sequence: [name, [t0, dt, nsteps]]")
    _validate_components(mwmodel)

    simulator.clear()
    simulator.setinitialconditions(*initialconditions)
    simulator.setscheme(*scheme)
    for comp in mwmodel["components"]:
        simulator.add_component(comp["name"], comp["parameters"])
    simulator.setbackwardorbit()
    simulator.trim_orbits(trim_orbits)
    simulator.finalize()
    simulator.run()
    orbits=simulator.orbits.copy()
    timestamps=simulator.orbits_timestamps.copy()
    simulator.clear()
    orbits = orbits[:,:,0]
    return timestamps, orbits.T


def generate_vanilla_stream(simulator, mwmodel, scheme, initialconditions, hostcluster_kinematics, hostcluster_structure, hostcluster_structure_parameter_tables=None, trim_orbits=100):
    """with the option to update the parameters """
    _require_methods(simulator,
                     ("clear",
                      "setinitialconditions",
                      "setscheme",
                      "add_component",
                      "setbackwardorbit",
                      "finalize",
                      "run",
                      "configure_hostcluster_structure",
                      "configure_hostcluster_kinematics",
                      "configure_hostcluster_structure_parameter_table"),
                     "simulator",) 
    _validate_components(mwmodel)

    simulator.clear()
    simulator.setinitialconditions(*initialconditions)
    simulator.setscheme(*scheme)
    for comp in mwmodel['components']:
        simulator.add_component(comp['name'],comp['parameters'])
    simulator.configure_hostcluster_kinematics(*hostcluster_kinematics)
    simulator.configure_hostcluster_structure(*hostcluster_structure)
    _configure_hostcluster_structure_parameter_tables(
        simulator,
        hostcluster_structure_parameter_tables,
    )


    simulator.trim_orbits(trim_orbits)
    simulator.finalize()
    simulator.run()
    return simulator


# helper functions
def _downsample_time_series(timestamps, values, max_steps):
    """Downsample along time axis 0 while preserving endpoints."""
    if not isinstance(max_steps, (int, np.integer)) or int(max_steps) < 2:
        raise ValueError("max_steps must be an integer >= 2")
    timestamps = np.asarray(timestamps)
    if timestamps.ndim != 1:
        raise ValueError("timestamps must be a 1D array")

    nsteps = timestamps.shape[0]
    if nsteps <= max_steps:
        return timestamps, values
    indices = np.linspace(0, nsteps - 1, max_steps, dtype=np.int64)
    return timestamps[indices], values[indices, ...]


# checks
def _normalize_hostcluster_structure_parameter_tables(hostcluster_structure_parameter_tables):
    """Normalize None/single-table/list-of-tables into a list of tables."""
    if hostcluster_structure_parameter_tables is None:
        return []

    if _looks_like_single_parameter_table(hostcluster_structure_parameter_tables):
        return [hostcluster_structure_parameter_tables]

    if isinstance(hostcluster_structure_parameter_tables, (list, tuple)):
        return list(hostcluster_structure_parameter_tables)

    print("hostcluster_structure_parameter_tables value:", hostcluster_structure_parameter_tables)
    raise TypeError(
        "hostcluster_structure_parameter_tables must be None, a single parameter table, "
        "or a list/tuple of parameter tables; "
        f"got {type(hostcluster_structure_parameter_tables).__name__} "
        f"with value {hostcluster_structure_parameter_tables!r}"
    )


def _looks_like_single_parameter_table(value):
    """A parameter table is [index, timestamps, parameter_values]."""
    return isinstance(value, (list, tuple)) and len(value) == 3


def _validate_parameter_table(table, idx=None):
    """Validate one parameter table and return normalized values."""
    label = f"table[{idx}]" if idx is not None else "table"
    if not isinstance(table, (list, tuple)) or len(table) != 3:
        print(f"{label} value:", table)
        raise TypeError(
            f"{label} must be a 3-item list/tuple: [index, timestamps, parameter_values]; "
            f"got {type(table).__name__} with value {table!r}"
        )

    index, timestamps, parameter_values = table

    if not isinstance(index, (int, np.integer)):
        print(f"{label} index value:", index)
        raise TypeError(
            f"{label}[0] (index) must be an integer; "
            f"got {type(index).__name__} with value {index!r}"
        )

    timestamps_arr = np.asarray(timestamps)
    parameter_values_arr = np.asarray(parameter_values)

    if timestamps_arr.ndim != 1:
        print(f"{label} timestamps value:", timestamps)
        raise ValueError(
            f"{label}[1] (timestamps) must be 1D; "
            f"got array with shape {timestamps_arr.shape}"
        )

    if parameter_values_arr.ndim != 1:
        print(f"{label} parameter_values value:", parameter_values)
        raise ValueError(
            f"{label}[2] (parameter_values) must be 1D; "
            f"got array with shape {parameter_values_arr.shape}"
        )

    if timestamps_arr.shape[0] != parameter_values_arr.shape[0]:
        print(f"{label} timestamps value:", timestamps)
        print(f"{label} parameter_values value:", parameter_values)
        raise ValueError(
            f"{label} timestamps and parameter_values must have the same length; "
            f"got {timestamps_arr.shape[0]} and {parameter_values_arr.shape[0]}"
        )

    return int(index), timestamps_arr, parameter_values_arr


def _configure_hostcluster_structure_parameter_tables(simulator, hostcluster_structure_parameter_tables):
    """Normalize, validate, and upload host cluster structure parameter tables."""
    tables = _normalize_hostcluster_structure_parameter_tables(hostcluster_structure_parameter_tables)
    for i, table in enumerate(tables):
        index, timestamps, parameter_values = _validate_parameter_table(table, idx=i)
        try:
            simulator.configure_hostcluster_structure_parameter_table(
                index,
                timestamps,
                parameter_values,
            )
        except Exception as exc:
            print(f"failed to upload hostcluster_structure_parameter_table at table[{i}]")
            print("table value:", table)
            raise RuntimeError(
                f"Failed to configure hostcluster_structure_parameter_table at table[{i}] "
                f"(index={index})"
            ) from exc

def _require_methods(obj, methods, obj_name):
    missing = [name for name in methods if not hasattr(obj, name)]
    if missing:
        raise TypeError(f"{obj_name} is missing required methods: {missing}")


def _validate_components(mwmodel):
    if not isinstance(mwmodel, dict) or "components" not in mwmodel:
        raise TypeError("mwmodel must be a dictionary with a 'components' key")
    if not isinstance(mwmodel["components"], (list, tuple)):
        raise TypeError("mwmodel['components'] must be a list or tuple")

    for i, comp in enumerate(mwmodel["components"]):
        if not isinstance(comp, dict) or "name" not in comp or "parameters" not in comp:
            raise TypeError(
                f"mwmodel['components'][{i}] must be a dict with 'name' and 'parameters'"
            )

