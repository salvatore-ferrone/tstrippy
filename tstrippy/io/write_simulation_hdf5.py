"""
tstrippy/io/write_simulation_hdf5.py

Convert a completed tstrippy simulation to a standardized HDF5 file.

HDF5 schema v1.0
----------------
/meta                       attrs: schema_version, created_utc, code
/config                     attrs: nparticles
    timestamps              float64 (Nt,)   — full integration timestamps
    backwardorbit           bool
    /initialconditions
        x, y, z, vx, vy, vz float64 (Np,)
    /scheme
        method              str
        parameters          float64 (Nparams,)
/snapshots                  attrs: phase_labels, shape_convention, nsaved
    state                   float32 (Nsaved, 6, Np)  — [time, phase, particle]
    time                    float64 (Nsaved,)
/orbits  (optional)         attrs: phase_labels, shape_convention
    state                   float32 (No, 6, Nt)      — [object, phase, time]
    time                    float64 (Nt,)
    object_names            str     (No,)  optional
"""

import datetime
import glob
import os

import h5py
import numpy as np

SCHEMA_VERSION = "1.0"
_PHASE_LABELS = ["x", "y", "z", "vx", "vy", "vz"]


# ---------------------------------------------------------------------------
# Public entry point
# ---------------------------------------------------------------------------

def write_simulation_hdf5(simulator, filename,
                          objectnames=None, delete_temp_binaries=False):
    """Write a completed tstrippy simulation to a standardized HDF5 file.

    Parameters
    ----------
    simulator : tstrippy.simulator module
        The simulator module after run() has completed successfully.
    filename : str
        Path to the output HDF5 file.
    objectnames : list of str, optional
        Display names for orbit objects.  Length must match number of tracked
        orbit objects.
    delete_temp_binaries : bool
        If True, delete temporary snapshot/orbit binary files after ingestion.
    """
    if not bool(simulator.run_success):
        raise RuntimeError(
            "simulator.run_success is False. "
            "Complete a successful run() before writing output."
        )

    out_dir = os.path.dirname(os.path.abspath(filename))
    os.makedirs(out_dir, exist_ok=True)

    with h5py.File(filename, "w") as f:
        _write_meta_group(f, simulator)
        _write_config_group(f, simulator)
        _write_snapshots_group(f, simulator, delete_temp_binaries)
        _write_orbits_group(f, simulator, objectnames, delete_temp_binaries)


# ---------------------------------------------------------------------------
# Group writers
# ---------------------------------------------------------------------------

def _write_meta_group(f, simulator):
    import platform
    grp = f.create_group("meta")
    grp.attrs["schema_version"] = SCHEMA_VERSION
    grp.attrs["created_utc"] = datetime.datetime.utcnow().isoformat()
    grp.attrs["code"] = "tstrippy"

    # Timing data
    timing_grp = grp.create_group("timing")
    timing_grp.create_dataset("run_seconds", data=float(simulator.timer_run_seconds))
    timing_grp.create_dataset("finalize_seconds", data=float(simulator.timer_finalize_seconds))
    timing_grp.create_dataset("scheme_seconds", data=float(simulator.timer_scheme_seconds))
    timing_grp.create_dataset("write_snapshots_seconds", data=float(simulator.timer_write_snapshots_seconds))
    timing_grp.create_dataset("write_orbits_seconds", data=float(simulator.timer_write_orbits_seconds))

    # Machine info
    machine_grp = grp.create_group("machine")
    machine_grp.create_dataset("architecture", data=platform.machine().encode('ascii'))
    machine_grp.create_dataset("processor", data=platform.processor().encode('ascii'))
    machine_grp.create_dataset("system", data=platform.system().encode('ascii'))
    machine_grp.create_dataset("release", data=platform.release().encode('ascii'))
    machine_grp.create_dataset("hostname", data=platform.node().encode('ascii'))


def _write_config_group(f, simulator):
    grp = f.create_group("config")
    grp.attrs["nparticles"] = int(simulator.nparticles)
    grp.attrs["backwardorbit"] = bool(getattr(simulator, "backward_orbit_enabled", False))

    timestamps = np.asarray(simulator.timestamps, dtype=np.float64)
    grp.create_dataset("timestamps", data=timestamps)

    ic_grp = grp.create_group("initialconditions")
    x0 = _safe_array(simulator, "x_initial")
    y0 = _safe_array(simulator, "y_initial")
    z0 = _safe_array(simulator, "z_initial")
    vx0 = _safe_array(simulator, "vx_initial")
    vy0 = _safe_array(simulator, "vy_initial")
    vz0 = _safe_array(simulator, "vz_initial")

    # Fallback: first in-memory orbit snapshot if dedicated initial arrays are unavailable.
    if x0 is None or y0 is None or z0 is None or vx0 is None or vy0 is None or vz0 is None:
        orbits = _safe_orbits(simulator)
        if orbits is not None and orbits.ndim >= 3 and orbits.shape[0] > 0:
            x0 = orbits[0, 0, :]
            y0 = orbits[0, 1, :]
            z0 = orbits[0, 2, :]
            vx0 = orbits[0, 3, :]
            vy0 = orbits[0, 4, :]
            vz0 = orbits[0, 5, :]

    if x0 is not None and y0 is not None and z0 is not None and vx0 is not None and vy0 is not None and vz0 is not None:
        ic_grp.create_dataset("x", data=np.asarray(x0, dtype=np.float64))
        ic_grp.create_dataset("y", data=np.asarray(y0, dtype=np.float64))
        ic_grp.create_dataset("z", data=np.asarray(z0, dtype=np.float64))
        ic_grp.create_dataset("vx", data=np.asarray(vx0, dtype=np.float64))
        ic_grp.create_dataset("vy", data=np.asarray(vy0, dtype=np.float64))
        ic_grp.create_dataset("vz", data=np.asarray(vz0, dtype=np.float64))

    scheme_grp = grp.create_group("scheme")
    method = _fortran_str(getattr(simulator, "scheme_method", ""))
    if method == "":
        method = "unknown"
    scheme_grp.create_dataset("method", data=method.encode("ascii", errors="replace"))

    nparams = int(getattr(simulator, "scheme_nparams", 0))
    params_raw = np.asarray(getattr(simulator, "scheme_parameters", np.array([], dtype=np.float64)), dtype=np.float64)
    if nparams > 0 and params_raw.ndim == 1 and params_raw.size >= nparams:
        params = params_raw[:nparams]
    else:
        params = params_raw
    scheme_grp.create_dataset("parameters", data=params.astype(np.float64))


# def _write_gravity_group(f, gravity):
#     grp = f.create_group("gravity")
#     grp.attrs["G"] = float(gravity.gravity_g)
#     ncomp = int(gravity.gravity_ncomp)
#     grp.attrs["ncomponents"] = ncomp

#     for i in range(1, ncomp + 1):
#         # component_model_names is a (16,) char array; index i-1 (0-based)
#         raw_name = gravity.component_model_names[i - 1]
#         if hasattr(raw_name, 'decode'):
#             name = raw_name.decode("ascii").strip()
#         else:
#             name = str(raw_name).strip()
#         # Derive nparams from gravity_params: count non-zero entries in the column
#         col = np.asarray(gravity.gravity_params[:, i - 1], dtype=np.float64)
#         # Use getcomponentnparams if available, otherwise fall back to GRAVITY_MAX_PARAMS
#         try:
#             nparams = int(gravity.getcomponentnparams(i))
#         except AttributeError:
#             # count non-zero trailing entries
#             nparams = int(gravity.gravity_max_params)
#         params = col[:nparams]

#         cgrp = grp.create_group(f"component_{i:02d}")
#         cgrp.attrs["model_name"] = name
#         cgrp.attrs["component_index"] = i
#         cgrp.create_dataset("params", data=params)


def _write_snapshots_group(f, simulator, delete_temp_binaries=False):
    grp = f.create_group("snapshots")

    snapshots = []
    times = []

    if bool(simulator.did_write_snapshots):
        directory = _fortran_str(simulator.directory_snapshots)
        basename  = _fortran_str(simulator.basename_snapshots)
        pattern   = os.path.join(directory, f"{basename}_*.bin")
        files     = sorted(glob.glob(pattern))

        for fpath in files:
            t, x, y, z, vx, vy, vz = read_snapshot_binary(fpath)
            snapshots.append(np.stack([x, y, z, vx, vy, vz], axis=0))  # (6, Np)
            times.append(t)

        if delete_temp_binaries:
            for fpath in files:
                os.remove(fpath)

    # Always ensure the final state is the last entry
    final_time = float(np.asarray(simulator.timestamps)[-1])
    final_phase = np.stack([
        np.asarray(simulator.x,  dtype=np.float32),
        np.asarray(simulator.y,  dtype=np.float32),
        np.asarray(simulator.z,  dtype=np.float32),
        np.asarray(simulator.vx, dtype=np.float32),
        np.asarray(simulator.vy, dtype=np.float32),
        np.asarray(simulator.vz, dtype=np.float32),
    ], axis=0)  # (6, Np)

    if len(times) == 0 or not np.isclose(times[-1], final_time, rtol=1e-6):
        snapshots.append(final_phase)
        times.append(final_time)

    # Stack → (Nsaved, 6, Np) = [time, phase, particle]
    state_array = np.stack(snapshots, axis=0)
    grp.create_dataset("state", data=state_array, compression="gzip", compression_opts=4)
    grp.create_dataset("time",  data=np.array(times, dtype=np.float64))
    grp.attrs["phase_labels"]      = _PHASE_LABELS
    grp.attrs["shape_convention"]  = "[time, phase, particle]"
    grp.attrs["nsaved"]            = len(times)


def _write_orbits_group(f, simulator, objectnames=None, delete_temp_binaries=False):
    orbits      = _safe_orbits(simulator)
    orb_times   = _safe_orbits_timestamps(simulator)
    has_memory  = (orbits is not None and orbits.ndim >= 3 and 
                   orbits.shape[2] > 0)
    has_files   = bool(simulator.did_write_orbits)

    if not has_memory and not has_files:
        return

    grp = f.create_group("orbits")

    if has_files:
        directory = _fortran_str(simulator.directory_orbits)
        basename  = _fortran_str(simulator.basename_orbits)
        pattern   = os.path.join(directory, f"{basename}_particle_*.bin")
        files     = sorted(glob.glob(pattern))

        all_states = []
        common_times = None

        for fpath in files:
            times, states = read_orbit_binary(fpath)  # states: (Nt, 6)
            all_states.append(states)
            if common_times is None:
                common_times = times

        if delete_temp_binaries:
            for fpath in files:
                os.remove(fpath)

        if all_states:
            # (No, Nt, 6) → transpose → (No, 6, Nt) = [object, phase, time]
            state_array = np.transpose(np.stack(all_states, axis=0), (0, 2, 1))
            grp.create_dataset("state", data=state_array.astype(np.float32),
                               compression="gzip", compression_opts=4)
            grp.create_dataset("time", data=common_times.astype(np.float64))

    elif has_memory:
        # orbits shape: (Nsaved, 6, Nobjects)
        # transpose → (Nobjects, 6, Nsaved) = [object, phase, time]
        state_array = np.transpose(orbits, (2, 1, 0))
        grp.create_dataset("state", data=state_array.astype(np.float32),
                           compression="gzip", compression_opts=4)
        grp.create_dataset("time", data=orb_times.astype(np.float64))

    grp.attrs["phase_labels"]     = _PHASE_LABELS
    grp.attrs["shape_convention"] = "[object, phase, time]"

    if objectnames is not None:
        nobjects = grp["state"].shape[0]
        if len(objectnames) != nobjects:
            raise ValueError(
                f"objectnames length ({len(objectnames)}) != "
                f"number of orbit objects ({nobjects})"
            )
        dt = h5py.string_dtype(encoding="ascii")
        grp.create_dataset("object_names",
                           data=np.array(objectnames, dtype=object), dtype=dt)


# ---------------------------------------------------------------------------
# Binary file readers (numpy only, no scipy)
# ---------------------------------------------------------------------------

def read_snapshot_binary(filename):
    """Read one Fortran sequential-unformatted snapshot binary file.

    Written by simulator.write_snapshot_file.  Layout::

        record 1 : int32 Nparticles, float32 time
        record 2 : float32[Nparticles]  x
        record 3 : float32[Nparticles]  y
        record 4 : float32[Nparticles]  z
        record 5 : float32[Nparticles]  vx
        record 6 : float32[Nparticles]  vy
        record 7 : float32[Nparticles]  vz

    Each Fortran record is wrapped with 4-byte integer markers (record length
    in bytes).

    Returns
    -------
    time : float
    x, y, z, vx, vy, vz : np.ndarray float32, shape (Nparticles,)
    """
    def _read_record(fh, dtype):
        marker = int(np.frombuffer(fh.read(4), dtype=np.int32)[0])
        data   = np.frombuffer(fh.read(marker), dtype=dtype).copy()
        fh.read(4)  # trailing marker
        return data

    with open(filename, "rb") as fh:
        # Record 1: Nparticles (int32) + time (float32) packed in one record
        marker = int(np.frombuffer(fh.read(4), dtype=np.int32)[0])
        nparticles = int(np.frombuffer(fh.read(4), dtype=np.int32)[0])
        time       = float(np.frombuffer(fh.read(4), dtype=np.float32)[0])
        fh.read(4)  # trailing marker

        x  = _read_record(fh, np.float32)
        y  = _read_record(fh, np.float32)
        z  = _read_record(fh, np.float32)
        vx = _read_record(fh, np.float32)
        vy = _read_record(fh, np.float32)
        vz = _read_record(fh, np.float32)

    return time, x, y, z, vx, vy, vz


def read_orbit_binary(filename):
    """Read one Fortran sequential-unformatted per-particle orbit file.

    Written by simulator.write_orbit_records.  Each record contains seven
    float32 values: time, x, y, z, vx, vy, vz (28 bytes of payload).

    Returns
    -------
    times  : np.ndarray float64, shape (Nt,)
    states : np.ndarray float32, shape (Nt, 6)  — columns: x y z vx vy vz
    """
    PAYLOAD_BYTES = 7 * 4  # 7 × float32
    times  = []
    states = []

    with open(filename, "rb") as fh:
        while True:
            marker_bytes = fh.read(4)
            if len(marker_bytes) < 4:
                break
            marker = int(np.frombuffer(marker_bytes, dtype=np.int32)[0])
            if marker != PAYLOAD_BYTES:
                break
            data = np.frombuffer(fh.read(PAYLOAD_BYTES), dtype=np.float32).copy()
            fh.read(4)  # trailing marker
            times.append(data[0])
            states.append(data[1:])  # 6 values: x y z vx vy vz

    return np.array(times, dtype=np.float64), np.array(states, dtype=np.float32)


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

def _fortran_str(fortran_char_array):
    """Convert a f2py CHARACTER array to a stripped Python string."""
    try:
        arr = np.asarray(fortran_char_array)
        if arr.dtype.kind == "S":
            if arr.ndim == 0:
                return arr.item().decode("ascii", errors="replace").strip()
            return b"".join(arr.tolist()).decode("ascii", errors="replace").strip()
        if arr.dtype.kind == "U":
            if arr.ndim == 0:
                return arr.item().strip()
            return "".join(arr.tolist()).strip()
        raw = bytes(arr).decode("ascii", errors="replace")
        return raw.strip()
    except Exception:
        return str(fortran_char_array).strip()


def _safe_array(simulator, name):
    """Return simulator.<name> as ndarray, or None when unavailable."""
    try:
        arr = np.asarray(getattr(simulator, name))
        if arr.size == 0:
            return None
        return arr
    except Exception:
        return None


def _safe_orbits(simulator):
    """Return simulator.orbits as ndarray, or None if unallocated / empty."""
    try:
        arr = np.asarray(simulator.orbits)
        if arr.size == 0:
            return None
        return arr
    except Exception:
        return None


def _safe_orbits_timestamps(simulator):
    """Return simulator.orbits_timestamps as ndarray, or None if unallocated."""
    try:
        arr = np.asarray(simulator.orbits_timestamps)
        if arr.size == 0:
            return None
        return arr
    except Exception:
        return None
