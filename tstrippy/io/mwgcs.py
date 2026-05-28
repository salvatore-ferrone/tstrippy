from pathlib import Path
import re
import numpy as np


_DATA_DIR = (Path(__file__).parent / ".." / "data").resolve()
_KIN_PATH = _DATA_DIR / "mwgcs-kinematics"
_STR_PATH = _DATA_DIR / "mwgcs-structural"

_CACHE = None

baumgardt_to_astropy_headers = {
    "ra": "ra",
    "dec": "dec",
    "rsun": "distance",
    "mualpha": "pm_ra_cosdec",
    "mu_delta": "pm_dec",
    "rv": "radial_velocity"
}


_baumgardt_structural_columns_with_errors = {
    "Mass": "DM",
    "v": "Delta_V",
    "M/L_V": "DM/L",
    "MF": "Delta_MF"
}



def _normalize_name(name):
    return re.sub(r"[^A-Za-z0-9]", "", str(name)).lower()


def _normalize_column_name(name):
    cleaned = str(name).strip()
    cleaned = cleaned.replace("<", "").replace(">", "")
    cleaned = re.sub(r"[^A-Za-z0-9]+", "_", cleaned).strip("_").lower()
    return cleaned or "col"


def _safe_float(values):
    out = np.full(len(values), np.nan, dtype=float)
    for i, v in enumerate(values):
        try:
            out[i] = float(v)
        except (TypeError, ValueError):
            out[i] = np.nan
    return out


def _column_array(values):
    as_float = _safe_float(values)
    if np.all(np.isnan(as_float)):
        return np.array([str(v) for v in values], dtype=str)

    # Keep true name-like columns as strings.
    as_str = [str(v) for v in values]
    if any(not _is_number_token(v) for v in as_str):
        return np.array(as_str, dtype=str)
    return as_float


def _is_number_token(text):
    try:
        float(text)
        return True
    except (TypeError, ValueError):
        return False


def _parse_units_line(line):
    text = line.lstrip("#").strip()
    if not text:
        return []

    out = []
    for token in text.split():
        # Expand glued bracket tokens like "[mas/yr][mas/yr][mas/yr]".
        grouped = re.findall(r"\[[^\]]*\]", token)
        if grouped and "".join(grouped) == token:
            out.extend(grouped)
        else:
            out.append(token)

    # Normalize to plain unit labels without wrappers/quotes.
    cleaned = []
    for unit in out:
        s = str(unit).strip().strip("\"'")
        if s.startswith("[") and s.endswith("]"):
            s = s[1:-1]
        cleaned.append(s)
    return cleaned


def _parse_ascii_catalog(path):
    lines = Path(path).read_text(encoding="utf-8").splitlines()
    if len(lines) < 3:
        raise ValueError(f"Catalog file has insufficient header lines: {path}")

    raw_headers = lines[0].lstrip("#").strip().split()
    raw_units = _parse_units_line(lines[1])

    if len(raw_units) < len(raw_headers):
        raw_units = raw_units + [""] * (len(raw_headers) - len(raw_units))
    else:
        raw_units = raw_units[: len(raw_headers)]

    rows = []
    for raw in lines[3:]:
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        tokens = line.split()
        if len(tokens) < len(raw_headers):
            tokens = tokens + ["nan"] * (len(raw_headers) - len(tokens))
        else:
            tokens = tokens[: len(raw_headers)]
        rows.append(tokens)

    if not rows:
        raise ValueError(f"No data rows found in {path}")

    canon_headers = []
    seen = {}
    for h in raw_headers:
        base = _normalize_column_name(h)
        seen[base] = seen.get(base, 0) + 1
        if seen[base] > 1:
            canon_headers.append(f"{base}_{seen[base]}")
        else:
            canon_headers.append(base)

    table = np.array(rows, dtype=object)
    columns = {}
    for j, key in enumerate(canon_headers):
        columns[key] = _column_array(table[:, j])

    raw_to_canon = {}
    for raw, canon in zip(raw_headers, canon_headers):
        raw_to_canon.setdefault(raw, []).append(canon)

    canon_to_raw = {canon: raw for raw, canon in zip(raw_headers, canon_headers)}

    return {
        "path": str(path),
        "raw_headers": raw_headers,
        "canon_headers": canon_headers,
        "raw_units": raw_units,
        "columns": columns,
        "raw_to_canon": raw_to_canon,
        "canon_to_raw": canon_to_raw,
    }


def _build_cache():
    kin = _parse_ascii_catalog(_KIN_PATH)
    st = _parse_ascii_catalog(_STR_PATH)

    kin_names = kin["columns"]["cluster"].astype(str)
    st_names = st["columns"]["cluster"].astype(str)

    kin_index = {_normalize_name(n): i for i, n in enumerate(kin_names)}
    st_index = {_normalize_name(n): i for i, n in enumerate(st_names)}

    return {
        "kinematics": kin,
        "structural": st,
        "kin_index": kin_index,
        "st_index": st_index,
    }


def _get_cache():
    global _CACHE
    if _CACHE is None:
        _CACHE = _build_cache()
    return _CACHE


def _get_table(table):
    if table not in ("kinematics", "structural"):
        raise ValueError("table must be 'kinematics' or 'structural'")
    return _get_cache()[table]


def _resolve_clusters(clusters, table="kinematics"):
    cache = _get_cache()
    if table == "kinematics":
        index = cache["kin_index"]
        names = cache["kinematics"]["columns"]["cluster"].astype(str)
    else:
        index = cache["st_index"]
        names = cache["structural"]["columns"]["cluster"].astype(str)

    if clusters is None:
        return np.arange(len(names), dtype=int)

    if isinstance(clusters, str):
        clusters = [clusters]

    out = []
    for name in clusters:
        key = _normalize_name(name)
        if key not in index:
            raise KeyError(f"Unknown {table} cluster: {name}")
        out.append(index[key])
    return np.array(out, dtype=int)


def _slice_dict_columns(table, idx, fields=None):
    cat = _get_table(table)
    if fields is None:
        fields = cat["canon_headers"]

    out = {}
    for field in fields:
        out[field] = np.asarray(cat["columns"][field])[idx]
    return out


def _resolve_column_key(table, header, occurrence=1):
    cat = _get_table(table)
    if header in cat["columns"]:
        return header

    keys = cat["raw_to_canon"].get(header)
    if not keys:
        raise KeyError(f"Unknown column '{header}' in {table}")
    if occurrence < 1 or occurrence > len(keys):
        raise IndexError(
            f"occurrence must be in 1..{len(keys)} for column '{header}'"
        )
    return keys[occurrence - 1]


def _rng(seed=None, rng=None):
    if seed is not None and rng is not None:
        raise ValueError("Pass either seed or rng, not both")
    if rng is not None:
        return rng
    return np.random.default_rng(seed)


def _normalize_icrs_sample_outtype(outtype):
    aliases = {
        "skycoords": "skycoords",
        "skycoordinates": "skycoords",
        "astropy": "skycoords",
        "dict": "dict",
        "dictionary": "dict",
        "array": "array",
    }
    key = str(outtype).strip().lower()
    if key not in aliases:
        raise ValueError(
            "outtype must be one of: 'dict', 'array', 'skycoords', "
            "'skycoordinates', 'astropy'"
        )
    return aliases[key]


def _icrs_dict_from_sample_matrix(sample_matrix):
    from astropy import units as u

    kin_units = units("kinematics")
    return {
        baumgardt_to_astropy_headers["ra"]: sample_matrix[:, 0] * u.Unit(kin_units["ra"]),
        baumgardt_to_astropy_headers["dec"]: sample_matrix[:, 1] * u.Unit(kin_units["dec"]),
        baumgardt_to_astropy_headers["rsun"]: sample_matrix[:, 2] * u.Unit(kin_units["rsun"]),
        baumgardt_to_astropy_headers["mualpha"]: sample_matrix[:, 3] * u.Unit(kin_units["mualpha"]),
        baumgardt_to_astropy_headers["mu_delta"]: sample_matrix[:, 4] * u.Unit(kin_units["mu_delta"]),
        baumgardt_to_astropy_headers["rv"]: sample_matrix[:, 5] * u.Unit(kin_units["rv"]),
    }


def _icrs_means(clusters=None):
    idx = _resolve_clusters(clusters, table="kinematics")
    cat = _get_table("kinematics")["columns"]
    means = np.column_stack(
        [
            np.asarray(cat["ra"])[idx],
            np.asarray(cat["dec"])[idx],
            np.asarray(cat["rsun"])[idx],
            np.asarray(cat["mualpha"])[idx],
            np.asarray(cat["mu_delta"])[idx],
            np.asarray(cat["rv"])[idx],
        ]
    )
    return idx, means


def _kinematic_covariance_for_index(i):
    cat = _get_table("kinematics")["columns"]

    s_dist = float(cat["delta_r"][i])
    s_pmra = float(cat["del_mu"][i])
    s_pmdec = float(cat["del_mu_2"][i])
    s_rv = float(cat["erv"][i])
    corr = float(cat["corr"][i])

    cov = np.zeros((6, 6), dtype=float)
    cov[2, 2] = s_dist * s_dist
    cov[3, 3] = s_pmra * s_pmra
    cov[4, 4] = s_pmdec * s_pmdec
    cov[5, 5] = s_rv * s_rv
    cov[3, 4] = corr * s_pmra * s_pmdec
    cov[4, 3] = cov[3, 4]
    return cov


def names():
    """Return all cluster names from the kinematics catalog."""
    return list(_get_cache()["kinematics"]["columns"]["cluster"].astype(str))


def columns(table):
    """Return canonical column names for a table."""
    return list(_get_table(table)["canon_headers"])


def units(table):
    """Return unit strings keyed by canonical column names."""
    cat = _get_table(table)
    return {
        canon: unit
        for canon, unit in zip(cat["canon_headers"], cat["raw_units"])
    }


def column(table, header, occurrence=1, clusters=None):
    """Return one column by raw header or canonical key.

    Duplicate raw headers can be selected with `occurrence`.
    """
    key = _resolve_column_key(table, header, occurrence=occurrence)
    idx = _resolve_clusters(clusters, table=table)
    return np.asarray(_get_table(table)["columns"][key])[idx]


def kinematics(clusters=None, fields=None):
    """Return selected kinematics columns as plain numpy arrays."""
    idx = _resolve_clusters(clusters, table="kinematics")
    return _slice_dict_columns("kinematics", idx, fields=fields)


def structural(clusters=None, fields=None):
    """Return selected structural columns as plain numpy arrays.

    If `fields` is None, all structural columns are returned.
    """
    idx = _resolve_clusters(clusters, table="structural")
    return _slice_dict_columns("structural", idx, fields=fields)


def icrs(clusters=None):
    """Return astropy-ready ICRS matrix from kinematics.

    Output shape is (N, 6) with columns in this order:
    [ra_deg, dec_deg, distance_kpc, pm_ra_masyr, pm_dec_masyr, rv_kms].
    """
    from astropy import units as u
    _, means = _icrs_means(clusters)
    # instead, return a dictionary ready for astropy 
    coordinates = {
        baumgardt_to_astropy_headers["ra"]: means[:,0] * u.Unit(units("kinematics")['ra']),
        baumgardt_to_astropy_headers["dec"]: means[:,1] * u.Unit(units("kinematics")['dec']),
        baumgardt_to_astropy_headers["rsun"]: means[:,2] * u.Unit(units("kinematics")['rsun']),
        baumgardt_to_astropy_headers["mualpha"]: means[:,3] * u.Unit(units("kinematics")['mualpha']),
        baumgardt_to_astropy_headers["mu_delta"]: means[:,4] * u.Unit(units("kinematics")['mu_delta']),
        baumgardt_to_astropy_headers["rv"]: means[:,5] * u.Unit(units("kinematics")['rv']),
    }    
    return coordinates


def kinematic_covariance(cluster):
    """Return 6x6 covariance matrix for (ra, dec, distance, pm_ra, pm_dec, rv)."""
    idx = _resolve_clusters(cluster, table="kinematics")
    if idx.size != 1:
        raise ValueError("covariance expects exactly one cluster")
    return _kinematic_covariance_for_index(int(idx[0]))


def icrs_sample(n_samples, clusters=None, seed=None, rng=None):
    """Sample ICRS phase-space vectors from per-cluster Gaussian errors.

    `seed` provides reproducible output. Pass `rng` for external RNG control.

    RETURN SHAPE: 
    ARRAY
        [Nsample X Ncuster X 6 (phase_space)]
    DICT
        {sample: {"ra": []*u.deg}}
    """

    if n_samples < 1:
        raise ValueError("n_samples must be >= 1")

    samples = []
    idx, means = _icrs_means(clusters)

    gen = _rng(seed=seed, rng=rng)
    for i in range(means.shape[0]):
        c = _kinematic_covariance_for_index(int(idx[i]))
        samples.append(gen.multivariate_normal(mean=means[i], cov=c, size=n_samples))

    out = np.stack(samples, axis=1)

    if means.shape[0] == 1:
        return out[:, 0, :]
    return out


def icrs_sample_catalog(n_samples, outtype="dict", seed=None, rng=None):
    """
        wrapper to icrs_sample and returns the entire cataloge

    """
    normalized_outtype = _normalize_icrs_sample_outtype(outtype)
    samples = icrs_sample(n_samples=n_samples, clusters=None, seed=seed, rng=rng)
    if normalized_outtype == "array":
        return samples

    out = {i: _icrs_dict_from_sample_matrix(samples[i, :, :]) for i in range(n_samples)}
    if normalized_outtype == "dict":
        return out

    from astropy import coordinates
    return [coordinates.SkyCoord(**out[i]) for i in range(n_samples)]


def icrs_sample_cluster(n_samples, cluster, outtype="dict", seed=None, rng=None):
    """
    shape will now be 
    """
    normalized_outtype = _normalize_icrs_sample_outtype(outtype)
    if not isinstance(cluster, str):
        raise TypeError("cluster must be a string")

    samples = icrs_sample(n_samples=n_samples, clusters=cluster, seed=seed, rng=rng)
    if normalized_outtype == "array":
        return samples

    out = _icrs_dict_from_sample_matrix(samples)
    if normalized_outtype == "dict":
        return out

    from astropy import coordinates
    return coordinates.SkyCoord(**out)


def sample_structural(n_samples, key, clusters=None, seed=None, rng=None):
    """the structural parameters are not correlated"""

    _valid_keys_with_errors = _structural_columns_with_errors.keys()
    key = _resolve_column_key("structural", key)
    if not isinstance(key, str):
        raise ValueError("key must be string but was instead {:s}".format(str(type(key))))
    if not key in _valid_keys_with_errors:
        raise ValueError("{:s} does not have an assocaited uncertainty ".format(key))
    error_key = _structural_columns_with_errors[key]
    value = column("structural", key, occurrence=1, clusters=clusters)
    uncertainty = column("structural", error_key, occurrence=1,clusters=clusters)
    gen = _rng(seed=seed, rng=rng)
    samples=gen.normal(value,uncertainty, size=(n_samples, len(value)))
    return samples


_structural_columns_with_errors = {}
for key in _baumgardt_structural_columns_with_errors:
    _structural_columns_with_errors[_normalize_column_name(key)] = _normalize_column_name(_baumgardt_structural_columns_with_errors[key])

__all__ = [
    "names",
    "columns",
    "units",
    "column",
    "kinematics",
    "structural",
    "icrs",
    "kinematic_covariance",
    "icrs_sample",
    "icrs_sample_catalog",
    "icrs_sample_cluster",
]
