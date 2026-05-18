import numpy as np

from tstrippy.io import mwgcs


def test_names_and_column_metadata():
    nms = mwgcs.names()
    assert len(nms) > 100
    assert "NGC_104" in nms

    kin_cols = mwgcs.columns("kinematics")
    str_cols = mwgcs.columns("structural")
    assert len(kin_cols) == 31
    assert len(str_cols) == 37
    assert "rv" in kin_cols
    assert "rh_m" in str_cols


def test_icrs_and_case_insensitive_lookup():
    a = mwgcs.icrs("NGC_104")
    b = mwgcs.icrs("ngc104")

    assert a.shape == (1, 6)
    assert b.shape == (1, 6)
    np.testing.assert_allclose(a, b)


def test_covariance_shape_and_sanity():
    cov = mwgcs.covariance("NGC104")
    assert cov.shape == (6, 6)
    assert np.all(np.isfinite(cov))
    assert np.all(np.diag(cov) >= 0.0)


def test_icrs_sampling_reproducible_seed():
    s1 = mwgcs.icrs_sample(32, "NGC_104", seed=123)
    s2 = mwgcs.icrs_sample(32, "ngc104", seed=123)

    assert s1.shape == (32, 6)
    np.testing.assert_allclose(s1, s2)
