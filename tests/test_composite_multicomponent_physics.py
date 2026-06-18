import numpy as np

import tstrippy


def _composite_components():
    return [
        ("miyamotonagai", [5.0e10, 6.5, 0.26]),
        ("exponentialdisk", [1.1e9, 2.8, 0.35]),
        ("ibata2024halo", [2.0e7, 1.2, 45.0, 0.85, 1.3, 3.4]),
        ("exponentialoblatehalo", [3.5e7, 8.0, 0.9]),
        ("exponentialdisk", [6.0e8, 1.4, 0.22]),
    ]


def _build_and_eval_force(component_list, x, y, z):
    g = tstrippy.gravity
    g.clear()

    for model_name, params in component_list:
        g.add_component(model_name, params)

    g.finalize()
    ax, ay, az = g.force(x, y, z)
    return ax, ay, az


def _stack_force(ax, ay, az):
    return np.column_stack((ax, ay, az))


def test_composite_force_component_sum_matches_net_force():
    g = tstrippy.gravity
    components = _composite_components()

    x = np.array([1.0, 2.5, 5.0, 9.0], dtype=float)
    y = np.array([0.0, 0.8, -0.4, 0.3], dtype=float)
    z = np.array([0.1, -0.2, 0.5, -0.7], dtype=float)

    g.clear()
    for model_name, params in components:
        g.add_component(model_name, params)
    g.finalize()

    ax_net, ay_net, az_net = g.force(x, y, z)
    ax_comp, ay_comp, az_comp = g.force_components(x, y, z)

    ncomp = int(g.gravity_ncomp)
    np.testing.assert_equal(ncomp, len(components))

    np.testing.assert_allclose(ax_comp[:ncomp, :].sum(axis=0), ax_net, rtol=1e-8, atol=1e-8)
    np.testing.assert_allclose(ay_comp[:ncomp, :].sum(axis=0), ay_net, rtol=1e-8, atol=1e-8)
    np.testing.assert_allclose(az_comp[:ncomp, :].sum(axis=0), az_net, rtol=1e-8, atol=1e-8)


def test_composite_force_is_order_independent():
    components = _composite_components()

    x = np.array([1.0, 2.5, 5.0, 9.0], dtype=float)
    y = np.array([0.0, 0.8, -0.4, 0.3], dtype=float)
    z = np.array([0.1, -0.2, 0.5, -0.7], dtype=float)

    ax_a, ay_a, az_a = _build_and_eval_force(components, x, y, z)
    ax_b, ay_b, az_b = _build_and_eval_force(list(reversed(components)), x, y, z)

    np.testing.assert_allclose(ax_a, ax_b, rtol=1e-8, atol=1e-8)
    np.testing.assert_allclose(ay_a, ay_b, rtol=1e-8, atol=1e-8)
    np.testing.assert_allclose(az_a, az_b, rtol=1e-8, atol=1e-8)


def test_composite_force_matches_sum_of_independent_components():
    g = tstrippy.gravity
    components = _composite_components()

    x = np.array([1.0, 2.5, 5.0, 9.0], dtype=float)
    y = np.array([0.0, 0.8, -0.4, 0.3], dtype=float)
    z = np.array([0.1, -0.2, 0.5, -0.7], dtype=float)

    g.clear()
    for model_name, params in components:
        g.add_component(model_name, params)
    g.finalize()
    ax_net, ay_net, az_net = g.force(x, y, z)
    force_net = _stack_force(ax_net, ay_net, az_net)

    force_sum = np.zeros_like(force_net)
    for model_name, params in components:
        g.clear()
        g.add_component(model_name, params)
        g.finalize()
        ax_i, ay_i, az_i = g.force(x, y, z)
        force_sum += _stack_force(ax_i, ay_i, az_i)

    np.testing.assert_allclose(force_sum, force_net, rtol=1e-8, atol=1e-8)
