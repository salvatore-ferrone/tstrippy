#!/usr/bin/env python3
"""Generate quick diagnostic plots for Bessel backend physics checks.

Usage:
    conda run -n tstrippy python tests/diagnostics/plot_bessel_physics_diagnostics.py

Optional args:
    --output-dir tests/diagnostics/figures
"""

from __future__ import annotations

import argparse
import math
from pathlib import Path
import sys

import matplotlib.pyplot as plt
import numpy as np

# Prefer the installed/runtime package first. Fall back to repo-root import.
try:
    import tstrippy
except ModuleNotFoundError:
    REPO_ROOT = Path(__file__).resolve().parents[2]
    if str(REPO_ROOT) not in sys.path:
        sys.path.insert(0, str(REPO_ROOT))
    import tstrippy


def configure_single_bessel_disk(sigma0=1.0, hR=4.0, hZ=0.8, nr=192, nz=160, nk=384):
    g = tstrippy.gravity
    g.clear()
    g.bessel_initialize(nr, nz, nk)
    g.add_component("exponentialdisk", [sigma0, hR, hZ])
    g.bessel_set_component_scales(1, hR, hZ)
    g.finalize()
    return g


def configure_two_bessel_disks(
    sigma0_1=2.0,
    hR_1=4.0,
    hZ_1=2.0,
    sigma0_2=1.0,
    hR_2=6.0,
    hZ_2=0.5,
    nr=256,
    nz=256,
    nk=256,
    scale_divisor=2.0,
):
    g = tstrippy.gravity
    g.clear()
    g.bessel_initialize(nr, nz, nk)
    g.add_component("exponentialdisk", [sigma0_1, hR_1, hZ_1])
    g.add_component("exponentialdisk", [sigma0_2, hR_2, hZ_2])
    g.bessel_set_component_scales(1, hR_1 / scale_divisor, hZ_1 / scale_divisor)
    g.bessel_set_component_scales(2, hR_2 / scale_divisor, hZ_2 / scale_divisor)
    g.finalize()
    return g


def disk_total_mass(sigma0, hR):
    return 2.0 * math.pi * sigma0 * hR**2


def expdisk_density(R, z, sigma0, hR, hZ):
    return (sigma0 / (2.0 * hZ)) * np.exp(-R / hR) * np.exp(-np.abs(z) / hZ)


def sample_alignment(g, radii, theta):
    x = radii * np.cos(theta)
    y = np.zeros_like(x)
    z = radii * np.sin(theta)

    ax, ay, az = g.force(x, y, z)
    force = np.column_stack((ax, ay, az))
    pos = np.column_stack((x, y, z))

    minus_rhat = -pos / np.linalg.norm(pos, axis=1, keepdims=True)
    force_norm = np.linalg.norm(force, axis=1)
    cos_align = np.sum(force * minus_rhat, axis=1) / np.maximum(force_norm, 1.0e-30)

    return cos_align


def far_field_inferred_masses(g, radii, theta):
    gconst = float(tstrippy.gravity.gravity_g_default)

    x = radii * np.cos(theta)
    y = np.zeros_like(x)
    z = radii * np.sin(theta)

    phi = g.potential(x, y, z)
    ax, ay, az = g.force(x, y, z)

    pos = np.column_stack((x, y, z))
    rhat = pos / np.linalg.norm(pos, axis=1, keepdims=True)
    acc = np.column_stack((ax, ay, az))
    a_r = np.sum(acc * rhat, axis=1)

    m_phi = -radii * phi / gconst
    m_a = -(radii**2) * a_r / gconst
    return m_phi, m_a


def poisson_residual_map(g, sigma0, hR, hZ, npoints=100):
    gconst = float(tstrippy.gravity.gravity_g_default)


    R = np.linspace(0.6, 8.0, npoints)
    Z = np.linspace(-3.0, 3.0, npoints)
    RR, ZZ = np.meshgrid(R, Z, indexing="ij")

    x = RR.ravel()
    y = np.zeros_like(x)
    z = ZZ.ravel()

    phi = g.potential(x, y, z).reshape(RR.shape)

    dR = R[1] - R[0]
    dZ = Z[1] - Z[0]

    phi_rr = (phi[2:, 1:-1] - 2.0 * phi[1:-1, 1:-1] + phi[:-2, 1:-1]) / (dR**2)
    phi_r = (phi[2:, 1:-1] - phi[:-2, 1:-1]) / (2.0 * dR)
    phi_zz = (phi[1:-1, 2:] - 2.0 * phi[1:-1, 1:-1] + phi[1:-1, :-2]) / (dZ**2)

    R_inner = RR[1:-1, 1:-1]
    Z_inner = ZZ[1:-1, 1:-1]

    laplacian = phi_rr + (phi_r / R_inner) + phi_zz
    rho = expdisk_density(R_inner, Z_inner, sigma0, hR, hZ)
    source = 4.0 * math.pi * gconst * rho

    residual = laplacian - source
    rel = np.abs(residual) / np.maximum(np.abs(source), 1.0e-12)

    return R_inner, Z_inner, rel, source


def plot_far_field_alignment(output_dir: Path):
    g_single = configure_single_bessel_disk()
    g_two = configure_two_bessel_disks()

    radii_single = np.array([20.0, 30.0, 40.0, 60.0, 80.0, 120.0], dtype=float)
    radii_two = np.array([60.0, 90.0, 120.0, 180.0, 240.0, 320.0], dtype=float)

    cos_single = sample_alignment(g_single, radii_single, theta=0.35)
    cos_two = sample_alignment(g_two, radii_two, theta=0.30)

    fig, axis = plt.subplots(1, 1, figsize=(7, 4.5))
    axis.plot(radii_single, cos_single, "o-", label="single disk", linewidth=2)
    axis.plot(radii_two, cos_two, "s-", label="two disks", linewidth=2)
    axis.axhline(1.0, color="k", linestyle="--", linewidth=1, alpha=0.5)
    axis.set_xlabel("radius r")
    axis.set_ylabel("cos(theta) with -rhat")
    axis.set_ylim(0.95, 1.002)
    axis.set_title("Far-field force directionality")
    axis.grid(alpha=0.3)
    axis.legend()
    fig.tight_layout()
    fig.savefig(output_dir / "far_field_alignment.png", dpi=220)
    plt.close(fig)

    print("[alignment] single min cos:", float(np.min(cos_single)))
    print("[alignment] two-disk min cos:", float(np.min(cos_two)))


def plot_inferred_mass_curves(output_dir: Path):
    sigma0, hR, hZ = 1.0, 4.0, 0.8
    g = configure_single_bessel_disk(sigma0=sigma0, hR=hR, hZ=hZ)

    radii = np.array([30.0, 40.0, 50.0, 70.0, 90.0, 120.0, 150.0, 180.0], dtype=float)
    m_phi, m_a = far_field_inferred_masses(g, radii, theta=0.33)
    m_true = disk_total_mass(sigma0, hR)

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5), gridspec_kw={"width_ratios": [2.2, 1.2]})

    axes[0].plot(radii, m_phi, "o-", label="M_phi = -r Phi / G", linewidth=2)
    axes[0].plot(radii, m_a, "s-", label="M_a = -r^2 a_r / G", linewidth=2)
    axes[0].axhline(m_true, color="k", linestyle="--", linewidth=1.2, label="M_true")
    axes[0].set_xlabel("radius r")
    axes[0].set_ylabel("inferred mass")
    axes[0].set_title("Far-field inferred mass curves")
    axes[0].grid(alpha=0.3)
    axes[0].legend()

    err_phi_pct = 100.0 * (m_phi - m_true) / m_true
    err_a_pct = 100.0 * (m_a - m_true) / m_true
    axes[1].plot(radii, err_phi_pct, "o-", linewidth=2, label="phi")
    axes[1].plot(radii, err_a_pct, "s-", linewidth=2, label="a_r")
    axes[1].axhline(0.0, color="k", linestyle="--", linewidth=1)
    axes[1].set_xlabel("radius r")
    axes[1].set_ylabel("mass error [%]")
    axes[1].set_title("Bias vs true mass")
    axes[1].grid(alpha=0.3)
    axes[1].legend()

    fig.tight_layout()
    fig.savefig(output_dir / "far_field_inferred_mass.png", dpi=220)
    plt.close(fig)

    print("[mass] M_true:", float(m_true))
    print("[mass] largest-r M_phi error [%]:", float(err_phi_pct[-1]))
    print("[mass] largest-r M_a error [%]:", float(err_a_pct[-1]))


def plot_poisson_residual_map(output_dir: Path):
    sigma0, hR, hZ = 1.0, 4.0, 0.8
    g = configure_single_bessel_disk(sigma0=sigma0, hR=hR, hZ=hZ)

    R, Z, rel, source = poisson_residual_map(g, sigma0, hR, hZ)

    source_mask = np.abs(source) > 1.0e-3 * np.nanmax(np.abs(source))
    dz = np.abs(Z[0, 1] - Z[0, 0])
    midplane_mask = np.abs(Z) > 2.0 * dz
    rel_masked = np.where(source_mask & midplane_mask, rel, np.nan)

    rel_plot = np.log10(np.maximum(rel_masked, 1.0e-6))

    fig, axis = plt.subplots(1, 1, figsize=(7, 4.8))
    mesh = axis.pcolormesh(R, Z, rel_plot, shading="auto", cmap="magma", vmin=-1.0, vmax=0.7)
    cbar = fig.colorbar(mesh, ax=axis)
    cbar.set_label("log10(relative Poisson residual)")

    axis.set_xlabel("R")
    axis.set_ylabel("z")
    axis.set_title("Interior Poisson residual map (masked)")
    axis.set_aspect("auto")
    fig.tight_layout()
    fig.savefig(output_dir / "poisson_residual_map.png", dpi=220)
    plt.close(fig)

    rel_vals = rel_masked[np.isfinite(rel_masked)]
    if rel_vals.size > 0:
        print("[residual] median:", float(np.nanmedian(rel_vals)))
        print("[residual] p95:", float(np.nanpercentile(rel_vals, 95.0)))
        print("[residual] p99:", float(np.nanpercentile(rel_vals, 99.0)))


def parse_args():
    parser = argparse.ArgumentParser(description="Plot Bessel backend diagnostics used by tests.")
    parser.add_argument(
        "--output-dir",
        default="tests/diagnostics/figures",
        help="Directory where PNG plots are saved.",
    )
    return parser.parse_args()


def main():
    args = parse_args()
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    plot_far_field_alignment(output_dir)
    plot_inferred_mass_curves(output_dir)
    plot_poisson_residual_map(output_dir)

    print(f"Saved diagnostics to: {output_dir}")
    print(" - far_field_alignment.png")
    print(" - far_field_inferred_mass.png")
    print(" - poisson_residual_map.png")


if __name__ == "__main__":
    main()
