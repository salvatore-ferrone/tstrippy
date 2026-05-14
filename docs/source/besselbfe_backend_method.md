# Bessel-BFE Backend: Poisson Solver Method (Axisymmetric)

This note documents the exact mathematical method implemented in `tstrippy/src/besselbfe.f90` for solving Poisson's equation for axisymmetric density profiles.

It is intended to be readable by both humans and AI tools, with explicit equations, implementation mapping, and numerical caveats.

## 1. Problem and assumptions

We solve Poisson's equation for an axisymmetric density field:

$$
\nabla^2 \Phi(R,z) = 4\pi G\,\rho(R,z),
$$

with no azimuthal dependence ($\partial/\partial\phi = 0$), so

$$
\frac{\partial^2 \Phi}{\partial R^2}
+ \frac{1}{R}\frac{\partial \Phi}{\partial R}
+ \frac{\partial^2 \Phi}{\partial z^2}
= 4\pi G\,\rho(R,z).
$$

Implementation assumptions:

1. Axisymmetry.
2. Density is sampled on $(R,z)$ with $z \ge 0$ and then extended evenly in $z$.
3. Potential and force are tabulated per component and queried by interpolation.

## 2. Why Bessel/Hankel transforms

For axisymmetric problems, the radial Laplacian is diagonalized by Bessel functions. Using the order-0 Hankel transform:

$$
\tilde f(k,z) = \int_0^\infty f(R,z) J_0(kR)\,R\,dR,
$$
$$
f(R,z) = \int_0^\infty \tilde f(k,z) J_0(kR)\,k\,dk.
$$

The transformed Poisson equation becomes a family of 1D equations in $z$:

$$
\frac{\partial^2 \tilde\Phi}{\partial z^2} - k^2\tilde\Phi = 4\pi G\,\tilde\rho(k,z).
$$

This is the key simplification: one 2D PDE becomes many 1D ODEs indexed by $k$.

## 3. Vertical Green-function solution

For each $k$, the Green-function solution is

$$
\tilde\Phi(k,z)
= -\frac{2\pi G}{k}
\int_{-\infty}^{\infty} \tilde\rho(k,z') e^{-k|z-z'|}\,dz'.
$$

With even extension in $z$ and evaluation for $z\ge 0$, this becomes

$$
\int_0^{\infty} \tilde\rho(k,z')
\left[e^{-k|z-z'|} + e^{-k(z+z')}\right] dz'.
$$

In code this quantity is accumulated as `zconv_kz(j, iZ)`.

## 4. Inverse transform for potential

The potential is reconstructed by

$$
\Phi(R,z) = -2\pi G \int_0^\infty J_0(kR)\,\mathcal{Z}(k,z)\,dk,
$$

where $\mathcal{Z}(k,z)$ is the vertical convolution above.

In implementation, this is approximated by bounded Gauss-Legendre quadrature over $k \in [0,k_{\max}]$:

$$
\Phi(R,z) \approx -2\pi G \sum_j w_j J_0(k_j R)\,\mathcal{Z}(k_j,z).
$$

Code mapping:

1. `k_grid(j)`, `quad_w(j)` hold $k_j, w_j$.
2. `BESSEL_TABLE_PHI(iR,iZ,ic)` stores tabulated $\Phi$.

## 5. Transform-space derivatives (force path)

The backend now computes force derivatives in transform space, not by finite-differencing the potential table.

### 5.1 Radial derivative

Using

$$
\frac{d}{dR}J_0(kR) = -kJ_1(kR),
$$

we get

$$
\frac{\partial \Phi}{\partial R}
= +2\pi G \int_0^\infty k J_1(kR)\,\mathcal{Z}(k,z)\,dk
\approx +2\pi G \sum_j w_j k_j J_1(k_j R)\,\mathcal{Z}(k_j,z).
$$

Stored in `BESSEL_TABLE_DPHI_DR`.

### 5.2 Vertical derivative

Differentiate the vertical kernel:

$$
\frac{\partial}{\partial z} e^{-k|z-z'|} = -k\,\mathrm{sgn}(z-z')\,e^{-k|z-z'|},
$$
$$
\frac{\partial}{\partial z} e^{-k(z+z')} = -k e^{-k(z+z')}.
$$

Define

$$
\mathcal{Z}_z(k,z) = \frac{\partial}{\partial z}\mathcal{Z}(k,z),
$$

then

$$
\frac{\partial \Phi}{\partial z}
= -2\pi G \int_0^\infty J_0(kR)\,\mathcal{Z}_z(k,z)\,dk
\approx -2\pi G \sum_j w_j J_0(k_j R)\,\mathcal{Z}_z(k_j,z).
$$

Stored in `BESSEL_TABLE_DPHI_DZ`.

## 6. Grid and metadata

Per component, the backend stores metadata in `BESSEL_TABLE_META`:

1. `meta(1) = logR_min`
2. `meta(2) = dlogR`
3. `meta(3) = dz`
4. `meta(4) = mass_est` (table-estimated total mass for far-field closure)

Radial grid is logarithmic:

$$
R_i = \exp\left(\log R_{\min} + (i-1)\,d\log R\right),
$$

with

$$
R_{\min}=10^{-3} r_{\mathrm{scale}},\qquad
R_{\max}=2\times10^2\,r_{\mathrm{scale}}.
$$

Vertical grid is linear over $z\in[0,z_{\max}]$ with

$$
z_{\max}=2\times10^2\,z_{\mathrm{scale}}.
$$

## 7. Runtime interpolation and force assembly

At query time, the active component is interpolated on the local cell using bilinear interpolation for robustness.

Interpolate:

1. $\Phi$
2. $\partial\Phi/\partial R$
3. $\partial\Phi/\partial |z|$

Then assemble accelerations:

$$
a_R = -\frac{\partial\Phi}{\partial R},
$$
$$
a_x = a_R\frac{x}{R},\quad a_y = a_R\frac{y}{R},
$$
$$
a_z = -\operatorname{sgn}(z)\,\frac{\partial\Phi}{\partial |z|}.
$$

## 8. Far-field stabilization (monopole blend)

To suppress non-physical far-field oscillations from finite-domain and finite-$k$ truncation effects, the backend blends to monopole closure:

$$
\Phi_{\mathrm{mono}} = -\frac{G M_{\mathrm{est}}}{r},
\qquad
\mathbf{a}_{\mathrm{mono}} = -\frac{G M_{\mathrm{est}}}{r^3}\mathbf{r}.
$$

Blend is smoothstep in spherical radius $r=\sqrt{R^2+z^2}$:

$$
t = \mathrm{clamp}\!\left(\frac{r-r_0}{r_1-r_0},0,1\right),
\qquad
w(t)=t^2(3-2t).
$$

Current implementation uses

$$
r_0 = 8\,s_{\mathrm{ref}},\qquad r_1 = 20\,s_{\mathrm{ref}},
\qquad s_{\mathrm{ref}}=\max(r_{\mathrm{scale}},z_{\mathrm{scale}}).
$$

Blended fields:

$$
\Phi \leftarrow (1-w)\Phi + w\Phi_{\mathrm{mono}},
$$
$$
\mathbf{a} \leftarrow (1-w)\mathbf{a} + w\mathbf{a}_{\mathrm{mono}}.
$$

## 9. Component scale control

Each component has independent domain scales:

1. `BESSEL_COMPONENT_R_SCALE(ic)`
2. `BESSEL_COMPONENT_Z_SCALE(ic)`

This supports mixed-scale multi-component models in one backend instance.

## 10. Numerical complexity (high level)

For each component, table build is dominated by:

1. Density sampling over $(N_R, N_Z)$
2. Hankel transform loops over $(N_k, N_Z, N_R)$
3. Vertical convolution loops over $(N_k, N_Z, N_Z)$
4. Final reconstruction over $(N_R, N_Z, N_k)$

So asymptotically it is expensive at build time, but cheap at runtime queries due to interpolation.

## 11. Known caveats

1. Finite $k$ truncation can still bias very large-radius behavior without stabilization.
2. Accuracy depends on consistent scale choices and query envelope relative to table domain.
3. Mixed-component extreme scale separation can require larger grids and/or tighter diagnostics.

## 12. Practical validation checks

Recommended checks for each model configuration:

1. Force-potential consistency via finite differences at interior points.
2. Symmetry/parity checks across $z\to -z$.
3. Far-field directionality:
   $$
   \cos\theta = \frac{\mathbf{a}\cdot(-\hat{\mathbf r})}{\|\mathbf a\|} \to 1.
   $$
4. Far-field monopole consistency:
   $$
   \Phi \sim -\frac{GM}{r},\quad a_r\sim -\frac{GM}{r^2}.
   $$

## 13. Mapping to implementation symbols

Primary file: `tstrippy/src/besselbfe.f90`

Key routines:

1. `initialize`, `default_initialize`, `allocate_component_tables`
2. `set_component_scales`, `project_density`, `load_component`
3. `build_table` (all transforms and derivative-table generation)
4. `component_force_potential` (runtime interpolation + force assembly + far-field blend)

---

If you share this method externally, include both this note and a pinned commit hash so reviewers can match equations to code exactly.