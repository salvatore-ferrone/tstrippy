Absolutely. Here is a clean Markdown version with KaTeX-friendly equations that should render properly.

## Basis-Expansion Feature Roadmap

### 1. Start On The Feature Branch
```bash
git fetch origin
git checkout basisExpansion
git pull --ff-only origin basisExpansion
```

If you have local uncommitted changes:
```bash
git stash -u
git checkout basisExpansion
git pull --ff-only origin basisExpansion
git stash pop
```

### 2. Define Scope For v1 (Axisymmetric First)
For a first implementation, use only axisymmetric terms ($m=0$). This keeps the method stable and avoids full spherical-harmonic complexity in v1.

Represent density as:
$$
\rho(r,\mu)=\sum_{\ell=0}^{\ell_{\max}}\rho_\ell(r)P_\ell(\mu),
\qquad
\mu=\cos\theta
$$

Compute coefficients:
$$
\rho_\ell(r)=\frac{2\ell+1}{2}\int_{-1}^{1}\rho(r,\mu)\,P_\ell(\mu)\,d\mu
$$

Potential coefficients:
$$
\Phi_\ell(r)= -\frac{4\pi G}{2\ell+1}
\left[
r^{-(\ell+1)}\int_0^r \rho_\ell(r')\,r'^{\ell+2}\,dr'
+
r^\ell\int_r^\infty \rho_\ell(r')\,r'^{1-\ell}\,dr'
\right]
$$

Total potential:
$$
\Phi(r,\theta)=\sum_{\ell=0}^{\ell_{\max}}\Phi_\ell(r)\,P_\ell(\cos\theta)
$$

Then compute accelerations from $a_r=-\partial\Phi/\partial r$ and $a_\theta=-(1/r)\partial\Phi/\partial\theta$, and convert to $(a_x,a_y,a_z)$.

### 3. Python Proof Of Concept Before Fortran
Build a PoC that:
1. Computes $\rho_\ell(r)$ on a radial grid.
2. Computes $\Phi_\ell(r)$ and $d\Phi_\ell/dr$.
3. Evaluates $\Phi,\mathbf{a}$ at random 3D points.
4. Validates gradient consistency: $-\nabla\Phi \approx \mathbf{a}$ (finite differences).
5. Sweeps $(\ell_{\max},N_r)$ to pick an accuracy/runtime target.

### 4. Fortran Design Split
Use mathutils.f90 for numerical primitives:
1. Legendre polynomial recurrence for $P_\ell(\mu)$.
2. Optional derivative helper for $dP_\ell/d\mu$.
3. Robust radial interpolation utilities.

Use potentials.f90 for model evaluation:
1. Add a new potential routine (for example, `flattenedexp_bfe`).
2. Add one initializer routine to precompute/load radial coefficient tables.
3. Runtime path should only interpolate and evaluate (no expensive radial integrals per particle per step).

### 5. Integrator Wiring
In integrator.f90:
1. Add name dispatch in `setstaticgalaxy` for the new model.
2. Add a guard so the code errors early if BFE tables are not initialized.
3. Keep the same potential interface signature as existing models.

### 6. Parser And Data Plumbing
1. Add a parser entry for BFE parameters in Python.
2. Add a config/data file for BFE settings (for example: $r_{\min},r_{\max},N_r,\ell_{\max}$, flattening parameters, table path if persisted).
3. Keep units consistent with the existing package convention.

### 7. Validation And Tests
Add tests for:
1. Force consistency: finite-difference $-\nabla\Phi$ vs returned acceleration.
2. Symmetry checks for axisymmetry.
3. Convergence vs $\ell_{\max}$ and $N_r$.
4. Integrator smoke test using the new potential name.
5. Basic runtime benchmark to catch regressions.

### 8. Milestones
1. M1: Python PoC validated with error plots.
2. M2: `mathutils` Legendre/interpolation utilities complete.
3. M3: `potentials` BFE evaluator complete.
4. M4: Integrator and parser wiring complete.
5. M5: Tests and docs example complete.

If you want, next I can give you a ready-to-paste plan.md template version of this with checkboxes.