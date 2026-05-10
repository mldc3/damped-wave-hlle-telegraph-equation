# Numerical Method and Implementation Workflow

This document explains how the damped wave / telegraph equation workflows are implemented in this repository and how each implementation stage maps to the theoretical model.

The main execution script is `intentofinalpractica6.py`. The file in `src/damped_wave_telegraph_methods.py` is preserved as uploaded and is not modified by this documentation rebuild.

## 1. Spatial grid

The script defines a one-dimensional interval $[x_{\min},x_{\max}]$ and uses a uniform mesh with `Nx` points:

$$
x_j=x_{\min}+j\Delta x,\qquad \Delta x=\frac{x_{\max}-x_{\min}}{N_x-1}.
$$

In code, this is represented by `x = np.linspace(x_min, x_max, Nx)` and reused by all methods for consistent comparisons.

## 2. Time grid

The simulation horizon is set by `t_max` with `Nt` temporal points and

$$
\Delta t=\frac{t_{\max}}{N_t-1}.
$$

All four numerical strategies use the same grid parameters, enabling method-to-method interpretation under matched physical and numerical settings.

## 3. Initialization of $u$ and $u_t$

The displacement is initialized with a sinusoidal mode through `initial_u(n)`, equivalent to

$$
u(x,0)=\sin(n\pi x).
$$

A startup state at one timestep is then generated in `u_at_dt(u0)`, using a second-order expansion based on the PDE at rest-like initial velocity. This provides the two temporal levels required by the explicit second-order recurrence.

## 4. Explicit finite-difference update

Functions:

- `explicit_method_dirichlet(...)`
- `explicit_method_neumann(...)`

The update is a three-level recurrence over interior points using the centered spatial Laplacian. The dimensionless coefficient

$$
r=\left(\frac{c\Delta t}{\Delta x}\right)^2
$$

controls wave transport strength in the recurrence.

Boundary handling differs by function:

- Dirichlet: endpoint values are fixed to zero each step.
- Neumann setup: right boundary is imposed by copying the adjacent interior value (zero-gradient approximation).

## 5. Construction of the first-order matrix system

Functions:

- `build_system_matrix_dirichlet(...)`
- `build_system_matrix_neumann(...)`

The second-order PDE is reformulated with state $\mathbf{W}=(\mathbf{u},\mathbf{v})^T$, producing

$$
\frac{d\mathbf{W}}{dt}=\mathbf{A}\mathbf{W}.
$$

The code builds block matrices with discrete Laplacian terms, damping block $-2(\kappa/\rho)\mathbf{I}$, and boundary-specific row constraints.

## 6. Matrix exponential workflow

Function:

- `solve_matrix_exponential(...)`

For each boundary type, the script constructs $\mathbf{A}$ and computes

$$
\mathbf{M}=e^{\mathbf{A}\Delta t}
$$

via `scipy.linalg.expm`. The state vector is then advanced by repeated matrix multiplication. This acts as a high-fidelity reference for the same semi-discrete operator.

## 7. Crank–Nicolson workflow

Functions:

- `solve_crank_nicolson_dirichlet(...)`
- `solve_crank_nicolson_neumann(...)`

The code forms

$$
\mathbf{B}=\mathbf{I}-\frac{\Delta t}{2}\mathbf{A},\qquad
\mathbf{C}=\mathbf{I}+\frac{\Delta t}{2}\mathbf{A},
$$

and solves

$$
\mathbf{B}\mathbf{W}^{n+1}=\mathbf{C}\mathbf{W}^{n}
$$

at each step using LU factorization (`lu_factor` / `lu_solve`). Boundary rows are explicitly fixed in the linear system setup.

## 8. Conservative HLL/HLLE workflow

Functions:

- `compute_eigenvalues(...)`
- `compute_hll_flux(...)`
- `solve_hll_dirichlet(...)`
- `solve_hll_neumann(...)`

This branch uses the conservative-like state $\mathbf{U}=(u,v,w)^T$ with $w\approx u_x$. At each step:

1. Reconstruct $w$ from spatial differences.
2. Compute interface fluxes with HLL speed bounds.
3. Update interior cells with flux divergence plus source terms.
4. Apply boundary conditions on the updated state.

The implementation uses characteristic speeds tied to $\pm c$ and includes damping/source effects in the source vector.

## 9. Dirichlet boundary implementation

Dirichlet behavior is enforced by fixed displacement at constrained ends:

- Explicit method sets endpoint values directly to zero.
- Matrix-system methods impose boundary equations by replacing corresponding rows.
- HLL branch sets edge state components associated with fixed displacement/velocity to zero.

This corresponds to fixed-end conditions in the project comparisons.

## 10. Neumann boundary implementation

The Neumann comparison uses a zero-gradient/free-end treatment at the right boundary in the implementation:

- Explicit method uses endpoint copying from the neighboring interior value.
- Matrix-system Neumann builder modifies the right-boundary Laplacian row with a ghost-point style closure.
- HLL branch sets right-edge gradient accordingly and copies right boundary state from interior in the update step.

The left boundary remains fixed in the script’s mixed setup for the Neumann experiments.

## 11. GIF generation

After computing solution arrays, the script builds six Matplotlib animations (`FuncAnimation`) and writes GIF files using the Pillow writer. Frame skipping is controlled by `SKIP_FRAMES` to reduce runtime and file size.

The six generated outputs correspond to:

1. Explicit (Dirichlet vs Neumann)
2. Matrix exponential (Dirichlet vs Neumann)
3. Crank–Nicolson (Dirichlet vs Neumann)
4. HLL (Dirichlet vs Neumann)
5. Multi-method comparison under Dirichlet
6. Multi-method comparison under Neumann

## 12. Mapping from code structure to theory

- **Second-order PDE viewpoint**: explicit recurrence directly advances $u$ from prior time levels.
- **First-order linear-system viewpoint**: matrix exponential and Crank–Nicolson operate on $\mathbf{W}=(u,v)^T$.
- **Conservative flux viewpoint**: HLL/HLLE branch advances $\mathbf{U}=(u,v,w)^T$ by interface fluxes and source terms.
- **Boundary physics viewpoint**: Dirichlet and Neumann closures are implemented separately for every numerical philosophy.

This separation makes the repository a method-comparison study of the same physical equation under consistent parameters and shared visualization outputs.
