# Damped Wave / Telegraph Equation: Explicit, Crank–Nicolson and HLL Solvers

This repository is a cleaned computational-physics portfolio project focused on the one-dimensional damped wave (telegraph) equation. It consolidates the numerical methods, boundary-condition experiments, and visual outputs into a reproducible and presentation-ready structure.

The model equation studied is

$$
\frac{\partial^2 u}{\partial t^2}-c^2\frac{\partial^2 u}{\partial x^2}=-2\frac{\kappa}{\rho}\frac{\partial u}{\partial t}+a u,
$$

where propagation, damping, and linear source effects interact.

## What is implemented

The project compares four numerical philosophies for the same PDE:

- **Explicit finite-difference integration** for the second-order-in-time form.
- **Matrix exponential time evolution** for the semi-discrete first-order linear system.
- **Crank–Nicolson time evolution** as an implicit second-order approximation to the matrix-exponential dynamics.
- **HLL/HLLE conservative flux method** after reformulation into a first-order hyperbolic-like system.

The simulations are shown under both:

- **Dirichlet boundary conditions** (fixed boundary displacement), and
- **Neumann boundary conditions** (zero-gradient/free-end behavior at the right boundary in the implemented setup).

A central interpretation theme is the difference between:

- **physical damping** from the PDE term $-2(\kappa/\rho)u_t$, and
- **numerical damping** introduced by discretization and flux stabilization.

## Available GIF outputs

Only the following six animations are part of this repository:

![Explicit finite-difference method](figures/explicit/explicit_finite_difference.gif)

![Matrix exponential method](figures/matrix_exponential/matrix_exponential_method.gif)

![Crank–Nicolson method](figures/crank_nicolson/crank_nicolson_method.gif)

![HLL flux method](figures/hll/hll_flux_method.gif)

![Dirichlet method comparison](figures/boundary_conditions/dirichlet_method_comparison.gif)

![Neumann method comparison](figures/boundary_conditions/neumann_method_comparison.gif)

## How to run the code

1. Create and activate a Python virtual environment.
2. Install dependencies:

```bash
pip install -r requirements.txt
```

3. Run the implementation script that contains the numerical workflows and GIF generation:

```bash
python intentofinalpractica6.py
```

The repository keeps the uploaded source file in `src/` unchanged as requested.

## Repository structure

- `README.md` — project overview and usage.
- `docs/theory.md` — detailed theoretical chapter.
- `docs/numerical_method.md` — implementation workflow and method mapping.
- `docs/results_summary.md` — detailed interpretation of the six generated animations.
- `docs/sources_and_notes.md` — coursework provenance and curation notes.
- `intentofinalpractica6.py` — telegraph-equation implementation and animation pipeline.
- `src/damped_wave_telegraph_methods.py` — preserved uploaded Python file.
- `figures/` — six GIF outputs organized by method and boundary-condition comparison.

## Skills demonstrated

- PDE modeling for damped wave/telegraph dynamics.
- Spatial finite differences and CFL-aware explicit stepping.
- First-order system reformulation and matrix-based time integrators.
- Implicit linear solves with Crank–Nicolson.
- Conservative numerical fluxes with HLL/HLLE ideas.
- Boundary-condition modeling (Dirichlet vs Neumann) and comparative analysis.
- Scientific communication in reproducible portfolio format.

## Author

**mldc3** (GitHub repository owner and maintainer).

This repository is presented as a cleaned portfolio version preserving the existing numerical logic and assets.
