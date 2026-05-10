# Damped Wave / Telegraph Equation: Explicit, Crank–Nicolson and HLL Solvers

This repository contains a cleaned portfolio version of a computational physics project focused on Practice 6 and the damped one-dimensional wave / telegraph equation.

The project includes:
- explicit finite-difference integration,
- matrix exponential time evolution,
- Crank–Nicolson time evolution,
- HLL/HLLE conservative flux methods,
- Dirichlet and Neumann boundary-condition comparisons,
- physical damping versus numerical damping,
- method comparison through animations.

## Repository structure
- `src/` - Uploaded Python source file for the project portfolio.
- `docs/` - Theory, numerical-method notes, result summaries, and source notes.
- `figures/` - Animations organized by method and boundary-condition comparisons.
- `notes/original_report/` - Original report archive location.

## Animations
### Explicit method figures (`figures/explicit/`)
![Explicit method](figures/explicit/explicit_finite_difference.gif)

### Matrix exponential figures (`figures/matrix_exponential/`)
![Matrix exponential method](figures/matrix_exponential/matrix_exponential_method.gif)

### Crank–Nicolson figures (`figures/crank_nicolson/`)
![Crank–Nicolson method](figures/crank_nicolson/crank_nicolson_method.gif)

### HLL/HLLE figures (`figures/hll/`)
![HLL flux method](figures/hll/hll_flux_method.gif)

### Boundary comparisons (`figures/boundary_conditions/`)
![Dirichlet method comparison](figures/boundary_conditions/dirichlet_method_comparison.gif)
![Neumann method comparison](figures/boundary_conditions/neumann_method_comparison.gif)

### Method comparisons (`figures/comparison/`)
No figure files are currently present in `figures/comparison/`.

### CFL / stability / damping plots (`figures/stability/`)
No figure files are currently present in `figures/stability/`.

### Runtime plots (`figures/performance/`)
No figure files are currently present in `figures/performance/`.

## Requirements
Install dependencies with:

```bash
pip install -r requirements.txt
```
