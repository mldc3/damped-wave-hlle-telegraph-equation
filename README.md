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
### Explicit finite difference
![Explicit finite difference](figures/explicit/explicit_finite_difference.gif)

### Matrix exponential method
![Matrix exponential method](figures/matrix_exponential/matrix_exponential_method.gif)

### Crank–Nicolson method
![Crank–Nicolson method](figures/crank_nicolson/crank_nicolson_method.gif)

### HLL flux method
![HLL flux method](figures/hll/hll_flux_method.gif)

### Dirichlet boundary-condition comparison
![Dirichlet method comparison](figures/boundary_conditions/dirichlet_method_comparison.gif)

### Neumann boundary-condition comparison
![Neumann method comparison](figures/boundary_conditions/neumann_method_comparison.gif)

## Requirements
Install dependencies with:

```bash
pip install -r requirements.txt
```
