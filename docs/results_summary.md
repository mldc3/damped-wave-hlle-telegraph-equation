# Results Summary and Interpretation

This analysis is based strictly on the six uploaded GIF files listed below.

![Explicit finite-difference method](../figures/explicit/explicit_finite_difference.gif)

![Matrix exponential method](../figures/matrix_exponential/matrix_exponential_method.gif)

![Crank–Nicolson method](../figures/crank_nicolson/crank_nicolson_method.gif)

![HLL flux method](../figures/hll/hll_flux_method.gif)

![Dirichlet method comparison](../figures/boundary_conditions/dirichlet_method_comparison.gif)

![Neumann method comparison](../figures/boundary_conditions/neumann_method_comparison.gif)

## 1. Purpose of comparing four numerical methods

The objective is not only to produce stable numerical runs, but to understand how different discretization philosophies represent the same damped wave physics. The explicit method, matrix exponential method, Crank–Nicolson method, and HLL/HLLE conservative method each encode different compromises between cost, fidelity, and robustness.

Using common parameters and matched boundary-condition scenarios allows direct interpretation of phase behavior, attenuation rate, waveform shape preservation, and sensitivity to numerical viscosity.

## 2. Detailed interpretation of the explicit finite-difference GIF

In the explicit animation, oscillatory transport is visible together with progressive attenuation consistent with damping in the governing equation. The method captures propagation clearly, but compared with higher-fidelity references it can display slightly stronger smoothing and phase drift over long integration windows, especially when the CFL ratio is chosen for practicality rather than strict high-accuracy tracking.

The Dirichlet-versus-Neumann contrast in the same animation makes boundary influence easy to observe: fixed-end constraints enforce stronger endpoint anchoring, while free-end style behavior on the Neumann side alters reflection patterns and endpoint amplitude response.

## 3. Detailed interpretation of the matrix exponential GIF

The matrix exponential animation is the cleanest representation of the semi-discrete linear dynamics because temporal evolution follows $e^{\mathbf{A}\Delta t}$ for the built operator. This generally yields smooth phase progression and attenuation behavior that can be interpreted as a reference trajectory for the same spatial discretization.

Because the method is exact in time for the linear semi-discrete system, differences from other methods are informative: deviations seen elsewhere are mostly attributable to temporal approximation error or additional method-specific dissipation.

## 4. Detailed interpretation of the Crank–Nicolson GIF

The Crank–Nicolson animation closely follows matrix-exponential behavior in both phase and amplitude envelope. This is expected from its second-order implicit midpoint structure and its rational approximation relationship to the matrix exponential operator.

Across the run, Crank–Nicolson generally retains oscillatory structure better than strongly dissipative alternatives while remaining computationally more practical than repeated full matrix exponentials in larger settings.

## 5. Detailed interpretation of the HLL/HLLE GIF

The HLL/HLLE animation emphasizes robustness and conservative interface treatment. A typical visual signature is additional smoothing relative to matrix exponential and Crank–Nicolson, especially in regions where gradients develop.

This extra diffusion is consistent with numerical viscosity introduced by approximate Riemann flux construction. Even in smooth regimes, the stabilizing jump term can reduce amplitude and sharpen damping appearance beyond the physical damping term alone.

## 6. Detailed comparison under Dirichlet boundary conditions

Under Dirichlet constraints, all methods satisfy fixed-end behavior, but interior evolution differs.

The comparison GIF shows that matrix exponential and Crank–Nicolson tend to remain closely aligned over time, indicating similar phase and damping representation for this linear problem. The explicit curve remains physically consistent but can separate gradually in amplitude or phase due to timestep-discretization effects. The HLL/HLLE solution typically appears more damped and smoother, reflecting stronger numerical viscosity.

Dirichlet constraints also make reflection timing and amplitude near the boundaries especially diagnostic for method fidelity.

## 7. Detailed comparison under Neumann boundary conditions

In the Neumann comparison, endpoint-gradient handling changes reflection structure and can preserve larger endpoint displacement in the free-end direction compared with fixed-end cases. This boundary change highlights how methods respond when closure relies on derivative constraints rather than fixed values.

The same ordering trend is visible: matrix exponential and Crank–Nicolson show close agreement; explicit remains credible with somewhat larger long-time discretization effects; HLL/HLLE appears systematically more diffusive.

## 8. Physical damping versus numerical damping

All methods include the same physical damping mechanism from the PDE parameters, so shared amplitude decay trends indicate real model attenuation. Method-dependent spread around that shared trend indicates numerical damping differences.

When one method consistently attenuates faster under the same physical setup, that excess attenuation is interpreted as numerical viscosity rather than physical dissipation.

## 9. Why HLL/HLLE may appear more diffusive

HLL/HLLE uses bounded wave speeds and a single approximate interface flux formula to remain robust without full characteristic reconstruction. The stabilizing structure introduces dissipative terms proportional to state jumps across interfaces.

In smooth damped-wave simulations, this mechanism can still produce visible smearing and stronger envelope decay than less dissipative linear-system solvers. The tradeoff is robustness and conservative flux consistency.

## 10. Why matrix exponential and Crank–Nicolson should agree closely

Both methods evolve the same first-order semi-discrete linear system. Matrix exponential is exact in time for that ODE, while Crank–Nicolson is a second-order implicit approximation with strong phase properties for linear oscillatory dynamics.

Therefore close alignment in the GIFs is expected and acts as an internal consistency check for implementation quality.

## 11. Final conclusions and skills demonstrated

The six uploaded animations provide a coherent comparative study of the damped telegraph equation under two boundary-condition families. They show that method choice affects not only stability but also waveform fidelity and apparent damping.

Main conclusions:

- Matrix exponential provides a high-fidelity reference trajectory.
- Crank–Nicolson offers near-reference behavior with practical implicit stepping.
- Explicit finite differences remain effective but are more sensitive to timestep and dispersion/dissipation balance.
- HLL/HLLE is robust and conservative, with a characteristic increase in numerical diffusion.

Skills demonstrated include PDE discretization, boundary-condition modeling, implicit and explicit integration, conservative flux methods, and rigorous interpretation of physical versus numerical dissipation in computational physics workflows.
