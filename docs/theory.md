# Theory

This project studies the damped one-dimensional wave / telegraph equation in the context of computational physics.

The equation combines wave propagation with dissipation. In this setting, the damping term reduces amplitude over time while preserving the underlying transport behavior.

The repository is limited to Practice 6 scope:
- damped 1D wave / telegraph dynamics,
- boundary conditions (Dirichlet and Neumann),
- comparison of numerical behavior across multiple time-integration and flux formulations.

## 18. Comparison of the numerical philosophies

The project compares four distinct numerical ideas.

The explicit finite-difference scheme is the most direct. It discretizes the second-order equation and advances the wave using known previous states. It is easy to implement and fast per step, but it is limited by the CFL condition and can introduce artificial damping or phase error.

The matrix exponential method reformulates the system as a first-order linear ODE after spatial discretization. It then uses the exact propagator

$$
e^{\mathbf{A}\Delta t}.
$$

This makes it an excellent reference for the semi-discrete problem. The price is the computational cost of forming and applying the matrix exponential.

Crank–Nicolson also uses the first-order system, but replaces the exact exponential with a second-order implicit rational approximation. It is stable, accurate and usually close to the matrix exponential solution, while being more practical for many simulations.

HLL/HLLE uses a conservative finite-volume perspective. It evolves cell averages through interface fluxes and uses approximate Riemann problems. It is robust and physically natural for hyperbolic conservation laws, but introduces numerical viscosity that can be visible in smooth wave problems.

## 19. Expected behaviour of the solution

For a low sinusoidal mode and zero initial velocity, the wave starts from a smooth standing-wave-like configuration. Because of the damping term, the amplitude should gradually decay. Boundary conditions determine how the wave interacts with the domain endpoints. Dirichlet boundaries force the field to remain zero at the ends, while Neumann boundaries impose zero slope and therefore modify the allowed modal structure.

The methods should agree qualitatively:
- the wave should remain bounded,
- the amplitude should decrease,
- propagation should occur at approximately speed $c$,
- boundary conditions should be respected,
- higher modes should be more sensitive to numerical resolution.

They should differ quantitatively:
- explicit finite differences may be more dissipative,
- matrix exponential should be closest to exact time evolution for the discretized operator,
- Crank–Nicolson should closely follow matrix exponential,
- HLL/HLLE should be smoother and more damped because of numerical viscosity.

## 20. Why this practice is important

This practice is not only about one particular equation. It connects several major topics in numerical physics.

It shows how second-order PDEs can be turned into first-order systems. It shows why matrix formulations are useful. It shows how explicit and implicit methods differ in stability and accuracy. It introduces the idea that hyperbolic systems should often be treated through characteristic speeds and numerical fluxes. It also demonstrates that conservation, stability and accuracy are separate concepts.

A method may be conservative but diffusive. A method may be stable but inaccurate. A method may be accurate but expensive. A good computational physicist must understand these trade-offs instead of choosing algorithms mechanically.

The main conceptual lesson is that numerical methods are part of the model. They do not merely “solve” the equation; they impose assumptions about propagation, dissipation, stability and boundary behaviour. The purpose of this repository is to make those assumptions visible through code, figures and method comparison.
