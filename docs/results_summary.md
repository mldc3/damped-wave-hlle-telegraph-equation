# Results Summary and Physical Interpretation

## 1. Purpose of the comparison

This project compares four numerical strategies applied to the same damped one-dimensional wave / telegraph equation:

$$
\frac{\partial^2 u}{\partial t^2}
-
c^2
\frac{\partial^2 u}{\partial x^2}
=
-2\frac{\kappa}{\rho}
\frac{\partial u}{\partial t}
+
a u.
$$

The purpose of the comparison is not simply to produce several animations of the same system. Each method represents a different numerical philosophy. The explicit finite-difference method follows directly from the discretized second-order equation. The matrix exponential method evolves the spatially discretized linear system almost exactly in time. Crank–Nicolson provides a stable implicit approximation to the same first-order system. The HLL/HLLE method rewrites the problem in conservative form and evolves it through numerical fluxes at cell interfaces.

Because the physical equation is the same in all cases, differences between the animations should be interpreted as numerical effects. These include artificial damping, phase error, boundary treatment and numerical viscosity.

---

## 2. Explicit finite-difference method

![Explicit finite-difference method](../figures/explicit/explicit_finite_difference.gif)

The explicit finite-difference method is the most direct way to approach the damped wave equation. It discretizes the second derivative in time, the second derivative in space and the damping term, then solves algebraically for the next time level.

The method is pedagogically valuable because every term in the update formula can be traced back to the original PDE. The future value $u_j^{n+1}$ is computed from neighbouring values at the current time level and values from the previous time level. This makes the method intuitive and computationally cheap per step.

However, the explicit method is limited by the CFL condition,

$$
\frac{c\Delta t}{\Delta x}\leq 1.
$$

Even when this condition is satisfied, the method can introduce numerical damping and phase error. Therefore, in the animation, the qualitative behaviour should be correct: the wave evolves, respects the imposed boundary treatment and decays due to damping. But if the wave appears to lose amplitude more strongly than in the matrix exponential or Crank–Nicolson animations, that extra decay should be interpreted as numerical dissipation rather than new physics.

The explicit method is therefore a good baseline: simple, transparent and fast, but not necessarily the most accurate for long-time damped wave evolution.

---

## 3. Matrix exponential method

![Matrix exponential method](../figures/matrix_exponential/matrix_exponential_method.gif)

The matrix exponential method starts by rewriting the second-order PDE as a first-order system,

$$
\frac{d\mathbf{W}}{dt}
=
\mathbf{A}\mathbf{W}.
$$

After spatial discretization, the exact time evolution of this semi-discrete linear system is

$$
\mathbf{W}^{n+1}
=
e^{\mathbf{A}\Delta t}
\mathbf{W}^{n}.
$$

This makes the matrix exponential method a strong reference for the other schemes. It does not remove spatial discretization error, but it treats the time evolution of the discretized linear system with very high fidelity.

In the animation, the decay of the wave should mainly reflect the physical damping present in the equation. Compared with the explicit method, the matrix exponential solution is expected to show cleaner phase behaviour and less artificial numerical dissipation.

This method is especially useful for validation. If another method produces a visibly different decay rate or phase evolution, the matrix exponential animation helps identify whether the difference is likely due to the time integrator rather than the physical model.

Its main drawback is computational cost. Computing a matrix exponential can be expensive, so although it is excellent as a reference, it may not be the best practical method for very large simulations.

---

## 4. Crank–Nicolson method

![Crank–Nicolson method](../figures/crank_nicolson/crank_nicolson_method.gif)

Crank–Nicolson is applied to the first-order matrix system through the update

$$
\left(
\mathbf{I}
-
\frac{\Delta t}{2}\mathbf{A}
\right)
\mathbf{W}^{n+1}
=
\left(
\mathbf{I}
+
\frac{\Delta t}{2}\mathbf{A}
\right)
\mathbf{W}^{n}.
$$

This method can be interpreted as a rational approximation to the matrix exponential. It is second order in time and generally much more stable than a simple explicit method.

The animation should be compared closely with the matrix exponential result. If both methods show similar amplitude decay and similar phase evolution, this is strong evidence that Crank–Nicolson is capturing the correct damped wave dynamics.

Crank–Nicolson is often the best practical compromise in this project. It is more expensive per step than the explicit method because it requires solving a linear system, but it is usually more accurate and stable. It is less exact than the matrix exponential method, but more scalable and more practical.

---

## 5. HLL/HLLE flux method

![HLL flux method](../figures/hll/hll_flux_method.gif)

The HLL/HLLE method approaches the same physical equation from a conservative flux perspective. Instead of directly updating the second-order equation or applying a matrix time propagator, the method writes the system in a form resembling

$$
\frac{\partial \mathbf{U}}{\partial t}
+
\frac{\partial \mathbf{F}(\mathbf{U})}{\partial x}
=
\mathbf{S}(\mathbf{U}).
$$

The update is based on numerical fluxes at cell interfaces. At each interface, the method uses a left state and a right state and estimates the fastest waves moving left and right. The HLL flux is

$$
\mathbf{F}_{HLL}
=
\frac{
S_R\mathbf{F}_L
-
S_L\mathbf{F}_R
+
S_LS_R
(\mathbf{U}_R-\mathbf{U}_L)
}{
S_R-S_L
}.
$$

The important point is that HLL/HLLE includes numerical viscosity. This viscosity stabilizes the method and makes it robust for hyperbolic systems, especially when discontinuities or sharp gradients are present. In a smooth linear damped wave problem, however, that same viscosity may make the solution look more diffusive than the matrix exponential or Crank–Nicolson methods.

Therefore, if the HLL animation appears smoother or more strongly damped, this should not be interpreted as a mistake. It is a signature of the method’s conservative and robust flux formulation.

---

## 6. Dirichlet method comparison

![Dirichlet method comparison](../figures/boundary_conditions/dirichlet_method_comparison.gif)

The Dirichlet comparison shows the behaviour of the numerical methods under fixed-end boundary conditions,

$$
u(0,t)=0,
\qquad
u(L,t)=0.
$$

These conditions force the field to remain zero at the endpoints. Physically, this resembles a vibrating string fixed at both ends. The allowed modes must be compatible with the endpoints, and reflections occur according to the fixed-boundary constraint.

This comparison is important because it places all methods under the same boundary condition. Differences between the curves are therefore not due to different physical endpoints, but to the numerical method itself.

The matrix exponential and Crank–Nicolson methods are expected to agree closely. The explicit method may show more damping or slightly different phase behaviour. The HLL method may appear smoother because of numerical viscosity. This single comparison therefore summarizes the central lesson of the repository: the same PDE and the same boundary condition can produce visibly different numerical behaviour depending on the algorithm.

---

## 7. Neumann method comparison

![Neumann method comparison](../figures/boundary_conditions/neumann_method_comparison.gif)

The Neumann comparison uses derivative boundary conditions, typically of the form

$$
\frac{\partial u}{\partial x}=0
$$

at the boundary. This represents a zero-gradient or free-end condition. It is physically different from Dirichlet because the endpoint is not forced to remain at zero displacement.

The animation should therefore be interpreted as a different physical problem, not just as a different implementation detail. Neumann conditions change the allowed modes and the reflection behaviour at the boundary.

This comparison is especially useful because derivative boundary conditions are often more delicate to implement than fixed-value conditions. They may require ghost-point relations, one-sided finite differences or modified boundary rows in matrix methods. If all methods produce compatible behaviour under Neumann conditions, this supports the consistency of the boundary implementation.

---

## 8. Main interpretation

Across the six animations, the central result is the comparison between numerical philosophies.

The explicit method is simple and direct, but it is limited by stability and may introduce artificial damping.

The matrix exponential method is the cleanest reference for the semi-discrete linear system because it evolves the first-order system through $e^{\mathbf{A}\Delta t}$.

Crank–Nicolson is a strong practical compromise. It approximates the matrix exponential well while avoiding the full cost of computing an exact propagator.

HLL/HLLE is robust and conservative. Its numerical viscosity is useful for hyperbolic systems, but in this smooth damped wave problem it can make the solution look more diffusive.

The Dirichlet and Neumann comparisons show that boundary conditions are part of the physical model. They change the modal structure and the way waves interact with the domain endpoints.

---

## 9. Conclusions

This project demonstrates that numerical solution of a PDE is not only a matter of coding the equation. The choice of numerical method changes the observed behaviour through stability limits, phase error, artificial damping, solver structure and flux viscosity.

The practice shows competence in:
- damped wave / telegraph equations,
- finite-difference discretization,
- first-order system reformulation,
- matrix exponential propagation,
- Crank–Nicolson integration,
- conservative flux methods,
- HLL/HLLE approximate Riemann solvers,
- Dirichlet and Neumann boundary conditions,
- physical damping versus numerical damping,
- scientific Python programming and animated visualization.

The six GIFs provide a compact but meaningful comparison of the main numerical approaches.
