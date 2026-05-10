# Theory: The Damped Wave / Telegraph Equation

## 1. Physical motivation for the damped wave / telegraph equation

Many physical systems propagate disturbances as waves but also dissipate energy while they propagate. Ideal strings, ideal acoustic media, and ideal electromagnetic wave models are useful abstractions, yet real materials and environments include frictional, resistive, or viscous effects. The damped wave (telegraph) equation is a compact model that captures this mixed behavior: finite-speed wave transport together with amplitude decay.

Historically, the telegraph equation appears in electrical transmission-line modeling, where distributed inductance and capacitance support wave propagation while resistance and conductance attenuate signals. In mechanical language, the same mathematical structure appears when displacement fields evolve under restoring forces and damping. The model is therefore a bridge between pure oscillatory behavior and diffusive attenuation.

Because it retains second-order time dynamics, the equation keeps an inertial wave character that first-order diffusion models do not capture. This is exactly why it is a good testbed for comparing explicit schemes, implicit midpoint-type schemes, and conservative flux-based formulations.

## 2. Relation with the ideal wave equation

The ideal one-dimensional wave equation is

$$
\frac{\partial^2 u}{\partial t^2}-c^2\frac{\partial^2 u}{\partial x^2}=0.
$$

Its key property is nondissipative propagation: in the ideal continuous model, modal energy is redistributed in space but not systematically destroyed. By contrast, the damped/telegraph form introduces linear terms in $u_t$ and $u$:

$$
\frac{\partial^2 u}{\partial t^2}-c^2\frac{\partial^2 u}{\partial x^2}
=-2\frac{\kappa}{\rho}\frac{\partial u}{\partial t}+a u.
$$

If $\kappa=0$ and $a=0$, the ideal wave equation is recovered. If damping is present, oscillations still travel but their amplitude envelope decays. The equation is therefore not a different physical category; it is an enriched wave model with attenuation and linear reaction/restoring effects.

## 3. Meaning of each parameter

For

$$
\frac{\partial^2 u}{\partial t^2}
-
c^2
\frac{\partial^2 u}{\partial x^2}
=
-2\frac{\kappa}{\rho}
\frac{\partial u}{\partial t}
+
a u,
$$

- $u(x,t)$ is the field variable (for example displacement or voltage-like quantity).
- $c$ is the wave speed in the ideal limit, controlling characteristic propagation velocity.
- $\kappa/\rho$ is the effective damping ratio in the first-order time derivative term.
- $a$ is a linear source/restoring coefficient that shifts modal dynamics.
- $x$ is the spatial coordinate and $t$ is time.

The sign and magnitude of $a$ matter. Depending on convention and regime, it may reinforce restoration or alter effective modal frequency. In this project, parameters are fixed in code and are not altered by documentation changes.

## 4. Physical damping versus numerical damping

Physical damping is part of the PDE itself. It is independent of mesh resolution and should remain when the numerical grid is refined. If resolution changes and damping behavior converges toward a common trajectory, that behavior is likely physical.

Numerical damping is introduced by discretization choices. Upwinding, approximate Riemann solvers, low-order temporal schemes, and coarse grids can remove amplitude even when the continuous equation would not. This artificial attenuation can be useful for robustness, but it can also distort phase and amplitude if excessive.

A core objective of method comparison is therefore to separate attenuation caused by $-2(\kappa/\rho)u_t$ from attenuation caused by scheme-dependent viscosity.

## 5. Need for two initial conditions

Because the equation is second order in time, one must prescribe both displacement and velocity information at $t=0$. In continuous terms, this means specifying:

$$
u(x,0)=u_0(x), \qquad \frac{\partial u}{\partial t}(x,0)=v_0(x).
$$

In explicit second-order finite differences, this is often implemented as $u^0$ and $u^1$ (or equivalently $u^0$ plus $u_t^0$ and a startup relation). Without this second piece of data, the temporal recurrence is underdetermined.

## 6. Sinusoidal initial modes and the meaning of mode number

A standard initialization is modal:

$$
u(x,0)=\sin(n\pi x),
$$

with integer mode number $n$. The mode number counts spatial oscillations; larger $n$ means shorter wavelength and larger curvature, so high modes are generally more sensitive to spatial resolution and numerical dissipation.

For boundary-value wave problems, modal initial data provide clear interpretation: each mode has a known shape, and the time evolution reveals damping, phase shift, and boundary reflections in a transparent way.

## 7. Dirichlet and Neumann boundary conditions

Boundary conditions close the PDE and strongly influence dynamics.

- Dirichlet boundaries prescribe field values, typically fixed ends such as $u=0$.
- Neumann boundaries prescribe spatial derivatives, such as $u_x=0$, corresponding to zero-flux or free-end behavior.

Physically, a fixed end imposes displacement constraints and tends to reflect with sign change depending on context. A free end imposes slope constraints and reflects differently. Comparing both in the same numerical framework clarifies how boundary physics and numerical treatment interact.

## 8. Finite-difference grid and notation

Let $x_j=x_{\min}+j\Delta x$, $j=0,\dots,N_x-1$, and $t^n=n\Delta t$, $n=0,\dots,N_t-1$. Denote the discrete solution by $u_j^n\approx u(x_j,t^n)$.

The centered second derivative in space is

$$
\frac{\partial^2 u}{\partial x^2}(x_j,t^n)
\approx
\frac{u_{j+1}^n-2u_j^n+u_{j-1}^n}{\Delta x^2}.
$$

The second derivative in time for explicit schemes is approximated with a three-level stencil involving $u_j^{n+1}$, $u_j^n$, and $u_j^{n-1}$.

## 9. Explicit finite-difference scheme and CFL condition

A common update for interior points is obtained by discretizing the PDE and isolating $u_j^{n+1}$. The recurrence combines:

- the two previous time levels,
- the discrete Laplacian,
- damping contribution,
- linear source term.

A dimensionless control quantity is

$$
r=\left(\frac{c\Delta t}{\Delta x}\right)^2.
$$

For wave-like explicit methods, stability requires a CFL-type restriction linking $\Delta t$ to $\Delta x$ and $c$. The precise limit depends on full discretization details (including damping and source terms), but the principle is universal: increasing temporal step size beyond the admissible range causes oscillatory blow-up or severe nonphysical behavior.

## 10. Why stability is not the same as accuracy

A stable method only guarantees bounded growth under a given norm and timestep condition. It does not guarantee correct phase speed, amplitude decay rate, or dispersion properties.

One can be stable and still inaccurate if:

- the grid is too coarse to resolve wavelengths,
- the timestep is too large for phase fidelity,
- numerical viscosity dominates physical damping.

Hence meaningful comparison must examine qualitative and quantitative signal behavior, not merely whether the run avoids divergence.

## 11. Reformulation as a first-order system using $v=u_t$

Define

$$
v=\frac{\partial u}{\partial t}.
$$

Then

$$
\frac{\partial u}{\partial t}=v,
$$

and

$$
\frac{\partial v}{\partial t}=c^2\frac{\partial^2 u}{\partial x^2}-2\frac{\kappa}{\rho}v+a u.
$$

This conversion turns one second-order equation into coupled first-order equations, which is convenient for matrix evolution and implicit integrators.

## 12. Matrix formulation with state vector $\mathbf{W}=(u,v)^T$

After spatial discretization, collect grid values into vectors $\mathbf{u}$ and $\mathbf{v}$. The semi-discrete state

$$
\mathbf{W}=\begin{pmatrix}\mathbf{u}\\ \mathbf{v}\end{pmatrix}
$$

satisfies a linear ODE

$$
\frac{d\mathbf{W}}{dt}=\mathbf{A}\mathbf{W},
$$

with block structure

$$
\mathbf{A}=
\begin{pmatrix}
\mathbf{0} & \mathbf{I}\\
c^2\mathbf{L}+a\mathbf{I} & -2(\kappa/\rho)\mathbf{I}
\end{pmatrix},
$$

where $\mathbf{L}$ is the discrete Laplacian (modified at boundaries according to Dirichlet or Neumann treatment).

## 13. Matrix exponential as high-fidelity reference

For constant $\mathbf{A}$, the exact time map of the semi-discrete linear system is

$$
\mathbf{W}^{n+1}=e^{\mathbf{A}\Delta t}\mathbf{W}^n.
$$

This does not approximate time integration at the ODE level; it is exact for the chosen spatial discretization and boundary-imposed operator. Therefore it is an excellent reference when assessing phase and damping behavior of other integrators applied to the same semi-discrete model.

Its drawback is cost and dense linear-algebra overhead as system size grows.

## 14. Crank–Nicolson as second-order implicit approximation

Crank–Nicolson applied to $d\mathbf{W}/dt=\mathbf{A}\mathbf{W}$ yields

$$
\left(\mathbf{I}-\frac{\Delta t}{2}\mathbf{A}\right)\mathbf{W}^{n+1}
=
\left(\mathbf{I}+\frac{\Delta t}{2}\mathbf{A}\right)\mathbf{W}^{n}.
$$

This is second order in time and can be interpreted as a rational approximation to $e^{\mathbf{A}\Delta t}$:

$$
e^{\mathbf{A}\Delta t}
\approx
\left(\mathbf{I}-\frac{\Delta t}{2}\mathbf{A}\right)^{-1}
\left(\mathbf{I}+\frac{\Delta t}{2}\mathbf{A}\right).
$$

That interpretation explains why Crank–Nicolson often tracks matrix-exponential results closely while remaining computationally more scalable.

## 15. Conservative formulation with $\mathbf{U}=(u,v,w)^T$

Introduce

$$
w=\frac{\partial u}{\partial x}.
$$

One can write a first-order system in conservative-like form

$$
\frac{\partial \mathbf{U}}{\partial t}+\frac{\partial \mathbf{F}(\mathbf{U})}{\partial x}=\mathbf{S}(\mathbf{U}),
$$

with representative state components displacement, temporal velocity, and spatial gradient. This representation enables finite-volume style updates and interface fluxes.

## 16. Fluxes, Jacobian, characteristic speeds, and wave directions

In the implemented HLL-style system, the physical flux has components analogous to

$$
\mathbf{F}(u,v,w)=\begin{pmatrix}0\\ -c^2 w\\ -v\end{pmatrix}.
$$

The flux Jacobian has characteristic speeds approximately $\lambda_\pm=\pm c$ for propagating modes. Positive speed corresponds to right-moving information; negative speed corresponds to left-moving information. The Riemann solver uses these wave-speed bounds to build a single consistent interface flux.

## 17. HLL/HLLE approximate Riemann solver

Given left and right interface states, HLL/HLLE replaces the full wave fan by two bounding speeds $s_-$ and $s_+$. The interface flux is then a piecewise formula that blends left and right physical fluxes plus a jump term proportional to state difference. This avoids full characteristic decomposition while preserving robust shock/gradient-capturing behavior in hyperbolic-like systems.

For smooth damped-wave dynamics, it functions as a robust conservative transport update with source coupling.

## 18. Numerical viscosity in HLL/HLLE

The same jump-stabilization mechanism that makes HLL/HLLE robust also introduces numerical viscosity. In smooth oscillatory regimes, this often manifests as stronger amplitude smearing relative to matrix-exponential or Crank–Nicolson trajectories.

This is not a bug; it is a known tradeoff in approximate Riemann methods. The method can be preferable when robustness and conservative interface treatment are priorities.

## 19. Detailed comparison of four numerical philosophies

The explicit second-order finite-difference method is direct and inexpensive but constrained by CFL and more exposed to timestep sensitivity. Matrix exponential is the most faithful in time for the semi-discrete linear operator but computationally heavier. Crank–Nicolson balances fidelity and efficiency through an implicit linear solve each step. HLL/HLLE reframes the PDE in conservative form with interface fluxes and tends to be more diffusive but very robust.

Under identical physical parameters, differences in attenuation and waveform sharpness primarily reveal each scheme’s numerical dissipation and dispersion signatures. Agreement between matrix exponential and Crank–Nicolson supports the expected second-order implicit approximation quality. Greater smoothing in HLL/HLLE indicates stronger numerical viscosity beyond physical damping.

## 20. Final theoretical interpretation

The damped telegraph equation is a wave system with memory, inertia, and dissipation. It requires careful handling of initial data, boundary physics, and discretization choices. The four methods studied are not merely alternative code paths; they represent distinct numerical philosophies:

- direct explicit recurrence,
- operator-exact linear evolution,
- implicit midpoint-style approximation,
- conservative interface-flux transport.

A rigorous interpretation of results must always separate what belongs to the physics of the PDE from what belongs to the numerics of approximation. This distinction is the central theoretical lesson of the project and the reason method comparison is essential in computational physics.

<!-- THEORY PART 0 END -->

---

## 1. Physical motivation

The starting point of the project is the damped one-dimensional wave equation

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

The unknown field is $u(x,t)$. It may be interpreted as a displacement-like quantity evolving along a one-dimensional domain. The parameter $c$ is the wave speed, so it controls how fast disturbances propagate. The ratio $\kappa/\rho$ controls physical damping, and the coefficient $a$ introduces a linear restoring or reaction term proportional to the field itself.

The left-hand side is the standard wave operator. If the right-hand side were zero, the equation would reduce to the ideal wave equation,

$$
\frac{\partial^2 u}{\partial t^2}
-
c^2
\frac{\partial^2 u}{\partial x^2}
=
0.
$$

This equation describes waves travelling with finite speed. A disturbance introduced at one point does not instantly affect the whole domain. Instead, information propagates along characteristic directions associated with speeds $+c$ and $-c$.

The right-hand side modifies the ideal wave equation in two ways. The term

$$
-2\frac{\kappa}{\rho}
\frac{\partial u}{\partial t}
$$

is a damping term. It acts against motion and removes energy from the wave. As a result, even if the initial condition is a clean sinusoidal mode, the amplitude should decay with time. This decay is physical because it is explicitly included in the differential equation.

The term

$$
a u
$$

is a linear contribution proportional to the field. Depending on the sign of $a$, it may act like an additional restoring term or like a source of growth. In the practice, the representative value is negative, so it contributes to bounded oscillatory behaviour rather than uncontrolled exponential growth.

This equation is useful because it contains several ingredients that are central in computational physics: wave propagation, damping, boundary conditions, stability restrictions, matrix formulations, implicit integration and conservative flux methods.

---

## 2. Why this is a hyperbolic problem

The ideal wave equation is a hyperbolic partial differential equation. Hyperbolic problems are characterized by finite-speed propagation and characteristic directions. This makes them very different from parabolic equations such as the heat equation.

For the heat equation, a perturbation spreads diffusively and smooths out. For a wave equation, a perturbation travels. This means that numerical methods must respect the direction and speed at which information propagates.

The damped telegraph equation keeps the hyperbolic character of the wave equation but adds decay. The wave still propagates through the domain, but its amplitude is reduced by physical damping. A good numerical method should reproduce both aspects: propagation and damping.

This is why it is important to separate physical damping from numerical damping. Physical damping comes from the equation itself. Numerical damping comes from the method used to approximate the equation. A method can be stable and still too dissipative if it damps the wave more strongly than the physical equation requires.

---

## 3. Initial conditions

Because the equation is second order in time, one initial condition is not enough. The future evolution is determined by both the initial field and the initial velocity:

$$
u(x,0)=u_0(x),
$$

$$
\frac{\partial u}{\partial t}(x,0)=v_0(x).
$$

This is analogous to classical mechanics. To determine the motion of a particle governed by a second-order equation, one must specify both position and velocity. Similarly, for a second-order wave equation, one must specify both the initial shape and the initial rate of change.

A common initial condition in this practice is a sinusoidal mode,

$$
u(x,0)=\sin\left(\frac{n\pi x}{L}\right),
$$

with zero initial velocity,

$$
\frac{\partial u}{\partial t}(x,0)=0.
$$

The integer $n$ is the mode number. Low values of $n$ correspond to long wavelengths and smooth spatial profiles. Higher values of $n$ correspond to shorter wavelengths and more oscillations over the same domain. High modes are harder to resolve numerically because they require a finer spatial grid.

This is important when comparing methods. A method may perform well for a smooth low mode but introduce stronger damping or phase error for high-frequency components.

<!-- THEORY PART 1 END -->

---

## 4. Boundary conditions

Boundary conditions are part of the physical model. They are not simply a technical coding detail. They determine which solutions are allowed and how waves interact with the endpoints of the domain.

For Dirichlet boundary conditions, the value of the field is fixed at the boundaries:

$$
u(0,t)=0,
\qquad
u(L,t)=0.
$$

Physically, this resembles a string fixed at both ends. The endpoints cannot move. A wave reaching the boundary reflects according to this fixed-end constraint.

For Neumann boundary conditions, the spatial derivative is prescribed. A homogeneous Neumann boundary condition has the form

$$
\frac{\partial u}{\partial x}=0
$$

at the boundary. This represents a zero-gradient or free-end condition. The endpoint is not forced to have zero displacement; instead, the slope is constrained.

The difference between Dirichlet and Neumann conditions is physically visible. They lead to different allowed modes and different reflection behaviour. Therefore, the comparison of Dirichlet and Neumann animations is not merely a numerical test. It is also a comparison between two different physical boundary models.

---

## 5. Discrete grid

To solve the equation numerically, the continuous domain is replaced by a finite grid. Let

$$
x_j=x_{\min}+j\Delta x,
\qquad
j=0,1,\ldots,N-1,
$$

and

$$
t^n=n\Delta t.
$$

The numerical solution is denoted by

$$
u_j^n \approx u(x_j,t^n).
$$

The goal of a finite-difference method is to replace derivatives by algebraic combinations of nearby grid values. This transforms the partial differential equation into a recurrence relation or a matrix system.

The second derivative in time is approximated by

$$
\frac{\partial^2 u}{\partial t^2}
\approx
\frac{u_j^{n+1}-2u_j^n+u_j^{n-1}}{\Delta t^2}.
$$

The second derivative in space is approximated by

$$
\frac{\partial^2 u}{\partial x^2}
\approx
\frac{u_{j+1}^{n}-2u_j^n+u_{j-1}^{n}}{\Delta x^2}.
$$

The damping term can be approximated using a centred time difference:

$$
\frac{\partial u}{\partial t}
\approx
\frac{u_j^{n+1}-u_j^{n-1}}{2\Delta t}.
$$

These approximations are natural because they are centred and second-order accurate in the variables they approximate. They also make clear why the method uses neighbouring spatial points and two previous time levels.

---

## 6. Explicit finite-difference update

Substituting the finite differences into the damped wave equation gives a direct update for $u_j^{n+1}$. Define

$$
k=\frac{\kappa}{\rho},
$$

and

$$
r=\left(\frac{c\Delta t}{\Delta x}\right)^2.
$$

A representative explicit update has the form

$$
u_j^{n+1}
=
\frac{
(2-k\Delta t+a\Delta t^2)u_j^n
-
(1-k\Delta t)u_j^{n-1}
+
r
\left(
u_{j+1}^{n}
-
2u_j^{n}
+
u_{j-1}^{n}
\right)
}{
1+k\Delta t
}.
$$

The precise arrangement of signs depends on the convention used to move terms from one side of the equation to the other, but the structure is the same: the future value is computed from known values at the current and previous time levels.

This explicit scheme is attractive because it is simple and cheap. No linear system has to be solved. Each time step is obtained by applying a local stencil.

However, explicit wave schemes are conditionally stable. The time step must be small enough relative to the spatial step. The relevant quantity is the Courant number,

$$
\frac{c\Delta t}{\Delta x}.
$$

For the standard centred wave scheme, the stability condition is essentially

$$
\frac{c\Delta t}{\Delta x}\leq 1.
$$

This is the CFL condition. Its physical meaning is that information cannot be allowed to move across more grid cells in one time step than the numerical stencil can represent. If the time step is too large, the scheme may become unstable.

---

## 7. Stability is not accuracy

Satisfying the CFL condition prevents catastrophic instability, but it does not guarantee that the result is accurate. A stable method may still have phase error, amplitude error or artificial damping.

This distinction is central in this project. The physical equation already contains damping. Therefore, if a numerical method damps the wave too strongly, one must ask whether the damping is physical or artificial.

The explicit method is the clearest example. It may remain stable, but it can still introduce numerical dissipation. This means the amplitude may decay faster than in a more accurate reference method such as the matrix exponential.

The comparison between methods is therefore not just about whether the animations look reasonable. It is about understanding how each algorithm modifies the physical behaviour of the equation.
