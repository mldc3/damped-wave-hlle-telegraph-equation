# Theory

This project studies the damped one-dimensional wave / telegraph equation in the context of computational physics.

The equation combines wave propagation with dissipation. In this setting, the damping term reduces amplitude over time while preserving the underlying transport behavior.

The repository is limited to Practice 6 scope:
- damped 1D wave / telegraph dynamics,
- boundary conditions (Dirichlet and Neumann),
- comparison of numerical behavior across multiple time-integration and flux formulations.

## 13. Conservative formulation

Finite-difference methods approximate derivatives directly. For hyperbolic equations, however, it is often better to write the system in conservative form:

$$
\frac{\partial \mathbf{U}}{\partial t}
+
\frac{\partial \mathbf{F}(\mathbf{U})}{\partial x}
=
\mathbf{S}(\mathbf{U}).
$$

This formulation emphasizes fluxes through cell interfaces. Instead of asking only how a derivative behaves at a point, it asks how much quantity enters and leaves each computational cell. This is the natural setting of finite-volume methods.

For the wave problem, introduce another auxiliary variable

$$
w=\frac{\partial u}{\partial x}.
$$

The state vector can be written as

$$
\mathbf{U}
=
\begin{pmatrix}
u\\
v\\
w
\end{pmatrix},
$$

where $u$ is the field, $v=u_t$ is the velocity-like variable and $w=u_x$ is the spatial-gradient variable.

A possible conservative flux form for the wave part is

$$
\mathbf{F}(\mathbf{U})
=
\begin{pmatrix}
0\\
-c^2w\\
-v
\end{pmatrix}.
$$

The source terms then contain damping and reaction contributions, for example terms proportional to $-2(\kappa/\rho)v$ and $a u$.

## 14. Jacobian and characteristic speeds

The Jacobian of the flux is

$$
\mathbf{J}
=
\frac{\partial \mathbf{F}}{\partial \mathbf{U}}.
$$

For the wave system, the eigenvalues of this Jacobian are

$$
\lambda_1=0,
\qquad
\lambda_2=c,
\qquad
\lambda_3=-c.
$$

These eigenvalues are characteristic speeds. They tell us how information travels through the system. The two nonzero speeds correspond to waves propagating left and right. The zero speed corresponds to a non-propagating component associated with the chosen variables.

This characteristic structure is the key to upwind methods and approximate Riemann solvers. Instead of treating all directions symmetrically, upwind methods use the direction in which information actually travels.

## 15. Riemann problem at cell interfaces

In a finite-volume method, the domain is divided into cells. To update a cell, one needs the flux through its left and right interfaces. At each interface there is a left state and a right state:

$$
\mathbf{U}_L,
\qquad
\mathbf{U}_R.
$$

This defines a local Riemann problem: an initial-value problem with a discontinuity at the interface. The exact solution of the Riemann problem may contain several waves. In complicated systems, solving it exactly is expensive. Approximate solvers are therefore used.

The HLL/HLLE method replaces the full wave structure by an approximate two-wave model. It only keeps track of the fastest left-going and right-going waves.

The left and right signal speeds are estimated as

$$
S_L=\min(\lambda_{\min},0),
$$

$$
S_R=\max(\lambda_{\max},0).
$$

For the linear wave system, these are typically

$$
S_L=-c,
\qquad
S_R=c.
$$

## 16. HLL/HLLE numerical flux

The HLL flux is

$$
\mathbf{F}_{HLL}
=
\frac{
S_R\mathbf{F}_L
-
S_L\mathbf{F}_R
+
S_LS_R
\left(
\mathbf{U}_R-\mathbf{U}_L
\right)
}{
S_R-S_L
}.
$$

This formula has a clear interpretation.

If all waves move to the right, the interface flux should be the left flux. If all waves move to the left, the interface flux should be the right flux. If waves move in both directions, HLL uses an intermediate average that accounts for the jump between the two states.

The final term,

$$
S_LS_R(\mathbf{U}_R-\mathbf{U}_L),
$$

is responsible for much of the stabilizing numerical viscosity of the method. This term damps sharp jumps and suppresses spurious oscillations. That is why HLL/HLLE is robust.

## 17. Numerical viscosity in HLL/HLLE

Numerical viscosity is not necessarily a bug. In hyperbolic conservation laws, especially with discontinuities or shocks, numerical viscosity prevents nonphysical oscillations. However, in a smooth linear wave problem, the same viscosity can make the result look more damped than the physical model alone would predict.

This is why HLL/HLLE is expected to look different from matrix exponential or Crank–Nicolson. It is not simply solving the same problem in the same numerical philosophy. It is solving it with a conservative interface-flux method that deliberately adds robustness at the cost of extra smoothing.

In the animations, this may appear as:
- faster amplitude decay,
- smoother wave profiles,
- less sharp spatial structure,
- stronger damping of high-frequency components.

This behaviour is expected and should be interpreted correctly. HLL is particularly valuable because it prepares the code structure for more advanced hyperbolic systems, including nonlinear conservation laws.
