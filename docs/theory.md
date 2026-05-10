# Theory

This project studies the damped one-dimensional wave / telegraph equation in the context of computational physics.

The equation combines wave propagation with dissipation. In this setting, the damping term reduces amplitude over time while preserving the underlying transport behavior.

The repository is limited to Practice 6 scope:
- damped 1D wave / telegraph dynamics,
- boundary conditions (Dirichlet and Neumann),
- comparison of numerical behavior across multiple time-integration and flux formulations.

## 9. Matrix formulation

The first-order system can be written with a block matrix:

$$
\mathbf{A}
=
\begin{pmatrix}
\mathbf{0} & \mathbf{I}\\
c^2\mathbf{L}+a\mathbf{I} & -2(\kappa/\rho)\mathbf{I}
\end{pmatrix}.
$$

Here $\mathbf{I}$ is the identity matrix and $\mathbf{L}$ is the discrete Laplacian. The block structure has a direct physical interpretation.

The upper-right identity block represents

$$
\frac{d\mathbf{u}}{dt}=\mathbf{v}.
$$

The lower-left block represents the restoring dynamics due to spatial curvature and the linear source term:

$$
\frac{d\mathbf{v}}{dt}
=
(c^2\mathbf{L}+a\mathbf{I})\mathbf{u}
-
2(\kappa/\rho)\mathbf{v}.
$$

The lower-right block represents physical damping of the velocity-like variable.

The discrete Laplacian for a one-dimensional grid with centred second-order finite differences has the tridiagonal form

$$
\mathbf{L}
=
\frac{1}{\Delta x^2}
\begin{pmatrix}
-2 & 1 & 0 & \cdots & 0\\
1 & -2 & 1 & \cdots & 0\\
0 & 1 & -2 & \cdots & 0\\
\vdots & \vdots & \vdots & \ddots & 1\\
0 & 0 & 0 & 1 & -2
\end{pmatrix}.
$$

The boundary rows may be modified depending on whether Dirichlet or Neumann boundary conditions are imposed. This is one reason why boundary conditions must be handled carefully in matrix methods: they are encoded directly into the operator that generates the time evolution.

## 10. Matrix exponential method

For a linear system

$$
\frac{d\mathbf{W}}{dt}
=
\mathbf{A}\mathbf{W},
$$

with constant matrix $\mathbf{A}$, the exact solution over one time step is

$$
\mathbf{W}^{n+1}
=
e^{\mathbf{A}\Delta t}
\mathbf{W}^{n}.
$$

The matrix exponential method computes

$$
\mathbf{M}
=
e^{\mathbf{A}\Delta t}
$$

and then advances the solution by repeated multiplication:

$$
\mathbf{W}^{n+1}=\mathbf{M}\mathbf{W}^{n}.
$$

This method is conceptually very clean. Once the spatial grid and boundary conditions have been chosen, the semi-discrete linear system is solved exactly in time. Therefore, the matrix exponential provides a strong reference solution for the other time-integration schemes.

Its main disadvantage is computational cost. Computing a matrix exponential can be expensive for large systems. The method is excellent for moderate one-dimensional problems and for validation, but it may become impractical for very large multidimensional simulations.

## 11. Crank–Nicolson method

Crank–Nicolson is a semi-implicit second-order method. Applied to

$$
\frac{d\mathbf{W}}{dt}
=
\mathbf{A}\mathbf{W},
$$

it uses the average of the right-hand side at time levels $n$ and $n+1$:

$$
\frac{\mathbf{W}^{n+1}-\mathbf{W}^{n}}{\Delta t}
=
\frac{1}{2}
\mathbf{A}
\left(
\mathbf{W}^{n+1}
+
\mathbf{W}^{n}
\right).
$$

Rearranging gives

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

Thus each time step requires solving a linear system. The method can also be interpreted as a rational approximation to the matrix exponential:

$$
e^{\mathbf{A}\Delta t}
\approx
\left(
\mathbf{I}
-
\frac{\Delta t}{2}\mathbf{A}
\right)^{-1}
\left(
\mathbf{I}
+
\frac{\Delta t}{2}\mathbf{A}
\right).
$$

This explains why Crank–Nicolson often resembles the matrix exponential solution much more closely than an explicit first-order scheme. It is second order in time, stable for many linear problems, and usually has much less artificial damping than simple explicit or backward Euler methods.

## 12. Why Crank–Nicolson is a good compromise

The matrix exponential is theoretically elegant but expensive. The explicit method is cheap but conditionally stable and can be dissipative. Crank–Nicolson lies between them.

It is more expensive per time step than the explicit method because it solves a linear system. However, it permits larger time steps and gives much better long-time behaviour. It does not exactly reproduce the matrix exponential, but for this linear damped wave problem it should track the correct phase and damping very well when the time step is chosen reasonably.

In this repository, Crank–Nicolson is therefore an important method: it demonstrates how implicit and semi-implicit schemes can preserve the qualitative structure of wave propagation better than direct explicit updates while remaining more practical than a full matrix exponential approach.
