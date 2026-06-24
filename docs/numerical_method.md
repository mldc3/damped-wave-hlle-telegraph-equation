# Numerical Method and Implementation Notes

This document explains how the numerical methods in the repository are organized at the implementation level. The focus is not on re-deriving all the theory, but on connecting the mathematical formulation with the structure of the Python code.

## 1. Computational workflow

The project follows the same general workflow for all methods:

1. define the spatial grid,
2. define the time step and final time,
3. initialize the field $u(x,0)$,
4. initialize the velocity $u_t(x,0)$,
5. impose boundary conditions,
6. evolve the solution with a chosen method,
7. store frames for visualization,
8. compare methods through GIF animations.

The methods differ in how the time update is performed, but they represent the same underlying physical equation.

## 2. Explicit finite-difference method

The explicit method stores the solution at two previous time levels and computes the next one directly. The field is represented as an array of grid values,

$$
u_j^n \approx u(x_j,t^n).
$$

The method uses finite differences for the second time derivative, second spatial derivative and damping term. At each step, the code updates the interior points and then applies the selected boundary conditions.

This method is simple and fast because it does not require matrix factorization or solving a linear system. Its main restriction is the CFL condition,

$$
\frac{c\Delta t}{\Delta x}\leq 1.
$$

## 3. Matrix exponential method

The matrix exponential method rewrites the PDE as a first-order system using

$$
v=\frac{\partial u}{\partial t}.
$$

The numerical state is

$$
\mathbf{W}=(u_0,\ldots,u_{N-1},v_0,\ldots,v_{N-1})^T.
$$

The code constructs a matrix $\mathbf{A}$ such that

$$
\frac{d\mathbf{W}}{dt}
=
\mathbf{A}\mathbf{W}.
$$

The update is then

$$
\mathbf{W}^{n+1}
=
e^{\mathbf{A}\Delta t}
\mathbf{W}^{n}.
$$

This method is used as a high-quality reference for the time evolution of the semi-discrete system.

## 4. Crank--Nicolson method

Crank--Nicolson uses the same first-order matrix system, but advances the solution with

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

The implementation builds the left-hand and right-hand matrices once and then solves the linear system at each time step.

This method is more expensive than the explicit update but gives better stability and usually much better agreement with the matrix exponential solution.

## 5. HLL/HLLE method

The HLL/HLLE implementation uses a conservative representation of the wave problem. The state vector includes variables such as the field, its time derivative and its spatial derivative.

At each cell interface, the method constructs left and right states and computes an approximate HLL flux,

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

The solution is then updated by flux differences across cell interfaces. This is structurally different from the explicit finite-difference and matrix-based methods.

## 6. Boundary conditions

The code compares Dirichlet and Neumann boundary conditions.

For Dirichlet boundaries, the field is fixed at the endpoints:

$$
u(0,t)=0,
\qquad
u(L,t)=0.
$$

For Neumann boundaries, the derivative is constrained:

$$
\frac{\partial u}{\partial x}=0.
$$

The boundary conditions must be applied consistently in each method. In explicit methods, this usually means modifying endpoint values after each update. In matrix methods, it means constructing the operator with the correct boundary rows. In flux methods, it affects interface treatment near the edges.

## 7. GIF generation

The code produces GIF animations for:
- explicit finite-difference evolution,
- matrix exponential evolution,
- Crank--Nicolson evolution,
- HLL/HLLE evolution,
- Dirichlet method comparison,
- Neumann method comparison.

These GIFs are not only visual outputs. They are the main diagnostic figures used to compare the qualitative behaviour of the four numerical approaches.
