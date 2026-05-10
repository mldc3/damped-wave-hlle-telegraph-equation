## 5. Finite-difference grid and notation

The continuous domain is replaced by a discrete grid. Let

$$
x_j=x_{\min}+j\Delta x,
\qquad
t^n=n\Delta t.
$$

The numerical approximation to the continuous solution is denoted by

$$
u_j^n\approx u(x_j,t^n).
$$

This notation makes clear that the solution is evaluated at discrete spatial nodes and discrete time levels. The goal of a finite-difference method is to replace derivatives by algebraic combinations of neighbouring grid values.

For the second time derivative, a centred approximation gives

$$
\frac{\partial^2 u}{\partial t^2}
\approx
\frac{u_j^{n+1}-2u_j^n+u_j^{n-1}}{\Delta t^2}.
$$

For the second spatial derivative,

$$
\frac{\partial^2 u}{\partial x^2}
\approx
\frac{u_{j+1}^{n}-2u_j^n+u_{j-1}^{n}}{\Delta x^2}.
$$

For the damping term, a centred approximation in time is

$$
\frac{\partial u}{\partial t}
\approx
\frac{u_j^{n+1}-u_j^{n-1}}{2\Delta t}.
$$

Substituting these approximations into the damped wave equation gives a fully discrete recurrence relation. This is the most direct way to solve the problem, because it updates $u_j^{n+1}$ from already known values at previous time levels.

## 6. Explicit finite-difference scheme

After substituting the finite differences into the telegraph equation, the explicit update can be written in the form

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
},
$$

where

$$
k=\frac{\kappa}{\rho},
\qquad
r=\left(\frac{c\Delta t}{\Delta x}\right)^2.
$$

The exact sign convention of the $a\Delta t^2$ contribution depends on how the original equation is rearranged before discretization. What matters conceptually is that the explicit method computes the future field directly from the current and previous states.

The explicit scheme has two major advantages. First, it is easy to understand: each new value depends only on nearby spatial values and previous time levels. Second, it is computationally cheap per time step because no linear system has to be solved.

However, explicit hyperbolic schemes are conditionally stable. The time step cannot be chosen arbitrarily. For the centred wave equation, the Courant number must satisfy approximately

$$
\frac{c\Delta t}{\Delta x}\leq 1.
$$

Equivalently,

$$
r\leq 1.
$$

This is the Courant-Friedrichs-Lewy condition. Its physical meaning is that the numerical domain of dependence must contain the physical domain of dependence. If the time step is too large, the numerical method attempts to move information farther than the grid can represent in one step. The result is typically an unstable oscillation or complete divergence.

## 7. Stability is not the same as accuracy

A method satisfying the CFL condition may still be inaccurate. Stability only means that numerical errors do not grow without bound. It does not guarantee that phase, amplitude and damping are correct.

For wave problems, this distinction is essential. A numerical method can remain stable while introducing phase error, meaning the wave travels with a slightly wrong numerical speed. It can also introduce amplitude error, meaning that the wave decays too slowly or too quickly. In this project, the most important amplitude error is artificial numerical damping.

The physical equation already contains damping. Therefore, any extra damping from the numerical method must be interpreted carefully. If the explicit method produces a wave that disappears significantly faster than the matrix exponential or Crank–Nicolson methods, the difference is not necessarily a new physical effect. It may simply be numerical dissipation.

## 8. Why use first-order systems?

The original equation is second order in time. Many powerful numerical methods, however, are naturally written for first-order systems. To obtain such a system, introduce the auxiliary variable

$$
v=\frac{\partial u}{\partial t}.
$$

Then the second-order equation becomes

$$
\frac{\partial u}{\partial t}=v,
$$

$$
\frac{\partial v}{\partial t}
=
c^2\frac{\partial^2u}{\partial x^2}
-
2\frac{\kappa}{\rho}v
+
a u.
$$

The variable $u$ plays the role of a displacement-like field, while $v$ plays the role of a velocity-like field. This reformulation separates the kinematic relation $u_t=v$ from the dynamic equation for $v_t$.

After discretizing space, the unknowns can be collected in a single state vector:

$$
\mathbf{W}
=
\begin{pmatrix}
u_0\\
u_1\\
\vdots\\
u_{N-1}\\
v_0\\
v_1\\
\vdots\\
v_{N-1}
\end{pmatrix}.
$$

The semi-discrete problem can then be written as

$$
\frac{d\mathbf{W}}{dt}
=
\mathbf{A}\mathbf{W}.
$$

This form is extremely useful because it allows the use of matrix exponential time evolution and Crank–Nicolson time integration.
