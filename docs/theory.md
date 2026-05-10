# Theory Background: Damped 1D Wave / Telegraph Equation

## 1. Physical motivation

This project studies the numerical solution of a damped one-dimensional wave equation, also known in this context as a telegraph-type equation. The objective is not only to obtain an animation of a wave, but to understand how different numerical formulations represent propagation, damping, boundary reflections and numerical dissipation.

The equation considered is

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

Here $u(x,t)$ is the field being evolved, $c$ is the wave propagation speed, $\kappa/\rho$ controls physical damping, and $a$ is a linear reaction or source coefficient. In the implementation used for this project, representative parameters are

$$
c=1.0,\qquad \frac{\kappa}{\rho}=0.1,\qquad a=-2.0,
$$

on a one-dimensional spatial domain approximately $x\in[0,2]$ and over a finite simulation time.

The left-hand side is the standard wave operator. If the right-hand side were zero, the equation would describe ideal propagation at finite speed. The additional term proportional to $\partial u/\partial t$ introduces physical damping: it removes energy from the wave and makes the amplitude decay in time. The term proportional to $u$ modifies the local restoring behaviour of the system. Since in this project $a<0$, the source term behaves like an additional restoring contribution rather than an unstable growth term.

This equation is useful pedagogically because it combines several key ideas from computational physics: hyperbolic propagation, damping, finite-difference discretization, matrix time evolution, implicit schemes, conservative formulations and approximate Riemann solvers.

## 2. Relation with the ideal wave equation

The ideal one-dimensional wave equation is

$$
\frac{\partial^2 u}{\partial t^2}
-
c^2
\frac{\partial^2 u}{\partial x^2}
=
0.
$$

This is a hyperbolic partial differential equation. Hyperbolic equations describe finite-speed propagation of information. A disturbance introduced at one point does not instantly affect the entire domain; instead, it travels along characteristic directions. In one spatial dimension, the characteristic speeds of the ideal wave equation are related to $+c$ and $-c$, corresponding to waves moving to the right and to the left.

The damped telegraph equation keeps the wave-like character but adds a mechanism by which the oscillation loses amplitude. This is why the project is especially useful for comparing numerical methods: all methods should reproduce the same physical decay, but they may add different amounts of artificial numerical damping.

It is important to distinguish these two effects.

Physical damping is part of the differential equation:

$$
-2\frac{\kappa}{\rho}\frac{\partial u}{\partial t}.
$$

Numerical damping is introduced by the discretization method. It depends on the grid spacing, time step, boundary treatment and chosen numerical flux. A method can be stable but still too dissipative, meaning that it may make the wave disappear faster than the physical model requires.

## 3. Initial conditions

Because the equation contains a second derivative in time, a single initial condition is not sufficient. We must prescribe both the initial field and the initial velocity:

$$
u(x,0)=u_0(x),
$$

$$
\frac{\partial u}{\partial t}(x,0)=v_0(x).
$$

This is analogous to classical mechanics. A second-order equation for a particle requires initial position and initial velocity. A second-order wave equation requires the initial shape of the field and the initial rate of change of that field.

A typical initial condition used in the project is a sinusoidal mode,

$$
u(x,0)=\sin\left(\frac{n\pi x}{L}\right),
$$

with zero initial velocity,

$$
\frac{\partial u}{\partial t}(x,0)=0.
$$

The integer $n$ labels the spatial mode. Low values of $n$ correspond to long wavelengths and smooth solutions. Higher values contain more oscillations over the same domain and are therefore more sensitive to spatial resolution and numerical viscosity.

## 4. Boundary conditions

Boundary conditions are part of the physical problem. They are not a minor programming detail. They determine which modes are allowed, how waves reflect at the ends of the domain, and whether the numerical solution corresponds to the intended physical system.

For Dirichlet boundary conditions, the field is fixed at the boundaries:

$$
u(0,t)=0,\qquad u(L,t)=0.
$$

Physically, this resembles a string fixed at both ends. The endpoints cannot move, so the wave must vanish there.

For Neumann boundary conditions, the spatial derivative is prescribed. A homogeneous Neumann condition is

$$
\frac{\partial u}{\partial x}(L,t)=0.
$$

This corresponds to a zero-gradient or free-end condition. In a finite-difference method, this is usually implemented using a one-sided difference or a ghost-point relation. Dirichlet and Neumann conditions produce visibly different reflection behaviour, so comparing them is a meaningful part of the numerical study.
