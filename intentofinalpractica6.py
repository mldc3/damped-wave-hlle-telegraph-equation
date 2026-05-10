import numpy as np
import matplotlib.pyplot as plt
from scipy.sparse import diags, csr_matrix, lil_matrix
from scipy.sparse.linalg import spsolve
from scipy.linalg import expm
import time
import warnings
import matplotlib.animation as animation
import os
import scipy.linalg as sla

x_min = 0                               # Donde empieza el intervalo espacial
x_max = 2                               # Donde termina el intervalo espacial
t_max = 0                               # Tiempo inicial
t_max = 40                              # Tiempo final
kappa_rho = 0.1                         # Coef. amortiguamiento \frac{\kappa}{\rho}
a = -2                                  # Término fuente lineal
c = 1.0                                 # Velocidad de la onda

# Discretización
Nx = 201                               # Número de puntos espaciales
Nt = 4001                              # Número de puntos temporales

# Paso espacial y temporal
dx = (x_max - x_min) / (Nx - 1)
dt = t_max / (Nt - 1)

# Malla espacial
x = np.linspace(x_min, x_max, Nx)

state_size = 2 * Nx

#para que no tarde tanto el gif ponemos:
SKIP_FRAMES = 10




##############################PRIMER MÉTODO: explícito############################## 
def laplacian_interior(u, dx): #Calcula el Laplaciano discreto en los puntos interiores \frac{\partial^2 u}{\partial x^2} = \frac{u[i+1] - 2u[i] + u[i-1]}{dx^2}
    L = np.zeros_like(u)                                    # Inicializar el vector del Laplaciano
    L[1:-1] = (u[2:] - 2.0 * u[1:-1] + u[:-2]) / (dx**2)    # Aplicar fórmula de diferencia centrada en puntos interiores (1 a Nx-2): 
    return L

def explicit_method_dirichlet(u_initial, u_previous, c, kappa_rho, a, dx, dt, Nt):
    Nx = len(u_initial)                 # Número de puntos espaciales (malla en x)
    u = np.zeros((Nt, Nx))              # Matriz para almacenar la solución
    u[0] = u_initial                    # Condición inicial en t=0
    u[1] = u_previous                   # Condición inicial en t=dt (derivada del reposo)
    
    r = (c * dt / dx)**2                # Coeficiente de estabilidad (no se multiplica por 0.98)
    k = kappa_rho 
    
    for n in range(1, Nt-1):            # Iterar sobre el tiempo
        for j in range(1, Nx-1):        # Iterar sobre el espacio
            #u[n+1][j] = (2*u[n][j] - (1 - kappa_rho*dt)*u[n-1][j] + r*(u[n][j+1] - 2*u[n][j] + u[n][j-1]) - a*(dt**2)*u[n][j]) / (1 + kappa_rho*dt)
            numer = (2 - k*dt - a*(dt**2)) * u[n, j] \
            - u[n-1, j] \
            + r * (u[n, j+1] - 2*u[n, j] + u[n, j-1])
            u[n+1, j] = numer / (1 + k*dt)
        # Aplicar condiciones de contorno dirichet
        u[n+1][0] = 0.0                 # Condición de Dirichlet en el extremo izquierdo
        u[n+1][-1] = 0.0                # Condición de Dirichlet en el extremo derecho
    
    return u

def explicit_method_neumann(u_initial, u_previous, c, kappa_rho, a, dx, dt, Nt):
    Nx = len(u_initial)                 # Número de puntos espaciales (malla en x)
    u = np.zeros((Nt, Nx))              # Matriz para almacenar la solución en tiempo (Nt) y espacio (Nx)
    u[0] = u_initial                    # Condición inicial en t = 0 (desplazamiento inicial)
    u[1] = u_previous                   # Condición en t = dt (por ejemplo usando derivada temporal inicial / reposo)

    r = (c * dt / dx)**2                # r = (c*dt/dx)^2, coeficiente que multiplica al Laplaciano discreto (estabilidad/propagación)
    k = kappa_rho                       # Alias local para claridad: amortiguamiento combinado kappa/rho

    for n in range(1, Nt-1):            # Bucle temporal: para cada paso n (usamos u[n-1], u[n] para calcular u[n+1])

        for j in range(1, Nx-1):        # Bucle espacial: actualizamos sólo en puntos interiores j = 1 .. Nx-2
            numer = (2 - k*dt - a*(dt**2)) * u[n, j] \
            - u[n-1, j] \
            + r * (u[n, j+1] - 2*u[n, j] + u[n, j-1])
            u[n+1, j] = numer / (1 + k*dt)
            # División final por (1 + k*dt) que aparece por reorganizar el esquema (término de amortiguamiento en denominador)
            # Así se obtiene el valor de u en el tiempo n+1 para el punto espacial j.
        u[n+1][-1] = u[n+1][-2]         # Condición de contorno Neumann en el extremo derecho: \frac{\partial u}{\partial x} = 0 entinces u[-1] = u[-2]
    return u                            


##############################SEGUNDO MÉTODO: sistema de primer orden explícito############################## 

def build_system_matrix_dirichlet(c, kappa_rho, a, dx, Nx):

    O = np.zeros((Nx, Nx))
    I = np.eye(Nx)

    # Laplaciano espacial L con diagonales
    diag_main  = -2.0 * np.ones(Nx)
    diag_upper = 1.0  * np.ones(Nx - 1)
    diag_lower = 1.0  * np.ones(Nx - 1)
    L = (np.diag(diag_main) + np.diag(diag_upper, 1) + np.diag(diag_lower, -1)) / (dx**2)

    # Bloques inferiores
    abajoizq  = c**2 * L + a * I         
    abajoder = -2.0 * kappa_rho * I    

    A = np.block([[O, I],
                  [abajoizq, abajoder]])

    # Imponemos Dirichlet en u_0 y u_{N-1} y también fijamos las correspondientes ecuaciones para v (filas de A que hacen que u_0, u_{N-1}, v_0, v_{N-1} sean constantes)
    A[0, :] = 0.0
    A[0, 0] = 1.0

    A[Nx-1, :] = 0.0
    A[Nx-1, Nx-1] = 1.0

    A[Nx, :] = 0.0
    A[Nx, Nx] = 1.0

    A[2*Nx-1, :] = 0.0
    A[2*Nx-1, 2*Nx-1] = 1.0

    return A


def build_system_matrix_neumann(c, kappa_rho, a, dx, Nx):
    O = np.zeros((Nx, Nx))
    I = np.eye(Nx)

    # Construcción del Laplaciano L
    L = np.zeros((Nx, Nx))
    for i in range(1, Nx-1):
        L[i, i-1] = 1.0
        L[i, i]   = -2.0
        L[i, i+1] = 1.0
    L = L / (dx**2)

    # fila i=0: Dirichlet en x=0: anular fila (la condición la aplicaremos también en A)
    L[0, :] = 0.0

    # fila i = Nx-1: Neumann en extremo derecho con punto fantasma u_N = u_{N-2}
    L[Nx-1, :] = 0.0
    L[Nx-1, Nx-2] =  2.0 / (dx**2)
    L[Nx-1, Nx-1] = -2.0 / (dx**2)

    lower_left  = c**2 * L + a * I
    lower_right = -2.0 * kappa_rho * I

    A = np.block([[O, I],
                  [lower_left, lower_right]])

    # Fijamos Dirichlet en x=0: forzamos u_0 y v_0 
    A[0, :]   = 0.0
    A[0, 0]   = 1.0
    A[Nx, :]  = 0.0
    A[Nx, Nx] = 1.0
    
    return A



# Lo resolvemos mediante el método de exponencial matricial (método exacto)

def solve_matrix_exponential(u_initial, c, kappa_rho, a, dx, dt, Nt, boundary="dirichlet"):

    Nx = len(u_initial)

    if boundary.lower() == "dirichlet":
        A = build_system_matrix_dirichlet(c, kappa_rho, a, dx, Nx)
    elif boundary.lower() == "neumann":
        A = build_system_matrix_neumann(c, kappa_rho, a, dx, Nx)
    else:
        raise ValueError("boundary debe ser 'dirichlet' o 'neumann'")

    # Matriz de transición
    M = expm(A * dt)

    # estado inicial W = [u, v] (v = 0 por reposo)
    W = np.zeros(2 * Nx)
    W[:Nx] = u_initial

    u_solution = np.zeros((Nt, Nx))
    u_solution[0] = u_initial

    for n in range(1, Nt):
        W = M @ W
        u_solution[n] = W[:Nx]

    return u_solution


# Crank-Nicolson aplicado al sistema en primer orden

def solve_crank_nicolson_dirichlet(u_initial, c, kappa_rho, a, dx, dt, Nt):
    Nx = len(u_initial)
    A = build_system_matrix_dirichlet(c, kappa_rho, a, dx, Nx)

    I2N = np.eye(2 * Nx)
    B = I2N - 0.5 * dt * A   # (I - \frac{dt}{2} A)
    C = I2N + 0.5 * dt * A   # (I + \frac{dt}{2} A)

    # A estas alturas A ya tiene filas de frontera fijadas; sin embargo por seguridad
    # fijamos también esas mismas filas en B y C para que W^{n+1} mantenga los valores
    
    # Fila 0
    B[0, :] = 0.0
    B[0, 0] = 1.0
    C[0, :] = 0.0
    C[0, 0] = 1.0

    # Fila Nx-1
    B[Nx-1, :] = 0.0
    B[Nx-1, Nx-1] = 1.0
    C[Nx-1, :] = 0.0
    C[Nx-1, Nx-1] = 1.0

    # Fila Nx
    B[Nx, :] = 0.0
    B[Nx, Nx] = 1.0
    C[Nx, :] = 0.0
    C[Nx, Nx] = 1.0

    # Fila 2*Nx-1
    B[2*Nx-1, :] = 0.0
    B[2*Nx-1, 2*Nx-1] = 1.0
    C[2*Nx-1, :] = 0.0
    C[2*Nx-1, 2*Nx-1] = 1.0


    # Inicialización
    W = np.zeros(2 * Nx)
    W[:Nx] = u_initial
    u_solution = np.zeros((Nt, Nx))
    u_solution[0] = u_initial

    # Resolver paso a paso: B W^{n+1} = C W^n
    lu, piv = sla.lu_factor(B)              # Para eficiencia si Nt grande, factorizar B (por ej. LU)
    for n in range(1, Nt):
        rhs = C @ W
        W = sla.lu_solve((lu, piv), rhs)
        u_solution[n] = W[:Nx]

    return u_solution


def solve_crank_nicolson_neumann(u_initial, c, kappa_rho, a, dx, dt, Nt):

    Nx = len(u_initial)
    A = build_system_matrix_neumann(c, kappa_rho, a, dx, Nx)

    I2N = np.eye(2 * Nx)
    B = I2N - 0.5 * dt * A
    C = I2N + 0.5 * dt * A

    # Fijar filas correspondientes a x=0 (u_0 y v_0) para imponer Dirichlet en ese extremo. Solo x=0 (izquierda) en tu enunciado; no anular filas del extremo derecho
    # Fila 0
    B[0, :] = 0.0
    B[0, 0] = 1.0
    C[0, :] = 0.0
    C[0, 0] = 1.0

    # Fila Nx
    B[Nx, :] = 0.0
    B[Nx, Nx] = 1.0
    C[Nx, :] = 0.0
    C[Nx, Nx] = 1.0
    
    # Inicialización y resolución (factorizamos B)
    W = np.zeros(2 * Nx)
    W[:Nx] = u_initial
    u_solution = np.zeros((Nt, Nx))
    u_solution[0] = u_initial

    lu, piv = sla.lu_factor(B)
    for n in range(1, Nt):
        rhs = C @ W
        W = sla.lu_solve((lu, piv), rhs)
        u_solution[n] = W[:Nx]

    return u_solution


##############################CUARTO MÉTODO: HLL N############################## 
def compute_eigenvalues(W_state, c):
    # Dado un estado W_state, devolvemos los autovalores del sistema linealizado. En este problema son constantes: lambda_plus = c y lambda_minus = -c, no lo hacemos con la funcion de python porque sino tarda mas y los valores son constantes
    return c, -c # que representan las velocidades características hacia la derecha e izquierda.


def compute_hll_flux(W_L, W_R, c):
    # Autovalores (velocidades características) en el estado izquierdo
    lambda_plus_L, lambda_minus_L = compute_eigenvalues(W_L, c)
    # Autovalores (velocidades características) en el estado derecho
    lambda_plus_R, lambda_minus_R = compute_eigenvalues(W_R, c)
    

    s_plus = max(0.0, lambda_plus_L, lambda_plus_R)     # Velocidad de onda más rápida hacia la derecha: s_plus = max(0, lambda_plus_L, lambda_plus_R)
    s_minus = min(0.0, lambda_minus_L, lambda_minus_R)  # Velocidad de onda más rápida hacia la izquierda: s_minus = min(0, lambda_minus_L, lambda_minus_R)

    # Desempaquetamos las componentes del estado: u: desplazamiento, v: velocidad temporal, w: derivada espacial de u
    u_L, v_L, w_L = W_L
    u_R, v_R, w_R = W_R

    F_L = np.array([0.0, -c**2 * w_L, -v_L])    # Flujo físico en el lado izquierdo: F(u, v, w) = (0, -c^2 * w, -v)
    F_R = np.array([0.0, -c**2 * w_R, -v_R])    # Flujo físico en el lado derecho
    
    if abs(s_plus - s_minus) < 1e-14:           # Si s_plus = s_minus, el abanico de ondas colapsa (evitamos división por cero).
        F_HLL = 0.5 * (F_L + F_R)
    else:                                       # Fórmula del flujo HLL: F_HLL = ( s_plus * F_L - s_minus * F_R + s_plus * s_minus * (W_R - W_L) ) / ( s_plus - s_minus )
        F_HLL = (s_plus * F_L - s_minus * F_R + s_plus * s_minus * (W_R - W_L)) / (s_plus - s_minus)

    return F_HLL                                # Devolvemos el flujo numérico en la interfaz

def solve_hll_dirichlet(u_initial, c, kappa_rho, a, dx, dt, Nt):
    Nx = len(u_initial)                         # Número de puntos espaciales a partir del tamaño de u_initial
    W = np.zeros((Nx, 3))                       # W tiene dimensiones (Nx, 3) y almacena (u, v, w) en cada celda
    W[:,0] = u_initial                          # Inicializamos u en t = 0 con la condición inicial; v y w empiezan en cero
    u_solution = np.zeros((Nt, Nx))             # Matriz para guardar la evolución de u en el tiempo: u_solution[n, j] = u(x_j, t_n)
    u_solution[0] = u_initial

    for n in range(Nt-1):                       #Bucle principal en tiempo: avanzamos desde n = 0 hasta n = Nt - 2
        W[1:-1,2] = (W[2:,0] - W[:-2,0])/(2*dx) # Cálculo de w = \frac{\partial u}{\partial x} en los puntos interiores mediante diferencias centradas:w_i = ( u_{i+1} - u_{i-1} ) / (2 * dx)
        W[0,2] = (W[1,0] - W[0,0])/dx           # Extremo izquierdo: derivada hacia adelante: w_0 = ( u_1 - u_0 ) / dx
        W[-1,2] = (W[-1,0] - W[-2,0])/dx        # Extremo derecho: derivada hacia atrás w_{Nx-1} = ( u_{Nx-1} - u_{Nx-2} ) / dx

        F_interface = np.zeros((Nx+1,3)) 
        for i in range(1,Nx):                   # Recorremos las interfaces internas i = 1, ..., Nx-1
            F_interface[i] = compute_hll_flux(W[i-1], W[i], c)

        W_new = W.copy()                        # Creamos una copia del estado actual para construir el siguiente paso temporal
        for i in range(1,Nx-1):
            S = np.array([0.0, -2.0*kappa_rho*W[i,1] + a*W[i,0], 0.0])
            W_new[i] = W[i] + dt * (-(F_interface[i+1] - F_interface[i])/dx + S)

        # Condiciones de contorno Dirichlet
        W_new[0,0] = 0.0
        W_new[0,1] = 0.0
        W_new[-1,0] = 0.0
        W_new[-1,1] = 0.0

        W = W_new
        u_solution[n+1] = W[:,0]

    return u_solution

def solve_hll_neumann(u_initial, c, kappa_rho, a, dx, dt, Nt):
    Nx = len(u_initial)
    W = np.zeros((Nx,3))
    W[:,0] = u_initial
    u_solution = np.zeros((Nt,Nx))
    u_solution[0] = u_initial

    for n in range(Nt-1):
        # Calcular w = \frac{\partial u}{\partial x} 
        W[1:-1,2] = (W[2:,0] - W[:-2,0])/(2*dx)
        # Extremo izquierdo Dirichlet
        W[0,2] = (W[1,0] - W[0,0])/dx
        # Extremo derecho Neumann (libre): \frac{\partial u}{\partial x} = 0
        W[-1,2] = 0.0

        # Flujos HLL en las interfaces
        F_interface = np.zeros((Nx+1,3))
        for i in range(1,Nx):
            F_interface[i] = compute_hll_flux(W[i-1], W[i], c)

        W_new = W.copy()
        for i in range(1,Nx-1):
            S = np.array([0.0, -2.0*kappa_rho*W[i,1] + a*W[i,0], 0.0])
            W_new[i] = W[i] + dt * (-(F_interface[i+1] - F_interface[i])/dx + S)

        # Condiciones de contorno
        W_new[0,0] = 0.0            # izquierda fija (Dirichlet)
        W_new[0,1] = 0.0
        W_new[-1,0] = W_new[-2,0]   # derecha libre (Neumann), no lo he sustituido en la ecuación y he puesto esto ya que es lo que me dijiste de hacer en el primer método también
        W_new[-1,1] = W_new[-2,1]   # velocidad igual a la penúltima celda

        W = W_new
        u_solution[n+1] = W[:,0]

    return u_solution


###############################Condiciones iniciales y GIFs###############################


def initial_u(n):
    return np.sin(n * np.pi * x)


def u_at_dt(u0):
    # Calcular Laplaciano de u0
    L_u0 = laplacian_interior(u0, dx)
    
    # Aplicar fórmula
    return u0 + 0.5 * dt**2 * (c**2 * L_u0 + a * u0)



# Parámetro n para la condición inicial sin(n \pi x)
n_val = 1
    
# Calcular condición inicial
u0 = initial_u(n_val)
u1 = u_at_dt(u0)

# MÉTODO 1: EXPLÍCITO 

print("MÉTODO 1: ESQUEMA EXPLÍCITO")
print("")

print("DIRICHLET (extremos fijos):")
start_time = time.time()                          
U_explicit_dirichlet = explicit_method_dirichlet(u0, u1, c, kappa_rho, a, dx, dt, Nt)
dirichlet_explicit_time = time.time() - start_time 
print(f"Completo en {dirichlet_explicit_time:.3f}s\n")

print("NEUMANN (extremo izquierdo fijo, derecho libre):")
start_time = time.time()
U_explicit_neumann = explicit_method_neumann(u0, u1, c, kappa_rho, a, dx, dt, Nt)
neumann_explicit_time = time.time() - start_time
print(f"Completo en {neumann_explicit_time:.3f}s\n")

# MÉTODO 2: VARIABLES AUXILIARES 
print("MÉTODO 2: VARIABLES AUXILIARES")
start_time = time.time()
U_exp_dirichlet = solve_matrix_exponential(u0, c, kappa_rho, a, dx, dt, Nt, boundary="dirichlet")
dirichlet_exp_time = time.time() - start_time
print(f"Dirichlet completo en {dirichlet_exp_time:.3f}s\n")

start_time = time.time()
U_exp_neumann = solve_matrix_exponential(u0, c, kappa_rho, a, dx, dt, Nt, boundary="neumann")
neumann_exp_time = time.time() - start_time
print(f"Neumann completo en {neumann_exp_time:.3f}s\n")

print("MÉTODO 2b: VARIABLES AUXILIARES - CRANK-NICOLSON")
start_time = time.time()
U_cn_dirichlet = solve_crank_nicolson_dirichlet(u0, c, kappa_rho, a, dx, dt, Nt)
dirichlet_cn_time = time.time() - start_time
print(f"Dirichlet completo en {dirichlet_cn_time:.3f}s\n")

start_time = time.time()
U_cn_neumann = solve_crank_nicolson_neumann(u0, c, kappa_rho, a, dx, dt, Nt)
neumann_cn_time = time.time() - start_time
print(f"Neumann completo en {neumann_cn_time:.3f}s\n")


# MÉTODO 3: ESQUEMA HLL CONSERVATIVO
print("MÉTODO 3: ESQUEMA HLL CONSERVATIVO")
start_time = time.time()
U_hll_dirichlet = solve_hll_dirichlet(u0, c, kappa_rho, a, dx, dt, Nt)
dirichlet_hll_time = time.time() - start_time
print(f"Dirichlet completo en {dirichlet_hll_time:.3f}s\n")

start_time = time.time()
U_hll_neumann = solve_hll_neumann(u0, c, kappa_rho, a, dx, dt, Nt)
neumann_hll_time = time.time() - start_time
print(f"Neumann completo en {neumann_hll_time:.3f}s\n")


num_frames = (Nt - 1) // SKIP_FRAMES + 1  

# Figura 1: Métodos Explícitos 
print("GIF 1: Métodos Explícitos (Dirichlet vs Neumann)...")

fig1, ax1 = plt.subplots(figsize=(12, 5))
line_expl_dir, = ax1.plot(x, U_explicit_dirichlet[0], label='Explícito Dirichlet', linewidth=2.5)
line_expl_neu, = ax1.plot(x, U_explicit_neumann[0], label='Explícito Neumann', linewidth=2.5,linestyle='--')
ax1.set_ylim(-1.3, 1.3)
ax1.set_xlabel('Posición x', fontsize=12)
ax1.set_ylabel('u(x,t)', fontsize=12)
ax1.set_title(f'Métodos Explícitos', fontsize=13)
ax1.legend(fontsize=11)
ax1.grid(True, alpha=0.3)

def update_explicit(frame):
    n = SKIP_FRAMES * frame
    line_expl_dir.set_ydata(U_explicit_dirichlet[n])
    line_expl_neu.set_ydata(U_explicit_neumann[n])
    ax1.set_title(f'Métodos Explícitos')
    return line_expl_dir, line_expl_neu

ani1 = animation.FuncAnimation(fig1, update_explicit, frames=num_frames, interval=50, blit=True)
ani1.save("01_metodo_explicito.gif", writer='pillow', fps=12, dpi=80)
print("Guardado: 01_metodo_explicito.gif\n")
plt.close(fig1)

# Figura 2: Método Matriz Exponencial 
print("GIF 2: Matriz Exponencial...")

fig2, ax2 = plt.subplots(figsize=(12, 5))
line_matrix_dir, = ax2.plot(x, U_exp_dirichlet[0], label='Matriz Exponencial Dirichlet', linewidth=2.5)
line_matrix_neu, = ax2.plot(x, U_exp_neumann[0], label='Matriz Exponencial Neumann', linewidth=2.5,linestyle='--')
ax2.set_ylim(-1.3, 1.3)
ax2.set_xlabel('Posición x', fontsize=12)
ax2.set_ylabel('u(x,t)', fontsize=12)
ax2.set_title(f'Matriz Exponencial ', fontsize=13)
ax2.legend(fontsize=11)
ax2.grid(True, alpha=0.3)

def update_matrix(frame):
    n = SKIP_FRAMES * frame
    line_matrix_dir.set_ydata(U_exp_dirichlet[n])
    line_matrix_neu.set_ydata(U_exp_neumann[n])
    ax2.set_title(f'Matriz Exponencial')
    return line_matrix_dir, line_matrix_neu

ani2 = animation.FuncAnimation(fig2, update_matrix, frames=num_frames, interval=50, blit=True)
ani2.save("02_metodo_matriz_exponencial.gif", writer='pillow', fps=12, dpi=80)
print("Guardado: 02_metodo_matriz_exponencial.gif\n")
plt.close(fig2)

# Figura 3: Crank-Nicolson
print("GIF 3: Crank-Nicolson (Dirichlet vs Neumann)...")

fig3, ax3 = plt.subplots(figsize=(12, 5))
line_cn_dir, = ax3.plot(x, U_cn_dirichlet[0], label='Crank-Nicolson Dirichlet', linewidth=2.5)
line_cn_neu, = ax3.plot(x, U_cn_neumann[0], label='Crank-Nicolson Neumann', linewidth=2.5, linestyle='--')
ax3.set_ylim(-1.3, 1.3)
ax3.set_xlabel('Posición x', fontsize=12)
ax3.set_ylabel('u(x,t)', fontsize=12)
ax3.set_title(f'Crank-Nicolson', fontsize=13)
ax3.legend(fontsize=11)
ax3.grid(True, alpha=0.3)

def update_cn(frame):
    n = SKIP_FRAMES * frame
    line_cn_dir.set_ydata(U_cn_dirichlet[n])
    line_cn_neu.set_ydata(U_cn_neumann[n])
    ax3.set_title(f'Crank-Nicolsonn')
    return line_cn_dir, line_cn_neu

ani3 = animation.FuncAnimation(fig3, update_cn, frames=num_frames, interval=50, blit=True)
ani3.save("03_metodo_crank_nicolson.gif", writer='pillow', fps=12, dpi=80)
print("Guardado: 03_metodo_crank_nicolson.gif\n")
plt.close(fig3)


# Figura 3: Método HLL 
print("GIF 4: Método HLL Conservativo...")

fig3, ax3 = plt.subplots(figsize=(12, 5))
line_hll_dir, = ax3.plot(x, U_hll_dirichlet[0], label='HLL Dirichlet', linewidth=2.5)
line_hll_neu, = ax3.plot(x, U_hll_neumann[0], label='HLL Neumann', linewidth=2.5, linestyle='--')
ax3.set_ylim(-1.3, 1.3)
ax3.set_xlabel('Posición x', fontsize=12)
ax3.set_ylabel('u(x,t)', fontsize=12)
ax3.set_title(f'HLL', fontsize=13)
ax3.legend(fontsize=11)
ax3.grid(True, alpha=0.3)

def update_hll(frame):
    n = SKIP_FRAMES * frame
    line_hll_dir.set_ydata(U_hll_dirichlet[n])
    line_hll_neu.set_ydata(U_hll_neumann[n])
    ax3.set_title(f'HLL')
    return line_hll_dir, line_hll_neu

ani3 = animation.FuncAnimation(fig3, update_hll, frames=num_frames, interval=50, blit=True)
ani3.save("04_metodo_hll.gif", writer='pillow', fps=12, dpi=80)
print("Guardado: 04_metodo_hll.gif\n")
plt.close(fig3)

# Figura 4: Comparación de los tres métodos (Dirichlet) 
print("GIF 5: Comparación Métodos (Dirichlet)...")

fig4, ax4 = plt.subplots(figsize=(12, 5))
line_comp_expl, = ax4.plot(x, U_explicit_dirichlet[0], label='Explícito', linewidth=2.5)
line_comp_matrix, = ax4.plot(x, U_exp_dirichlet[0], label='Matriz Exponencial', linewidth=3.5)
line_comp_hll, = ax4.plot(x, U_hll_dirichlet[0], label='HLL', linewidth=2.5)
line_comp_cn, = ax4.plot(x, U_cn_dirichlet[0], label='Crank-Nicolson', linewidth=1.5)
ax4.set_ylim(-1.3, 1.3)
ax4.set_xlabel('Posición x', fontsize=12)
ax4.set_ylabel('u(x,t)', fontsize=12)
ax4.set_title(f'Comparación (Dirichlet)', fontsize=13)
ax4.legend(fontsize=11)
ax4.grid(True, alpha=0.3)

def update_comparison_dir(frame):
    n = SKIP_FRAMES * frame
    line_comp_expl.set_ydata(U_explicit_dirichlet[n])
    line_comp_matrix.set_ydata(U_exp_dirichlet[n])
    line_comp_hll.set_ydata(U_hll_dirichlet[n])
    line_comp_cn.set_ydata(U_cn_dirichlet[n])
    ax4.set_title(f'Comparación Dirichlet')
    return line_comp_expl, line_comp_matrix, line_comp_hll, line_comp_cn

ani4 = animation.FuncAnimation(fig4, update_comparison_dir, frames=num_frames, interval=50, blit=True)
ani4.save("05_comparacion_dirichlet.gif", writer='pillow', fps=12, dpi=80)
print("Guardado: 05_comparacion_dirichlet.gif\n")
plt.close(fig4)

# Figura 5: Comparación Métodos (Neumann) 
print("GIF 6: Comparación Métodos (Neumann)...")

fig5, ax5 = plt.subplots(figsize=(12, 5))
line_comp_expl_n, = ax5.plot(x, U_explicit_neumann[0], label='Explícito', linewidth=2.5)
line_comp_matrix_n, = ax5.plot(x, U_exp_neumann[0], label='Matriz Exponencial', linewidth=3.5)
line_comp_hll_n, = ax5.plot(x, U_hll_neumann[0], label='HLL', linewidth=2.5)
line_comp_cn_n, = ax5.plot(x, U_cn_neumann[0], label='Crank-Nicolson', linewidth=1.5)
ax5.set_ylim(-1.3, 1.3)
ax5.set_xlabel('Posición x', fontsize=12)
ax5.set_ylabel('u(x,t)', fontsize=12)
ax5.set_title(f'Comparación (Neumann)', fontsize=13)
ax5.legend(fontsize=11)
ax5.grid(True, alpha=0.3)

def update_comparison_neu(frame):
    n = SKIP_FRAMES * frame
    line_comp_expl_n.set_ydata(U_explicit_neumann[n])
    line_comp_matrix_n.set_ydata(U_exp_neumann[n])
    line_comp_hll_n.set_ydata(U_hll_neumann[n])
    line_comp_cn_n.set_ydata(U_cn_neumann[n])
    ax5.set_title(f'Comparación Neumann')
    return line_comp_expl_n, line_comp_matrix_n, line_comp_hll_n, line_comp_cn_n

ani5 = animation.FuncAnimation(fig5, update_comparison_neu, frames=num_frames, interval=50, blit=True)
ani5.save("06_comparacion_neumann.gif", writer='pillow', fps=12, dpi=80)
print("Guardado: 06_comparacion_neumann.gif\n")
plt.close(fig5)


    
    














