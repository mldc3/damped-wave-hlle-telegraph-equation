import numpy as np
import matplotlib.pyplot as plt
# ------------------------------- generación de la red (práctica anterior)
sigma=2.3151

def generar_red_fcc(nx, ny, nz, a):
    # Definimos los vectores primitivos de la celda fcc (en unidades de a)
    base = [
        (0.0, 0.0, 0.0),
        (0.5, 0.5, 0.0),
        (0.5, 0.0, 0.5),
        (0.0, 0.5, 0.5)
    ]
    posiciones = []
    for ix in range(nx):
        for iy in range(ny):
            for iz in range(nz):
                # Calculamos el origen de la celda actual
                ox = ix * a
                oy = iy * a
                oz = iz * a
                # Para cada átomo de la base, calculamos su posición absoluta
                for bx, by, bz in base:
                    x = ox + bx * a    # suma la posición relativa de la base en x
                    y = oy + by * a    # suma la posición relativa de la base en y
                    z = oz + bz * a    # suma la posición relativa de la base en z
                    posiciones.append((x, y, z))
    return posiciones


# ------------------------------- 
#queremos calcular V(r_ij) para cada átomo i y el resto de los átomos son j
def lennardJones(r, sigma=2.3151, epsilon=0.167, n=12, m=6):
    return 4*epsilon*(((sigma)/(r))**n-((sigma)/(r))**m)


# -------------------------------       
# Deginimos una función para calcular la matriz de distancias entre todos los átomos.
def calcular_distancias(posiciones):
    diff = posiciones[:, None, :] - posiciones[None, :, :]  # Calculamos el vector diferencia entre cada par de átomos
    distancias = np.linalg.norm(diff, axis=2)               # Calculamos la norma (magnitud) de cada vector diferencia para obtener la distancia
    np.fill_diagonal(distancias, np.inf)                    # Rellenamos la diagonal con infinito para evitar el cálculo de la interacción de un átomo consigo mismo (r_ii = 0)
    return distancias

    
# Definimos una función para encontrar los átomos vecinos dentro de un radio de corte
def obtener_primeros_vecinos(distancias, r_corte = 3*sigma):
    mask = distancias < r_corte
    return mask

# Función para graficar el potencial
def Plot3D_colormap(cristal, V_map, FC=3):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.set_title(f'Potencial por átomo (${FC} \sigma$)')
    p = ax.scatter(cristal[:,0], cristal[:,1], cristal[:,2], s = 120, c=V_map, cmap = plt.cm.viridis)
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    fig.colorbar(p, ax=ax, label= '$Potencial Lennard-Jones$')
    plt.show()
    return


# ------------------------------- Energía total considerando superficies libres
def calcular_energia_total(posiciones, sigma=2.3151, epsilon=0.167,n=12, m=6, r_corte=3*2.3151):

    distancias = calcular_distancias(posiciones)         # Calculamos la matriz de distancias entre todos los átomos.
    mask = obtener_primeros_vecinos(distancias, r_corte) # Obtienemos la máscara de vecinos que están dentro del radio de corte
    Vmat = np.zeros_like(distancias)                     # Creamos una matriz de ceros para almacenar los potenciales
    Vmat[mask] = lennardJones(distancias[mask], sigma=2.3151, epsilon=0.167, n=12, m=6) # Calculamos el potencial de Lennard-Jones solo para los pares de átomos dentro del radio de corte

    # Suma los potenciales para cada átomo (axis=1) y divide por 2 para no contar dos veces cada interacción.
    E_por_atomo = np.sum(Vmat, axis=1)
    E_total = 0.5*np.sum(E_por_atomo)

    return E_total, E_por_atomo




# ------------------------------- Energía total considerando condiciones periódicas

def Lennard_Jones_derivada(r, sigma=2.3151, epsilon=0.167, n=12, m=6):
    r_inv = sigma / r
    return 4 * epsilon * (n * r_inv**(n) - m * r_inv**(m)) / r


def calcular_distancias_periodicas(posiciones, L):
    L = np.asarray(L)
    diff = posiciones[:, None, :] - posiciones[None, :, :]
    diff = diff - L * np.round(diff / L)
    distancias = np.linalg.norm(diff, axis=2)
    unit_vectors = np.zeros_like(diff)
    mask = np.isfinite(distancias)
    unit_vectors[mask] = diff[mask] / distancias[mask][:, None]
    return distancias, unit_vectors


def calcular_energia_total_periodica(posiciones, L, sigma=2.3151, epsilon=0.167, n=12, m=6, r_corte=None):
    if r_corte is None:
        r_corte = 3 * sigma
    r_corte = min(r_corte, 0.5 * min(L) - 1e-9) #cambiado
        
        
        
    N = len(posiciones)
    distancias, unit_vectors = calcular_distancias_periodicas(posiciones, L)

    Vmat = np.zeros((N, N))
    F_modulo = np.zeros((N, N))

    for i in range(N):
        for j in range(i + 1, N):
            r = distancias[i, j]
            if r < r_corte:
                V = lennardJones(r, sigma, epsilon, n, m)
                dVdr = Lennard_Jones_derivada(r, sigma=2.3151, epsilon=0.167, n=12, m=6)
                Vmat[i, j] = V
                Vmat[j, i] = V
                F_modulo[i, j] = dVdr
                F_modulo[j, i] = dVdr

    E_por_atomo = np.sum(Vmat, axis=1)
    E_total = 0.5 * np.sum(E_por_atomo)
    return E_total, E_por_atomo, F_modulo



def Plot3D_colormapminmax(cristal, V_map, FC=3):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.set_title(f'Potencial por átomo (${FC} \sigma$)')
    p = ax.scatter(cristal[:,0], cristal[:,1], cristal[:,2], s = 120, vmin=-2, vmax=0, c=V_map, cmap = plt.cm.viridis)
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    fig.colorbar(p, ax=ax, label= '$Potencial Lennard-Jones$')
    plt.show()
    return



#################################dia 14 octubre#################################

#################################parte 1: cálculo de la fuerza#################################


def calcular_distancias_vector(posiciones):
    diff = posiciones[:, None, :] - posiciones[None, :, :]  # Calculamos el vector diferencia entre cada par de átomos
    distancias = np.linalg.norm(diff, axis=2)               # Calculamos la norma (magnitud) de cada vector diferencia para obtener la distancia
    np.fill_diagonal(distancias, np.inf)                    # Rellenamos la diagonal con infinito para evitar el cálculo de la interacción de un átomo consigo mismo (r_ii = 0)
    #añadimos lo siguiente para obtener los vectores unitarios:
    unit_vectors = np.zeros_like(diff)
    
    N = len(posiciones)
    for i in range(N):
        for j in range(N):
            if i != j:  # Evita el caso de la diagonal
                unit_vectors[i, j] = diff[i, j] / distancias[i, j]
    
    return distancias, unit_vectors


def calcular_distancias_periodicas_vector(posiciones, L):
    L = np.asarray(L)
    diff = posiciones[:, None, :] - posiciones[None, :, :]   # r_i - r_j
    # Aplicar mínima imagen (vectorial)
    diff = diff - L * np.round(diff / L)
    distancias = np.linalg.norm(diff, axis=2)
    unit_vectors = np.zeros_like(diff)
    mask = np.isfinite(distancias)  # True para todos excepto la diagonal
    unit_vectors[mask] = diff[mask] / distancias[mask][:, None]
    return distancias, unit_vectors


def calcular_energia_total_vector(posiciones, sigma=2.3151, epsilon=0.167,n=12, m=6, r_corte=3*2.3151):
    distancias = calcular_distancias(posiciones)         # Calculamos la matriz de distancias entre todos los átomos.
    mask = obtener_primeros_vecinos(distancias, r_corte) # Obtienemos la máscara de vecinos que están dentro del radio de corte
    Vmat = np.zeros_like(distancias)                     # Creamos una matriz de ceros para almacenar los potenciales
    Vmat[mask] = lennardJones(distancias[mask], sigma=2.3151, epsilon=0.167, n=12, m=6) # Calculamos el potencial de Lennard-Jones solo para los pares de átomos dentro del radio de corte

    # Suma los potenciales para cada átomo (axis=1) y divide por 2 para no contar dos veces cada interacción.
    E_por_atomo = np.sum(Vmat, axis=1)
    E_total = 0.5*np.sum(E_por_atomo)
    
    #ahora añadimos el calcular la fuerza en vector
    F_modulo = np.zeros_like(distancias) 
    F_modulo[mask] = Lennard_Jones_derivada(distancias[mask], sigma=2.3151, epsilon=0.167, n=12, m=6)

    return E_total, E_por_atomo, F_modulo


def calcular_energia_total_periodica_vector(posiciones, L, sigma=2.3151, epsilon=0.167,n=12, m=6, r_corte=None):

    if r_corte is None:
        r_corte = 3 * sigma
    r_corte = min(r_corte, 0.5 * min(L) - 1e-9)    #cambiado

    N = len(posiciones)
    distancias, unit_vectors = calcular_distancias_periodicas(posiciones, L)# Calcula la matriz de distancias usando la función periódica
    Vmat = np.zeros((N, N))
    F_modulo = np.zeros((N, N))
    for i in range(N):
        for j in range(i + 1, N):
            r = distancias[i, j]
            if r < r_corte:
                V = lennardJones(r, sigma, epsilon, n, m)
                F = Lennard_Jones_derivada(r, sigma, epsilon, n, m)  # fuerza = -dV/dr
                Vmat[i, j] = V
                Vmat[j, i] = V
                F_modulo[i, j] = F
                F_modulo[j, i] = F

    # Energía por átomo y total
    E_por_atomo = np.sum(Vmat, axis=1)
    E_total = 0.5 * np.sum(E_por_atomo)
              
    return E_total, E_por_atomo, F_modulo


def vector_fuerza(F_modulo, unit_vectors):
    # Asegurar que no haya NaN ni inf en la matriz de fuerzas
    F_modulo = np.nan_to_num(F_modulo, nan=0.0, posinf=0.0, neginf=0.0)
    unit_vectors = np.nan_to_num(unit_vectors, nan=0.0, posinf=0.0, neginf=0.0)

    # Calcular la matriz de fuerzas vectoriales (N,N,3)
    F_mat = F_modulo[:, :, None] * unit_vectors

    # Fuerza neta sobre cada átomo (suma sobre j)
    F_neta = np.nansum(F_mat, axis=1)  # ignora posibles NaN residuales

    return F_mat, F_neta

    

#################################parte 2: termodinámica#################################
kB = 8.617333262*10**(-5)  # eV/K (si estás trabajando en eV y Å)

def generar_velocidades_iniciales(posiciones, T_deseada, masas):
    #Dentro de una función se generará una distribución de números aleatorios (numpy.random_sample) centrada en 0
    N = len(posiciones)

    v = np.random.random_sample((N, 3)) - 0.5

    #se restará la velocidad del centro de masas a cada velocidad para no empezar con un centro de masas en movimiento.
    correc_x = np.sum(v[:, 0]) / N
    correc_y = np.sum(v[:, 1]) / N
    correc_z = np.sum(v[:, 2]) / N

    v[:, 0] -= correc_x
    v[:, 1] -= correc_y
    v[:, 2] -= correc_z
    #Se calculará la temperatura de esta distribución random y se calculará el cociente necesario para transformar esa temperatura aleatoria deseada (cuidado con las unidades)
    T_inicial = calcT(v, masas)
    factor = np.sqrt(T_deseada / T_inicial)
    v *= factor
    
    return v, T_inicial, factor

def calcT(v, masas):
    N = len(masas)
    Ec = calcula_Ec(v, masas) 
    return (2 * Ec) / (3 * N * kB)

def calcula_Ec(v, masas):
    return 0.5 * np.sum(masas[:, None] * v**2)




################################# PARTE 3: DINÁMICA MOLECULAR (VERSIÓN ESTABLE) #################################

def dinamica_verlet(posiciones_iniciales, velocidades_iniciales, L, N_pasos, dt, masas, sigma, epsilon, n, m):

    N = len(posiciones_iniciales)
    r_corte = 0.5 * min(L) - 1e-9  # radio de corte (solo se usa para calcular fuerzas)

    posiciones_t = np.zeros((N_pasos, N, 3))
    velocidades_t = np.zeros((N_pasos, N, 3))
    energia_potencial_t = np.zeros(N_pasos)

    posiciones = posiciones_iniciales.copy()
    velocidades = velocidades_iniciales.copy()

    posiciones_t[0] = posiciones
    velocidades_t[0] = velocidades

    # fuerzas iniciales
    E_pot_actual, _, F_mod = calcular_energia_total_vector(posiciones, sigma, epsilon, n, m, r_corte)
    unit_vectors = calcular_distancias_vector(posiciones)[1]
    F_vect_actual = vector_fuerza(F_mod, unit_vectors)[1]
    energia_potencial_t[0] = E_pot_actual

    # bucle principal
    for paso in range(1, N_pasos):
        print(paso)
        aceleracion_actual = F_vect_actual / masas[:, None]
        posiciones = posiciones + velocidades * dt + 0.5 * aceleracion_actual * dt**2
        E_pot_actual, _, F_mod_nuevo = calcular_energia_total_vector(posiciones, sigma, epsilon, n, m, r_corte)
        unit_vectors_nuevo = calcular_distancias_vector(posiciones)[1]
        F_vect_nuevo = vector_fuerza(F_mod_nuevo, unit_vectors_nuevo)[1]

        aceleracion_nueva = F_vect_nuevo / masas[:, None]
        velocidades = velocidades + 0.5 * (aceleracion_actual + aceleracion_nueva) * dt

        # guardar resultados
        posiciones_t[paso] = posiciones
        velocidades_t[paso] = velocidades
        energia_potencial_t[paso] = E_pot_actual

        F_vect_actual = F_vect_nuevo

    return posiciones_t, velocidades_t, energia_potencial_t


def graficar_termica_periodica(energia_cinetica, energia_potencial, temperatura):
    pasos = np.arange(len(energia_cinetica))
    plt.figure(figsize=(14, 6))
    
    # --- Gráfico de Temperatura ---
    plt.subplot(1, 2, 1)
    plt.plot(pasos, temperatura, color='crimson')
    plt.xlabel('Paso de Simulación')
    plt.ylabel('Temperatura (K)')
    plt.title('Evolución de la Temperatura')
    plt.grid(True)

    # --- Gráfico de Energías (Conforme al enunciado) ---
    plt.subplot(1, 2, 2)
    
    # Energía potencial relativa a la inicial
    E_pot_rel = energia_potencial - energia_potencial[0]
    
    # Energía total relativa a la inicial
    #E_total_rel = (E_cinetica + E_potencial) - (E_cinetica[0] + E_potencial[0])
    # que es lo mismo que:
    E_total_rel = (energia_cinetica - energia_cinetica[0]) + E_pot_rel
    
    plt.plot(pasos, energia_cinetica, label='E Cinética (Absoluta)')
    plt.plot(pasos, E_pot_rel, label='E Potencial (respecto a inicial)')
    plt.plot(pasos, E_total_rel, label='E Total (respecto a inicial)', linewidth=2.5, color='black')
    
    plt.xlabel('Paso de Simulación')
    plt.ylabel('Energía (eV)')
    plt.title('Conservación de la Energía (relativa al paso 0)')
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.show()
    

def graficar_termica(energia_cinetica, energia_potencial, temperatura):
    pasos = np.arange(len(energia_cinetica))
    plt.figure(figsize=(14, 6))
    
    # --- Gráfico de Temperatura ---
    plt.subplot(1, 2, 1)
    plt.plot(pasos, temperatura, color='crimson')
    plt.xlabel('Paso de Simulación')
    plt.ylabel('Temperatura (K)')
    plt.title('Evolución de la Temperatura')
    plt.grid(True)

    # --- Gráfico de Energías ---
    plt.subplot(1, 2, 2)
    
    # Energía potencial relativa a la inicial
    E_pot_rel = energia_potencial - energia_potencial[0]
    
    # Energía total relativa a la inicial
    E_total_rel = (energia_cinetica - energia_cinetica[0]) + E_pot_rel
    
    plt.plot(pasos, energia_cinetica, label='E Cinética (absoluta)', color='royalblue')
    plt.plot(pasos, E_pot_rel, label='E Potencial (relativa)', color='darkorange')
    plt.plot(pasos, E_total_rel, label='E Total (relativa)', linewidth=2.5, color='black')
    
    plt.xlabel('Paso de Simulación')
    plt.ylabel('Energía (eV)')
    plt.title('Conservación de la Energía (relativa al paso 0)')
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.show()






def save_to_film(r,t,N, name = 'Film'):
    r = r.round(5)
    out = open(name+'.xyz', 'w')
    for i in range(len(t)//2): # Ajustado para guardar todos los pasos
        out.write(f'{N}\n')
        out.write(f'Atoms. Timestep: {int(t[2*i])}\n')
        for j in range(N):
            out.write(f'Cu      {r[i][j][0]}        {r[i][j][1]}        {r[i][j][2]}\n')
    out.close()
    return


def Plot3D_quiver(cristal, vectors, V_map, label_cb = '$Potencial Lennard-Jones$'):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.set_title('Potencial por atomo')
    p = ax.scatter(cristal[:,0], cristal[:,1], cristal[:,2], c = V_map, cmap = plt.cm.viridis)
    ax.quiver(cristal[:,0], cristal[:,1], cristal[:,2], vectors[:,0], vectors[:,1], vectors[:,2], color = 'deepskyblue', length=1, normalize=True)
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    fig.colorbar(p, ax=ax, label= label_cb)
    plt.show()
    return


nx = 5          #número de celdas en la dirección x
ny = 5          #número de celdas en la dirección y
nz = 5          #número de celdas en la dirección z
a = 3.603       #distancia entre puntos
sigma=2.3151
epsilon=0.167
n=12
m=6

amu_a_kg      = 1.66053906660e-27   # 1 uma en kg
eV_a_J        = 1.602176634e-19     # 1 eV en J
angstrom_a_m  = 1e-10               # 1 Å en m
posiciones_fcc = generar_red_fcc(nx, ny, nz, a)
posiciones_fcc = np.array(posiciones_fcc)
print(f"Numero total de posiciones generadas fcc: {len(posiciones_fcc)}")

# graficamos lo obtenido para el caso fcc
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.scatter(posiciones_fcc[:,0], posiciones_fcc[:,1], posiciones_fcc[:,2], s=40, c='red')
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')
ax.set_title('Red cúbica centrada en las caras (FCC)')
E_total, E_por_atomo = calcular_energia_total(posiciones_fcc, sigma, epsilon, r_corte=3*sigma)
Plot3D_colormap(posiciones_fcc, E_por_atomo, FC=3)
print("Energía total:", E_total)
print("Energía total/atomo:", E_total/len(posiciones_fcc))
print("Energía total/atomo:", np.mean(E_por_atomo))


L = [nx * a, ny * a, nz * a]
# Parámetros
N_pasos = 3000
dt = 1e-15  # 1 femtosegundo (10^-15 s)
Temp_deseada = 300  # Kelvin

# Cálculo de la masa en unidades de la simulación
masa_u_en_kg = 1.660539e-27
eV_en_J = 1.602176e-19
angstrom_en_m = 1e-10
masa_Cu_kg = 63.546 * masa_u_en_kg
masa_Cu = masa_Cu_kg * (angstrom_en_m**2) / eV_en_J
masas = np.full(len(posiciones_fcc), masa_Cu)

# (c) Inicializar velocidades
velocidades_iniciales = generar_velocidades_iniciales(posiciones_fcc, Temp_deseada, masas)[0]

# (f) Verificación de la temperatura inicial
T_verificada = calcT(velocidades_iniciales, masas)
print(f"Temperatura objetivo: {Temp_deseada} K")
print(f"Temperatura real tras inicialización: {T_verificada:.2f} K\n")

# Ejecutar la simulación
print("Iniciando simulación de dinámica molecular...")
pos_t, vel_t, Epot_t = dinamica_verlet(posiciones_fcc, velocidades_iniciales, L, N_pasos, dt, masas, sigma, epsilon, n, m)
print("Simulación finalizada.")

# (e) Calcular las otras propiedades para cada paso
Ecin_t = np.array([calcula_Ec(v, masas) for v in vel_t])
Temp_t = np.array([calcT(v, masas) for v in vel_t])

# (f, g) Representar los resultados
graficar_termica(Ecin_t, Epot_t, Temp_t)
E_total, E_por_atomo, F_modulo = calcular_energia_total_periodica_vector(posiciones_fcc, L)
unit_vectors = calcular_distancias_periodicas_vector(posiciones_fcc, L)[1]
F_vect, F_neta= vector_fuerza(F_modulo, unit_vectors)
Plot3D_quiver(posiciones_fcc, F_neta, E_por_atomo)

# Guardar trayectoria para visualización
save_to_film(pos_t, np.arange(N_pasos), len(posiciones_fcc), name='Film_Cu_establepasos500T800')
print("Trayectoria guardada en 'Film_Cu_estable.xyz'")







"""
# gráficas dependencia temperatura
def graficar_comparativa_termica(lista_resultados, lista_temperaturas):

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(18, 8), sharex=True)
    pasos = np.arange(len(lista_resultados[0]['Ecin_t']))
    estilos_por_temp = {300: {'linestyle': '-', 'label_suffix': '1600 K (Líq.)'},
        100: {'linestyle': '--', 'label_suffix': '800 K (Sól.)'},
        50: {'linestyle': ':', 'label_suffix': '100 K (Sól.)'}}
    # El color representa el TIPO DE ENERGÍA
    colores_por_energia = {'Epot': 'blue','Ecin': 'orange','Etot': 'black'}

    temp_plot_colors = {'100 K (Sól.)': 'blue', '800 K (Sól.)': 'orange', '1600 K (Líq.)': 'red'}
    for i, res in enumerate(lista_resultados):
        temp_actual = lista_temperaturas[i]
        style = estilos_por_temp.get(temp_actual)
        label = style['label_suffix']
        ax1.plot(pasos, res['Temp_t'], color=temp_plot_colors[label], linestyle=style['linestyle'], label=label)

    ax1.set_title('Evolución Comparativa de la Temperatura')
    ax1.set_xlabel('Paso de Simulación')
    ax1.set_ylabel('Temperatura (K)')
    ax1.grid(True)
    ax1.legend()

    for i, res in enumerate(lista_resultados):
        temp_actual = lista_temperaturas[i]
        style = estilos_por_temp.get(temp_actual)

        Ecin_t = res['Ecin_t']
        Epot_t = res['Epot_t']
        E_pot_rel = Epot_t - Epot_t[0]
        E_total_rel = (Ecin_t - Ecin_t[0]) + E_pot_rel
        ax2.plot(pasos, E_pot_rel, color=colores_por_energia['Epot'], linestyle=style['linestyle'], 
                 label=f'E Potencial ({style["label_suffix"]})')
        ax2.plot(pasos, Ecin_t, color=colores_por_energia['Ecin'], linestyle=style['linestyle'], 
                 label=f'E Cinética ({style["label_suffix"]})')
        ax2.plot(pasos, E_total_rel, color=colores_por_energia['Etot'], linestyle=style['linestyle'], linewidth=2, label=f'E Total ({style["label_suffix"]})')

    ax2.set_title('Evolución Energética Comparativa')
    ax2.set_xlabel('Paso de Simulación')
    ax2.set_ylabel('Energía (eV)')
    ax2.grid(True)
    ax2.legend(bbox_to_anchor=(1.04, 1), loc="upper left")
    
    plt.suptitle('Análisis Comparativo del Sistema a Diferentes Temperaturas', fontsize=16)
    fig.tight_layout(rect=[0, 0, 0.85, 0.95]) 
    plt.show()
"""
# gráficas dependencia temperatura
def graficar_comparativa_termica(lista_resultados, lista_temperaturas):

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(18, 8), sharex=True)
    pasos = np.arange(len(lista_resultados[0]['Ecin_t']))

    # Diccionario de estilos según la temperatura simulada
    estilos_por_temp = {
        300: {'linestyle': '-', 'label_suffix': '300 K (Líq.)'},
        100: {'linestyle': '--', 'label_suffix': '100 K (Sól.)'},
        50: {'linestyle': ':', 'label_suffix': '50 K (Sól.)'}
    }

    # Colores para cada tipo de energía
    colores_por_energia = {
        'Epot': 'blue',
        'Ecin': 'orange',
        'Etot': 'black'
    }

    # Colores usados para las temperaturas en el gráfico de temperatura
    temp_plot_colors = {
        '50 K (Sól.)': 'blue',
        '100 K (Sól.)': 'orange',
        '300 K (Líq.)': 'red'
    }

    # --- Gráfico 1: Evolución de la temperatura ---
    for i, res in enumerate(lista_resultados):
        temp_actual = int(lista_temperaturas[i])  # aseguramos tipo int
        style = estilos_por_temp.get(
            temp_actual,
            {'linestyle': '-', 'label_suffix': f'{temp_actual} K'}  # estilo por defecto
        )

        label = style['label_suffix']
        color = temp_plot_colors.get(label, 'gray')  # color por defecto si falta
        ax1.plot(pasos, res['Temp_t'], color=color, linestyle=style['linestyle'], label=label)

    ax1.set_title('Evolución Comparativa de la Temperatura')
    ax1.set_xlabel('Paso de Simulación')
    ax1.set_ylabel('Temperatura (K)')
    ax1.grid(True)
    ax1.legend()

    # --- Gráfico 2: Energías comparativas ---
    for i, res in enumerate(lista_resultados):
        temp_actual = int(lista_temperaturas[i])
        style = estilos_por_temp.get(
            temp_actual,
            {'linestyle': '-', 'label_suffix': f'{temp_actual} K'}
        )

        Ecin_t = res['Ecin_t']
        Epot_t = res['Epot_t']
        E_pot_rel = Epot_t - Epot_t[0]
        E_total_rel = (Ecin_t - Ecin_t[0]) + E_pot_rel

        ax2.plot(pasos, E_pot_rel, color=colores_por_energia['Epot'],
                 linestyle=style['linestyle'], label=f'E Potencial ({style["label_suffix"]})')
        ax2.plot(pasos, Ecin_t, color=colores_por_energia['Ecin'],
                 linestyle=style['linestyle'], label=f'E Cinética ({style["label_suffix"]})')
        ax2.plot(pasos, E_total_rel, color=colores_por_energia['Etot'],
                 linestyle=style['linestyle'], linewidth=2, label=f'E Total ({style["label_suffix"]})')

    ax2.set_title('Evolución Energética Comparativa')
    ax2.set_xlabel('Paso de Simulación')
    ax2.set_ylabel('Energía (eV)')
    ax2.grid(True)
    ax2.legend(bbox_to_anchor=(1.04, 1), loc="upper left")

    plt.suptitle('Análisis Comparativo del Sistema a Diferentes Temperaturas', fontsize=16)
    fig.tight_layout(rect=[0, 0, 0.85, 0.95])
    plt.show()



def main():

    print("--- Ejecutando Experimento 2: Influencia de la Temperatura ---")
    nx, ny, nz, a = 5, 5, 5, 3.603
    sigma, epsilon, n, m = 2.3151, 0.167, 12, 6
    N_pasos = 3000
    dt = 1e-15

    posiciones_fcc = np.array(generar_red_fcc(nx, ny, nz, a))
    L = [nx * a, ny * a, nz * a]
    
    masa_u_en_kg = 1.660539e-27
    eV_en_J = 1.602176e-19
    masa_Cu_kg = 63.546 * masa_u_en_kg
    masa_Cu = masa_Cu_kg * (1e-10**2) / eV_en_J
    masas = np.full(len(posiciones_fcc), masa_Cu)
    
    temperaturas_a_simular = [300, 100, 50]
    resultados_completos = []

    for temp in temperaturas_a_simular:
        print(f"\n--- Iniciando simulación a {temp} K ---")
        
        velocidades_iniciales= generar_velocidades_iniciales(posiciones_fcc, temp, masas)[0]
        
        pos_t, vel_t, Epot_t = dinamica_verlet(posiciones_fcc, velocidades_iniciales, L, N_pasos, dt, masas, sigma, epsilon, n, m)
        
        Ecin_t = np.array([calcula_Ec(v, masas) for v in vel_t])
        Temp_t = np.array([calcT(v, masas) for v in vel_t])
        
        resultados_completos.append({'Ecin_t': Ecin_t, 'Epot_t': Epot_t,'Temp_t': Temp_t})
        
        # Guardar el archivo .xyz para esta temperatura
        nombre_archivo = f'Film_Cu_{temp}K.xyz'
        save_to_film(pos_t, np.arange(N_pasos), len(posiciones_fcc), name=nombre_archivo)
        print(f"Simulación a {temp} K finalizada. Trayectoria guardada en '{nombre_archivo}'")


    print("\n--- Generando Gráficos Comparativos ---")
    graficar_comparativa_termica(resultados_completos, temperaturas_a_simular)


if __name__ == "__main__":
    main()



# --- CÓDIGO PARA EL EXPERIMENTO 3: t ---

def graficar_estabilidad_dt(lista_resultados):

    plt.figure(figsize=(10, 6))
    
    estilos = {
        0.5e-15: {'color': 'purple', 'linestyle': '-'},
        1e-15: {'color': 'green', 'linestyle': '--'},
        10e-15: {'color': 'red', 'linestyle': ':'}
    }

    for res in lista_resultados:
        dt = res['dt']
        style = estilos.get(dt)
        pasos = np.arange(len(res['Ecin_t']))
        
        # Calcular energía total relativa
        energia_pot_rel = res['Epot_t'] - res['Epot_t'][0]
        energia_total_rel = (res['Ecin_t'] - res['Ecin_t'][0]) + energia_pot_rel
        
        if np.any(np.isnan(energia_total_rel)):
            print(f"Advertencia: Se encontraron NaNs en la simulación con dt = {dt*1e15:.1f} fs.")
            energia_total_rel[np.isnan(energia_total_rel)] = np.nan
        
        plt.plot(pasos, energia_total_rel, 
                 label=f'$\\Delta t$ = {dt*1e15:.1f} fs', 
                 color=style['color'], 
                 linestyle=style['linestyle'])

    plt.xlabel('Paso de Simulación')
    plt.ylabel('Energía Total Relativa (eV)')
    plt.title('Estabilidad del Integrador de Verlet vs. Tamaño de Paso')
    plt.legend()
    plt.grid(True)
    plt.ylim(-1, 10) 
    plt.show()


def main_experimento_dt():

    print("--- Ejecutando Experimento 3: Estabilidad del Integrador ---")
    
    # Parámetros
    nx, ny, nz, a = 5, 5, 5, 3.603
    sigma, epsilon, n, m = 2.3151, 0.167, 12, 6
    N_pasos = 3000
    Temp_deseada = 100

    posiciones_fcc = np.array(generar_red_fcc(nx, ny, nz, a))
    L = np.array([nx * a, ny * a, nz * a])
    
    masa_u_en_kg = 1.660539e-27; eV_en_J = 1.602176e-19
    masa_Cu_kg = 63.546 * masa_u_en_kg
    masa_Cu = masa_Cu_kg * (1e-10**2) / eV_en_J
    masas = np.full(len(posiciones_fcc), masa_Cu)
    
    # Bucle de simulaciones
    pasos_de_tiempo_a_probar = [0.5e-15, 1e-15, 10e-15]
    resultados = []

    for dt in pasos_de_tiempo_a_probar:
        print(f"\n--- Iniciando simulación con dt = {dt*1e15:.1f} fs ---")
        velocidades_iniciales, _, _ = generar_velocidades_iniciales(posiciones_fcc, Temp_deseada, masas)
        pos_t, vel_t, Epot_t = dinamica_verlet(posiciones_fcc, velocidades_iniciales, L, N_pasos, dt, masas, sigma, epsilon, n, m)
        Ecin_t = np.array([calcula_Ec(v, masas) for v in vel_t])
        
        resultados.append({'dt': dt, 'Ecin_t': Ecin_t, 'Epot_t': Epot_t})
        print(f"Simulación con dt = {dt*1e15:.1f} fs finalizada.")

    # Visualización
    print("\n--- Generando Gráfico de Estabilidad ---")
    graficar_estabilidad_dt(resultados)



def graficar_dependencia_tamano_con_temperatura(lista_resultados):

    fig, axs = plt.subplots(2, 1, figsize=(15, 12), sharex=True)
    pasos = np.arange(len(lista_resultados[0]['Ecin_t']))

    colores = {
        32: 'blue',
        108: 'green',
        256: 'orange',
        500: 'red'
    }

    estilos_energia = {
        'Epot': '--', # Potencial
        'Ecin': ':',  # Cinética
        'Etot': '-'   # Total
    }

    for res in lista_resultados:
        N = res['N']
        color = colores.get(N, 'black')
        axs[0].plot(pasos, res['Temp_t'], color=color, label=f'N = {N} átomos')

    axs[0].set_title('Evolución de la Temperatura vs. Pasos de Simulación', fontsize=14)
    axs[0].set_ylabel('Temperatura (K)')
    axs[0].grid(True)
    axs[0].legend()

    for res in lista_resultados:
        N = res['N']
        color = colores.get(N, 'black')
        
        # Calcular energías
        Ecin_t = res['Ecin_t']
        Epot_t = res['Epot_t']
        Etot_t = Ecin_t + Epot_t
        
        # Graficar cada tipo de energía
        label_base = f'N={N}'
        axs[1].plot(pasos, Epot_t, color=color, linestyle=estilos_energia['Epot'], label=f'E Potencial ({label_base})')
        axs[1].plot(pasos, Ecin_t, color=color, linestyle=estilos_energia['Ecin'], label=f'E Cinética ({label_base})')
        axs[1].plot(pasos, Etot_t, color=color, linestyle=estilos_energia['Etot'], linewidth=2.5, label=f'E Total ({label_base})')

    axs[1].set_title('Evolución de las Energías vs. Pasos de Simulación', fontsize=14)
    axs[1].set_xlabel('Número de Pasos')
    axs[1].set_ylabel('Energía (eV)')
    axs[1].grid(True)
    axs[1].legend(bbox_to_anchor=(1.04, 1), loc="upper left")

    plt.suptitle('Análisis Termodinámico en Función del Tamaño del Sistema', fontsize=18)
    fig.tight_layout(rect=[0, 0, 0.9, 0.95]) 
    plt.show()



def main_experimento_tamano():
    print("--- Ejecutando Experimento: Análisis Termodinámico vs. Tamaño ---")
    
    # Parámetros base
    a = 3.603; sigma, epsilon, n, m = 2.3151, 0.167, 12, 6
    N_pasos = 3000
    dt = 1e-15
    Temp_deseada = 100 

    tamanos_a_simular = [2, 3, 4] 
    resultados_completos = []

    for tamano in tamanos_a_simular:
        nx = ny = nz = tamano
        N = 4 * (tamano**3)
        print(f"\n--- Iniciando simulación para sistema {nx}x{ny}x{nz} ({N} átomos) ---")
        
        posiciones_fcc = np.array(generar_red_fcc(nx, ny, nz, a))
        L = np.array([nx * a, ny * a, nz * a])
        
        masa_u_en_kg = 1.660539e-27; eV_en_J = 1.602176e-19
        masa_Cu_kg = 63.546 * masa_u_en_kg
        masa_Cu = masa_Cu_kg * (1e-10**2) / eV_en_J
        masas = np.full(len(posiciones_fcc), masa_Cu)

        velocidades_iniciales, _, _ = generar_velocidades_iniciales(posiciones_fcc, Temp_deseada, masas)
        pos_t, vel_t, Epot_t = dinamica_verlet(posiciones_fcc, velocidades_iniciales, L, N_pasos, dt, masas, sigma, epsilon, n, m)
        
        print("Calculando propiedades termodinámicas...")
        Ecin_t = np.array([calcula_Ec(v, masas) for v in vel_t])
        Temp_t = np.array([calcT(v, masas) for v in vel_t])
        
        resultados_completos.append({
            'N': N,
            'Temp_t': Temp_t,
            'Ecin_t': Ecin_t,
            'Epot_t': Epot_t
        })
        print(f"Simulación para N={N} finalizada.")

    print("\n--- Generando Gráficos de Dependencia con el Tamaño ---")
    graficar_dependencia_tamano_con_temperatura(resultados_completos)


if __name__ == "__main__":
    main_experimento_dt()
    main_experimento_tamano()
    


