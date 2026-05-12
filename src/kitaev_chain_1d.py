"""
kitaev_chain_1d.py

Numerical study of Majorana zero modes in one-dimensional superconducting systems.
This file was prepared for a GitHub repository from the original practice scripts.
The numerical model is intentionally explicit: dense BdG matrices are built in real space
so that the Hamiltonian, observables, and generated figures can be read directly.

For publication-quality figures, increase the site counts carefully. Dense diagonalization
scales poorly with system size.
"""

import numpy as np
import scipy.linalg as la
import matplotlib.pyplot as plt

###############################################################################
#################### Parámetros físicos (cadena de Kitaev) ####################
###############################################################################

# Tomamos el hopping como unidade de energía para que este normalizado
numero_sitios = 180     # demo-safe default; increase for final high-resolution figures          # Sitios de la cadena finita
hopping = 1.0                   # Hopping entre sitios vecinos
potencial_quimico = 0.0         # Potencial químico mu
pairing_p_wave = 0.5            # Pairing superconductivo p-wave entre vecinos

# Parámetros para los barridos
numero_energias = 400           # Número de energías para la LDOS (valores que vamos a resolver)
ensanchamiento_ldos = 0.04      # Anchura lorentziana artificial para visualizar la LDOS

# Casos representativos.
mu_trivial = 3.0                # Caso trivial: |mu| > 2t
mu_transicion = 2.0             # Transición: |mu| = 2t
mu_topologico = 0.0             # Caso topológico: |mu| < 2t

# Pairings alternativos para comparar cómo cambia la localización
pairing_debil = 0.2
pairing_medio = 0.5
pairing_fuerte = 0.9


###############################################################################
############################## Hamiltoniano BdG ##############################
###############################################################################

def construir_matrices_normal_y_pairing(numero_sitios, hopping, potencial_quimico, pairing_p_wave):
    # h_{ij}: parte potencial químico y hopping
    matriz_normal = np.zeros((numero_sitios, numero_sitios), dtype=complex)

    # Matriz Delta_{ij} (pairing): conecta partículas con huecos
    matriz_pairing = np.zeros((numero_sitios, numero_sitios), dtype=complex)

    # Término onsite: h_{jj} = -mu
    for j in range(numero_sitios):
        matriz_normal[j, j] = -potencial_quimico

    # Términos entre primeros vecinos
    for j in range(numero_sitios - 1):

        # Hopping normal: -t c_j^dagger c_{j+1} + h.c.
        matriz_normal[j, j + 1] = -hopping
        matriz_normal[j + 1, j] = -hopping

        # Pairing p-wave: debe ser antisimétrico: Delta_{j,j+1} = +Delta y Delta_{j+1,j} = -Delta
        matriz_pairing[j, j + 1] = pairing_p_wave
        matriz_pairing[j + 1, j] = -pairing_p_wave

    return matriz_normal, matriz_pairing


def construir_hamiltoniano_bdg(numero_sitios, hopping, potencial_quimico, pairing_p_wave):
    # Construimos las matrices N x N de la parte normal y superconductiva
    matriz_normal, matriz_pairing = construir_matrices_normal_y_pairing(numero_sitios, hopping, potencial_quimico, pairing_p_wave,)

    # Construimos natriz completa con los bloques
    hamiltoniano_bdg = np.block([[matriz_normal, matriz_pairing],[-matriz_pairing.conj(), -matriz_normal.conj()],])

    # Simetrizamos para eliminar posibles errores numéricos de redondeo y estar seguros que es hermítica
    hamiltoniano_bdg = 0.5 * (hamiltoniano_bdg + hamiltoniano_bdg.conj().T)

    return hamiltoniano_bdg


def diagonalizar_cadena(numero_sitios, hopping, potencial_quimico, pairing_p_wave):
    # Construimos el Hamiltoniano BdG
    hamiltoniano = construir_hamiltoniano_bdg(numero_sitios, hopping, potencial_quimico, pairing_p_wave,)

    # Diagonalizamos la matriz hermítica
    energias, vectores = la.eigh(hamiltoniano)

    # Ordenamos las energías de menor a mayor
    orden = np.argsort(energias)
    energias = energias[orden]
    vectores = vectores[:, orden]

    return energias, vectores, hamiltoniano


###############################################################################
############################# Observables físicos #############################
###############################################################################

def separar_electron_hueco(vector_estado, numero_sitios): #separamos componentes electrón y hueco
    # Consideramos los primeros N componentes como electrones
    componente_electron = vector_estado[:numero_sitios]

    # Y los últimos N componentes como huecos
    componente_hueco = vector_estado[numero_sitios:]

    return componente_electron, componente_hueco


def densidad_estado(vector_estado, numero_sitios): #calculamos densidad de estados
    # Separamos el estado en sus componentes de partícula y hueco
    componente_electron, componente_hueco = separar_electron_hueco(vector_estado, numero_sitios)

    # Aplicamos ecuación de densidad total BdG
    densidad = np.abs(componente_electron)**2 + np.abs(componente_hueco)**2

    # Normalizamos para que la suma total sea 1
    densidad = densidad / max(np.sum(densidad), 10**(-30))

    return densidad


def densidad_electron_y_hueco(vector_estado, numero_sitios): #calculamos densidad de electron y hueco
    componente_electron, componente_hueco = separar_electron_hueco(vector_estado, numero_sitios)

    norma = max(np.sum(np.abs(componente_electron)**2 + np.abs(componente_hueco)**2), 10**(-30))

    return componente_electron**2 / norma, componente_hueco**2 / norma


def indice_estado_mas_cercano_a_cero(energias): #calculamos energía cercana a 0
    return int(np.argmin(np.abs(energias))) #argmin devuelve el índice del valor mínimo, en este caso el valor mínimo de |E|, que es el estado más cercano a energía cero

def energia_minima_absoluta(energias): # Devuelve min |E|, útil para detectar estados casi cero
    return float(np.min(np.abs(energias))) # Devuelve el valor mínimo de |E|, que es la energía del estado más cercano a cero


def es_topologico(hopping, potencial_quimico, pairing_p_wave): #comprobamos si cumple la condición topológica
    return bool(abs(potencial_quimico) < 2.0 * abs(hopping) and abs(pairing_p_wave) > 0.0) #bool devuelve verdadero o falso según se cumple |mu| < 2|t| y que el pairing no sea cero


def transformar_particula_hueco(vector_estado, numero_sitios): #aplicamos simetría partícula-hueco a la base
    componente_electron, componente_hueco = separar_electron_hueco(vector_estado, numero_sitios)

    return np.concatenate([componente_hueco.conj(), componente_electron.conj()])


def alinear_fase(estado_referencia, estado_objetivo): #elimina fase global entre vectores
    solapamiento = np.vdot(estado_referencia, estado_objetivo)  # vdot calcula el producto escalar complejo entre dos vectores, devuelve un número complejo que representa el solapamiento entre ambos estados
    if abs(solapamiento) < 10**(-30):
        return estado_objetivo.copy()                           # Si el solapamiento es muy pequeño, no podemos alinear la fase, así que devolvemos el estado sin modificar
    return estado_objetivo * np.exp(-1j * np.angle(solapamiento))


def construir_majoranas_desde_pareja(energias, vectores, numero_sitios): #construir majoranas con estados cercanos a cero
    # Obtenemos índices de energía positiva y negativa
    indices_positivos = np.where(energias >= 0.0)[0]
    indices_negativos = np.where(energias < 0.0)[0]

    # Cogemos los estados cerca de cero positivo y negativo por separado
    indice_positivo = indices_positivos[np.argmin(np.abs(energias[indices_positivos]))]
    indice_negativo = indices_negativos[np.argmin(np.abs(energias[indices_negativos] + energias[indice_positivo]))]

    # Obtenemos los estados de las energías cercanas a cero
    estado_positivo = vectores[:, indice_positivo]
    estado_negativo = vectores[:, indice_negativo]

    # Aplicamos PHS al estado positivo
    estado_phs = transformar_particula_hueco(estado_positivo, numero_sitios)

    # Alineamos la fase del estado negativo con el compañero PHS
    estado_negativo = alinear_fase(estado_phs, estado_negativo)

    # Combinaciones tipo Majorana
    majorana_1 = (estado_positivo + estado_negativo) / np.sqrt(2.0)
    majorana_2 = -1j * (estado_positivo - estado_negativo) / np.sqrt(2.0)

    # Densidades de las dos combinaciones
    densidad_1 = densidad_estado(majorana_1, numero_sitios)
    densidad_2 = densidad_estado(majorana_2, numero_sitios)

    # Decidimos cuál vive más a la izquierda.
    mitad = numero_sitios // 2
    peso_izquierdo_1 = np.sum(densidad_1[:mitad])
    peso_izquierdo_2 = np.sum(densidad_2[:mitad])

    if peso_izquierdo_1 >= peso_izquierdo_2:
        densidad_izquierda = densidad_1
        densidad_derecha = densidad_2
    else:
        densidad_izquierda = densidad_2
        densidad_derecha = densidad_1

    return densidad_izquierda, densidad_derecha


def lorentziana(energia, centro, anchura): #Función lorentziana usada para suavizar la densidad local de estados
    return anchura / (np.pi * ((energia - centro)**2 + anchura**2))


def calcular_ldos(numero_sitios, energias, vectores, energias_ldos, anchura): #Calcula la LDOS(E,j) sumando contribuciones lorentzianas de los autestados
    ldos = np.zeros((len(energias_ldos), numero_sitios), dtype=float) 

    for n, energia_n in enumerate(energias):
        densidad_n = densidad_estado(vectores[:, n], numero_sitios)                 # Densidad espacial del autestado n
        for i, energia in enumerate(energias_ldos): # Sumamos su contribución a cada energía de la malla
            ldos[i, :] += densidad_n * lorentziana(energia, energia_n, anchura)
    return ldos


def bandas_bulk_kitaev(hopping, potencial_quimico, pairing_p_wave, numero_k=500): #Calcula las bandas bulk de la cadena de Kitaev infinita:E(k) = pm sqrt[(-mu - 2t cos k)^2 + (2 Delta sin k)^2]
    k = np.linspace(-np.pi, np.pi, numero_k)                 # Vector de momentos k en el intervalo [-pi, pi]

    xi_k = -potencial_quimico - 2.0 * hopping * np.cos(k)
    delta_k = 2.0 * pairing_p_wave * np.sin(k)               # Pairing p-wave en espacio de momentos

    energia_positiva = np.sqrt(xi_k**2 + np.abs(delta_k)**2) # Banda positiva E(k) = +sqrt[xi(k)^2 + |Delta(k)|^2]
    energia_negativa = -energia_positiva

    return k, energia_negativa, energia_positiva


###############################################################################
################################## Barridos ###################################
###############################################################################

def barrer_potencial_quimico(numero_sitios, hopping, pairing_p_wave, valores_mu): # Barrido de potencial químico para detectar modos de borde
    todas_las_energias = [] # Guardamos el espectro completo para cada valor de mu
    energias_minimas = [] # Guardamos solo la energía más cercana a cero para cada mu, que es el indicador de modos de borde

    for mu in valores_mu:
        # Diagonalizamos la cadena para este valor de mu
        energias, vectores, hamiltoniano = diagonalizar_cadena(numero_sitios,hopping, mu, pairing_p_wave,)
        # Guardamos el espectro y la mínima energía absoluta
        todas_las_energias.append(energias)
        energias_minimas.append(energia_minima_absoluta(energias))

    return np.array(todas_las_energias), np.array(energias_minimas)


def barrer_longitud_cadena(valores_longitud, hopping, potencial_quimico, pairing_p_wave): #Barre la longitud N para estudiar el splitting de los modos de Majorana
    energias_minimas = []

    for N in valores_longitud:
        # Diagonalizamos una cadena de longitud N
        energias, vectores, hamiltoniano = diagonalizar_cadena(int(N),hopping,potencial_quimico,pairing_p_wave,)
        # Guardamos la energía más cercana a cero
        energias_minimas.append(energia_minima_absoluta(energias))
    return np.array(energias_minimas)


def calcular_diagrama_fase_mu_delta(hopping, valores_mu, valores_delta): #Calcula un diagrama de fase ideal en el plano (mu, Delta)
    diagrama = np.zeros((len(valores_delta), len(valores_mu)))

    for i, Delta in enumerate(valores_delta):
        for j, mu in enumerate(valores_mu):
            diagrama[i, j] = 1.0 if es_topologico(hopping, mu, Delta) else 0.0  # Marcamos con 1 la fase topológica y con 0 la fase trivial

    return diagrama # Devuelve una matriz 2D con valores 1 para la fase topológica y 0 para la fase trivial, que se puede graficar como un mapa de colores en el plano (mu, Delta)


###############################################################################
################################## Graficar ###################################
###############################################################################

def preparar_estilo_graficas(): #Estilo sencillo para que las figuras sean claras en un informe
    plt.style.use("default")
    plt.rcParams.update({"figure.figsize": (9, 5.5),"axes.grid": True,"grid.alpha": 0.25,"font.size": 11,"axes.titlesize": 13,"axes.labelsize": 12,"legend.fontsize": 10,})



def graficar_bandas_bulk_tres_regimenes(hopping, pairing_p_wave):   #Grafica las bandas bulk para tres casos:trivial, transición y topológico
    casos = [("trivial", mu_trivial), ("transicion", mu_transicion), ("topologico", mu_topologico), ]

    figura, ejes = plt.subplots(1, 3, figsize=(13, 4), sharey=True) #Creamos una figura con 3 subplots 
    for eje, (nombre, mu) in zip(ejes, casos):
        k, energia_negativa, energia_positiva = bandas_bulk_kitaev(hopping, mu, pairing_p_wave,) # Calculamos las bandas bulk para este caso
        eje.plot(k / np.pi, energia_positiva, linewidth=2.0)        # Graficamos la banda positiva E(k) frente a k/pi
        eje.plot(k / np.pi, energia_negativa, linewidth=2.0)        # Graficamos la banda negativa E(k) frente a k/pi
        eje.axhline(0.0, color="black", linewidth=1.0, alpha=0.5)
        eje.set_title(f"{nombre}: $\\mu={mu}$")
        eje.set_xlabel(r"$k/\pi$")

    ejes[0].set_ylabel(r"Energía $E/t$")
    figura.suptitle("Bandas bulk de la cadena de Kitaev")
    figura.tight_layout()

    return figura



def graficar_espectro_finito_vs_mu(numero_sitios, hopping, pairing_p_wave): #Grafica el espectro finito frente a mu
    valores_mu = np.linspace(-4.0, 4.0, 121)
    espectros, energias_minimas = barrer_potencial_quimico( numero_sitios, hopping, pairing_p_wave,valores_mu,) # Barrido de mu para obtener espectros y energías mínimas

    figura, eje = plt.subplots()
    for i, mu in enumerate(valores_mu):
        energias_cerca_cero = espectros[i][np.abs(espectros[i]) < 1.5] # Filtramos solo las energías cercanas a cero para este valor de mu
        eje.scatter(np.full_like(energias_cerca_cero, mu), energias_cerca_cero, s=8, color="black", alpha=0.65,) # Graficamos estas energías como puntos

    eje.axvline(-2.0 * hopping, linestyle="--", linewidth=2.0, label=r"$\mu=-2t$") # Marcamos la transición esperada en mu=-2t
    eje.axvline(+2.0 * hopping, linestyle="--", linewidth=2.0, label=r"$\mu=+2t$") # Marcamos la transición esperada en mu=+2t
    eje.axhline(0.0, color="black", linewidth=1.0, alpha=0.5)
    eje.set_title("Espectro de la cadena finita frente al potencial químico")
    eje.set_xlabel(r"Potencial químico $\mu/t$")
    eje.set_ylabel(r"Energía $E/t$")
    eje.legend()

    figura.tight_layout()

    return figura



def graficar_energia_minima_vs_mu(numero_sitios, hopping, pairing_p_wave):         #Grafica min|E| frente a mu
    valores_mu = np.linspace(-4.0, 4.0, 121)
    espectros, energias_minimas = barrer_potencial_quimico(numero_sitios,hopping,pairing_p_wave, valores_mu,) # Barrido de mu para obtener espectros y energías mínimas

    figura, eje = plt.subplots()
    eje.plot(valores_mu, energias_minimas, linewidth=2.2)                          # Graficamos min|E| frente a mu
    eje.axvline(-2.0 * hopping, linestyle="--", linewidth=2.0, label=r"$\mu=-2t$") # Marcamos la transición esperada en mu=-2t
    eje.axvline(+2.0 * hopping, linestyle="--", linewidth=2.0, label=r"$\mu=+2t$") # Marcamos la transición esperada en mu=+2t
    eje.set_title(r"Energía más cercana a cero: $\min |E|$")
    eje.set_xlabel(r"Potencial químico $\mu/t$")
    eje.set_ylabel(r"$\min |E|/t$")
    eje.legend()
    figura.tight_layout()

    return figura



def graficar_ldos_tres_regimenes(numero_sitios, hopping, pairing_p_wave): #Grafica LDOS(E,j) para régimen trivial, transición y topológico
    casos = [("trivial", mu_trivial),("transicion", mu_transicion),("topologico", mu_topologico),]
    energias_ldos = np.linspace(-2.0, 2.0, numero_energias)
    figura, ejes = plt.subplots(1, 3, figsize=(13, 4), sharey=True)

    for eje, (nombre, mu) in zip(ejes, casos):
        energias, vectores, hamiltoniano = diagonalizar_cadena(numero_sitios,hopping,mu,pairing_p_wave,)
        ldos = calcular_ldos(numero_sitios,energias,vectores,energias_ldos, ensanchamiento_ldos,)
        imagen = eje.imshow(ldos,aspect="auto",origin="lower",extent=[0, numero_sitios - 1, energias_ldos[0], energias_ldos[-1]], cmap="magma",) # Graficamos la LDOS como una imagen con el mapa de colores "magma"
        eje.axhline(0.0, color="white", linewidth=1.0, alpha=0.8) # Marcamos la energía cero con una línea horizontal blanca
        eje.set_title(f"{nombre}: $\\mu={mu}$")
        eje.set_xlabel("Sitio j")
    ejes[0].set_ylabel(r"Energía $E/t$")
    figura.suptitle(r"LDOS$(E,j)$ de la cadena finita")
    figura.colorbar(imagen, ax=ejes.ravel().tolist(), label="LDOS")
    figura.tight_layout()

    return figura



def graficar_modo_cero_topologico_y_trivial(numero_sitios, hopping, pairing_p_wave): #Compara la densidad del estado más cercano a cero en un caso trivial y uno topológico
    casos = [ ("trivial", mu_trivial),("topologico", mu_topologico),]

    figura, ejes = plt.subplots(1, 2, figsize=(12, 4), sharey=True)
    for eje, (nombre, mu) in zip(ejes, casos):
        energias, vectores, hamiltoniano = diagonalizar_cadena(numero_sitios,hopping,mu,pairing_p_wave,) # Diagonalizamos la cadena para este caso
        indice = indice_estado_mas_cercano_a_cero(energias)
        estado = vectores[:, indice] # Extraemos el estado más cercano a cero
        densidad = densidad_estado(estado, numero_sitios) # Calculamos su densidad espacial
        eje.plot(np.arange(numero_sitios), densidad, linewidth=2.3)
        eje.set_title(f"{nombre}: $\\mu={mu}$, $E={energias[indice]:.2e}$")
        eje.set_xlabel("Sitio j")
        eje.set_ylabel("Densidad normalizada")
    figura.suptitle("Estado BdG más cercano a energía cero")
    figura.tight_layout()

    return figura



def graficar_majoranas_separadas(numero_sitios, hopping, potencial_quimico, pairing_p_wave): #Grafica las dos combinaciones tipo Majorana localizadas a izquierda y derecha
    energias, vectores, hamiltoniano = diagonalizar_cadena(numero_sitios,hopping,potencial_quimico,pairing_p_wave,) # Diagonalizamos la cadena para obtener energías y vectores propios
    densidad_izquierda, densidad_derecha = construir_majoranas_desde_pareja(energias,vectores, numero_sitios,) # Construimos las dos combinaciones tipo Majorana a partir de la pareja partícula-hueco cercana a cero

    figura, eje = plt.subplots()
    eje.plot(np.arange(numero_sitios), densidad_izquierda, linewidth=2.3, label="Majorana izquierda")
    eje.plot(np.arange(numero_sitios), densidad_derecha, linewidth=2.3, label="Majorana derecha")
    eje.set_title(r"Descomposición del modo casi cero en dos Majoranas")
    eje.set_xlabel("Sitio j")
    eje.set_ylabel("Densidad normalizada")
    eje.legend()
    figura.tight_layout()

    return figura



def graficar_particula_hueco_modo_cero(numero_sitios, hopping, potencial_quimico, pairing_p_wave): #Grafica la parte electrónica y de hueco del estado más cercano a cero para mostrar su carácter de Majorana
    energias, vectores, hamiltoniano = diagonalizar_cadena( numero_sitios,hopping, potencial_quimico, pairing_p_wave, ) # Diagonalizamos la cadena para obtener energías y vectores propios
    indice = indice_estado_mas_cercano_a_cero(energias)                                  # Encontramos el índice del estado más cercano a cero
    estado = vectores[:, indice]                                                         # Extraemos el estado más cercano a cero
    densidad_electron, densidad_hueco = densidad_electron_y_hueco(estado, numero_sitios) # Calculamos por separado la densidad de la parte electrónica y de la parte de hueco

    figura, eje = plt.subplots()
    eje.plot(np.arange(numero_sitios), densidad_electron, linewidth=2.2, label="parte electronica")
    eje.plot(np.arange(numero_sitios), densidad_hueco, linewidth=2.2, label="parte de hueco")
    eje.set_title(r"Contenido electron-hueco del estado mas cercano a cero")
    eje.set_xlabel("Sitio j")
    eje.set_ylabel("Densidad normalizada")
    eje.legend()
    figura.tight_layout()

    return figura



def graficar_splitting_vs_longitud(hopping, potencial_quimico, pairing_p_wave): #Grafica la energía mínima frente a la longitud N
    valores_longitud = np.arange(12, 121, 2) 
    energias_minimas = barrer_longitud_cadena(valores_longitud,hopping,potencial_quimico,pairing_p_wave,) # Barrido de longitud para obtener la energía mínima en cada caso
    figura, eje = plt.subplots()
    eje.plot(valores_longitud, energias_minimas, marker="o", markersize=3, linewidth=1.6)
    eje.set_yscale("log")
    eje.set_title("Splitting de los modos de borde frente a la longitud")
    eje.set_xlabel("Número de sitios N")
    eje.set_ylabel(r"$\min |E|/t$")
    figura.tight_layout()

    return figura



def graficar_comparacion_pairing_localizacion(numero_sitios, hopping, potencial_quimico): #Compara cómo cambia la localización de los Majoranas al variar Delta
    valores_delta = [pairing_debil, pairing_medio, pairing_fuerte]

    figura, eje = plt.subplots()

    for Delta in valores_delta:
        energias, vectores, hamiltoniano = diagonalizar_cadena(numero_sitios, hopping, potencial_quimico, Delta, )
        indice = indice_estado_mas_cercano_a_cero(energias)
        densidad = densidad_estado(vectores[:, indice], numero_sitios)
        eje.plot( np.arange(numero_sitios), densidad, linewidth=2.0, label=rf"$\Delta={Delta}$",)

    eje.set_title("Efecto de $\\Delta$ en la localización del modo de borde")
    eje.set_xlabel("Sitio j")
    eje.set_ylabel("Densidad normalizada")
    eje.legend()

    figura.tight_layout()

    return figura



def graficar_diagrama_fase_mu_delta(hopping): #Grafica el diagrama de fase ideal en el plano
    valores_mu = np.linspace(-4.0, 4.0, 250)
    valores_delta = np.linspace(0.0, 1.5, 180)

    diagrama = calcular_diagrama_fase_mu_delta(hopping,valores_mu,valores_delta,)

    figura, eje = plt.subplots()

    imagen = eje.imshow(diagrama,origin="lower",aspect="auto",extent=[valores_mu[0], valores_mu[-1], valores_delta[0], valores_delta[-1]],cmap="viridis",)

    eje.axvline(-2.0 * hopping, color="white", linestyle="--", linewidth=2.0)
    eje.axvline(+2.0 * hopping, color="white", linestyle="--", linewidth=2.0)

    eje.set_title(r"Diagrama de fase ideal de la cadena de Kitaev")
    eje.set_xlabel(r"Potencial químico $\mu/t$")
    eje.set_ylabel(r"Pairing $p$-wave $\Delta/t$")

    barra = figura.colorbar(imagen, ax=eje)
    barra.set_label("0 = trivial, 1 = topológico")

    figura.tight_layout()

    return figura



###############################################################################
################################## Graficas ###################################
###############################################################################

if __name__ == "__main__":
    preparar_estilo_graficas()

    graficar_bandas_bulk_tres_regimenes(hopping, pairing_p_wave)
    graficar_espectro_finito_vs_mu(numero_sitios, hopping, pairing_p_wave)
    graficar_energia_minima_vs_mu(numero_sitios, hopping, pairing_p_wave)
    graficar_ldos_tres_regimenes(numero_sitios, hopping, pairing_p_wave)
    graficar_modo_cero_topologico_y_trivial(numero_sitios, hopping, pairing_p_wave)
    graficar_majoranas_separadas(numero_sitios, hopping, mu_topologico, pairing_p_wave)
    graficar_particula_hueco_modo_cero(numero_sitios, hopping, mu_topologico, pairing_p_wave)
    graficar_splitting_vs_longitud(hopping, mu_topologico, pairing_p_wave)
    graficar_comparacion_pairing_localizacion(numero_sitios, hopping, mu_topologico)
    graficar_diagrama_fase_mu_delta(hopping)

    plt.show()

