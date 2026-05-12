"""
self_consistent_bdg_nanowire.py

Numerical study of Majorana zero modes in one-dimensional superconducting systems.
This file was prepared for a GitHub repository from the original practice scripts.
The numerical model is intentionally explicit: dense BdG matrices are built in real space
so that the Hamiltonian, observables, and generated figures can be read directly.

For publication-quality figures, increase the site counts carefully. Dense diagonalization
scales poorly with system size.
"""

import os
import numpy as np
import scipy.linalg as la
import matplotlib.pyplot as plt

###############################################################################
################ Constantes fisicas y parametros del nanohilo #################
###############################################################################

# Trabajamos en unidades practicas: energia en meV, longitud en nm y temperatura en K.
meV = 1.0
nm = 1.0
kelvin = 1.0

constante_boltzmann = 0.08617333262 * meV / kelvin              # k_B en meV/K.
hbar2_sobre_2m_e = 38.0998212 * meV * nm**2                     # hbar^2/(2m_e) en meV nm^2.

masa_efectiva = 0.015                                           # m*/m_e.
alpha_rashba = 100.0 * meV * nm                                 # Acoplo Rashba alpha_R.
numero_sitios = 180                                             # Valor para probar. Para figuras finales prueba 180, 220 o 260 si el Mac aguanta.
numero_sitios_mapa = 80                                         # Valor para mapas 2D. Para final puedes probar 48 o 56.
paso_espacial = 20.0 * nm
hopping = hbar2_sobre_2m_e / (masa_efectiva * paso_espacial**2) # t = hbar^2/(2m*a^2).
hopping_spin_orbita = alpha_rashba / (2.0 * paso_espacial)      # t_so = alpha_R/(2a).
potencial_quimico = 0.0 * meV                                   # mu.
temperatura = 0.02 * kelvin                                     # Temperatura para tanh(E/2kBT).
interaccion_atractiva = 4.0 * meV                               # V_eff positivo: mide la atraccion efectiva.
energia_corte = 4.0 * meV                                       # Corte energetico para la suma autoconsistente.
delta_semilla = 0.30 * meV                                      # Semilla inicial del gap; el valor final sale de la autoconsistencia.
delta_minima_superconductora = 10**(-3) * meV                   # Por debajo de esto consideramos que el gap colapso.

# Parametros numericos del bucle de autoconsistencia.
max_iteraciones = 120
tolerancia_delta = 2*10**(-5) * meV
mezcla_delta = 0.35

# Barridos y visualizacion.
numero_energias = 360
ensanchamiento_ldos = 0.035 * meV
numero_puntos_barrido = 31
numero_k_bulk = 450
carpeta_figuras = "figuras_autoconsistente"

# Casos representativos. El valor intermedio esta cerca de la zona donde el gap tiende a cerrarse.
zeeman_trivial = 0.0 * meV
zeeman_intermedio = 0.50 * meV
zeeman_topologico = 0.56 * meV

# Valores de V_eff para comparar como cambia el gap autoconsistente y la localizacion del modo.
interaccion_debil = 4.0 * meV
interaccion_media = 4.1 * meV
interaccion_fuerte = 4.2 * meV


def calcular_hoppings_discretos(alpha_rashba_usado, paso_espacial_usado): #Calcula t y t_so para un valor concreto de alpha_R y del paso a.
    hopping_calculado = hbar2_sobre_2m_e / (masa_efectiva * paso_espacial_usado**2)
    hopping_so_calculado = alpha_rashba_usado / (2.0 * paso_espacial_usado)
    return hopping_calculado, hopping_so_calculado


###################### MATRICES DE PAULI ######################

matriz_identidad_2 = np.eye(2, dtype=complex)
matriz_sigma_x_2 = np.array([[0.0, 1.0], [1.0, 0.0]], dtype=complex)
matriz_sigma_y_2 = np.array([[0.0, -1.0j], [1.0j, 0.0]], dtype=complex)
matriz_sigma_z_2 = np.array([[1.0, 0.0], [0.0, -1.0]], dtype=complex)

matriz_tau_x_2 = np.array([[0.0, 1.0], [1.0, 0.0]], dtype=complex)
matriz_tau_y_2 = np.array([[0.0, -1.0j], [1.0j, 0.0]], dtype=complex)
matriz_tau_z_2 = np.array([[1.0, 0.0], [0.0, -1.0]], dtype=complex)

# En la base (Nambu, espin) usamos productos de Kronecker tau \otimes sigma.
sigma_x = np.kron(matriz_identidad_2, matriz_sigma_x_2)
sigma_y = np.kron(matriz_identidad_2, matriz_sigma_y_2)
tau_x = np.kron(matriz_tau_x_2, matriz_identidad_2)
tau_y = np.kron(matriz_tau_y_2, matriz_identidad_2)
tau_z = np.kron(matriz_tau_z_2, matriz_identidad_2)
sigma_y_tau_z = sigma_y @ tau_z

# Operador unitario de la simetria particula-hueco en esta base: C psi = (tau_y sigma_y) psi^*.
operador_particula_hueco_local = tau_y @ sigma_y


###################### PERFILES ESPACIALES ######################

def convertir_a_perfil_delta(delta, numero_sitios): #Convierte un numero complejo o un array en un perfil Delta_i de longitud N.
    if np.ndim(delta) == 0:
        return np.full(numero_sitios, complex(delta), dtype=complex)
    return np.asarray(delta, dtype=complex).copy()


def perfil_delta_inicial_uniforme(numero_sitios, delta_inicial): #Semilla uniforme para iniciar el bucle autoconsistente.
    return np.full(numero_sitios, delta_inicial + 0.0j, dtype=complex)


def potencial_electrostatico_nulo(numero_sitios): #Devuelve V_i = 0 en todo el hilo.
    return np.zeros(numero_sitios, dtype=float)


###############################################################################
############################## Hamiltoniano BdG ###############################
###############################################################################

def construir_hamiltoniano_bdg(numero_sitios, potencial_quimico, energia_zeeman, perfil_delta, potencial, hopping_usado, hopping_spin_orbita_usado): #Construye la matriz BdG 4N x 4N del nanohilo.
    potencial = np.asarray(potencial, dtype=float).copy()
    perfil_delta = convertir_a_perfil_delta(perfil_delta, numero_sitios)
    if len(potencial) != numero_sitios:
        raise ValueError("El perfil V_i no tiene el mismo numero de sitios que la cadena.")
    if len(perfil_delta) != numero_sitios:
        raise ValueError("El perfil Delta_i no tiene el mismo numero de sitios que la cadena.")

    hamiltoniano = np.zeros((4 * numero_sitios, 4 * numero_sitios), dtype=complex)
    bloque_hopping = -hopping_usado * tau_z + 1.0j * hopping_spin_orbita_usado * sigma_y_tau_z

    for i in range(numero_sitios):
        bloque_i = slice(4 * i, 4 * (i + 1))
        delta_i = perfil_delta[i]
        bloque_local = (2.0 * hopping_usado - potencial_quimico + potencial[i]) * tau_z + energia_zeeman * sigma_x + delta_i.real * tau_x - delta_i.imag * tau_y
        hamiltoniano[bloque_i, bloque_i] = bloque_local
        if i < numero_sitios - 1:
            bloque_j = slice(4 * (i + 1), 4 * (i + 2))
            hamiltoniano[bloque_i, bloque_j] = bloque_hopping
            hamiltoniano[bloque_j, bloque_i] = bloque_hopping.conj().T

    hamiltoniano = 0.5 * (hamiltoniano + hamiltoniano.conj().T)
    return hamiltoniano


def diagonalizar_hamiltoniano(hamiltoniano): #Diagonaliza una matriz BdG hermitica y ordena los autovalores.
    energias, vectores = la.eigh(hamiltoniano, overwrite_a=True, check_finite=False, driver="evd")
    orden = np.argsort(energias)
    energias = energias[orden]
    vectores = vectores[:, orden]
    return energias, vectores


def error_hermiticidad(hamiltoniano): #Mide cuanto se aparta H de su conjugada hermitica.
    norma_h = max(la.norm(hamiltoniano), 10**(-30))
    return float(la.norm(hamiltoniano - hamiltoniano.conj().T) / norma_h)


def error_simetria_particula_hueco(hamiltoniano, numero_sitios): #Comprueba C H^* C^dagger = -H.
    operador_global = np.kron(np.eye(numero_sitios, dtype=complex), operador_particula_hueco_local)
    prueba = operador_global @ hamiltoniano.conj() @ operador_global.conj().T + hamiltoniano
    norma_h = max(la.norm(hamiltoniano), 10**(-30))
    return float(la.norm(prueba) / norma_h)


###############################################################################
################### Autoconsistencia del gap superconductivo ##################
###############################################################################

def fijar_fase_global_delta(perfil_delta): #Fija la fase global de Delta para que el promedio sea real y positivo.
    perfil_delta = np.asarray(perfil_delta, dtype=complex).copy()
    promedio = np.mean(perfil_delta)
    if abs(promedio) > 10**(-14):
        perfil_delta = perfil_delta * np.exp(-1.0j * np.angle(promedio))
    if np.mean(perfil_delta.real) < 0.0:
        perfil_delta = -perfil_delta
    return perfil_delta


def calcular_delta_autoconsistente(energias, vectores, numero_sitios, interaccion_atractiva, temperatura_usada, energia_corte_usada): #Calcula Delta_i desde las soluciones BdG.
    delta_nueva = np.zeros(numero_sitios, dtype=complex)
    for n in range(len(energias)):
        energia_n = energias[n]
        if abs(energia_n) > energia_corte_usada:
            continue
        psi = vectores[:, n].reshape(numero_sitios, 4)
        u_arriba = psi[:, 0]
        u_abajo = psi[:, 1]
        v_abajo = psi[:, 2]
        menos_v_arriba = psi[:, 3]
        correlacion_anomala = u_arriba * np.conj(v_abajo) + u_abajo * np.conj(menos_v_arriba)
        factor_termico = np.tanh(energia_n / (2.0 * constante_boltzmann * temperatura_usada))
        delta_nueva += 0.5 * interaccion_atractiva * correlacion_anomala * factor_termico
    delta_nueva = fijar_fase_global_delta(delta_nueva)
    return delta_nueva


def diagonalizar_autoconsistente(numero_sitios, potencial_quimico, energia_zeeman, interaccion_atractiva, delta_inicial, potencial, hopping_usado, hopping_spin_orbita_usado, temperatura_usada, energia_corte_usada, max_iteraciones_usadas, tolerancia_usada, mezcla_usada): #Resuelve el problema BdG autoconsistente.
    perfil_delta = fijar_fase_global_delta(convertir_a_perfil_delta(delta_inicial, numero_sitios))
    historial_error = []
    historial_delta_media = []
    convergido = False

    for iteracion in range(max_iteraciones_usadas):
        hamiltoniano = construir_hamiltoniano_bdg(numero_sitios, potencial_quimico, energia_zeeman, perfil_delta, potencial, hopping_usado, hopping_spin_orbita_usado)
        energias, vectores = diagonalizar_hamiltoniano(hamiltoniano)
        delta_calculada = calcular_delta_autoconsistente(energias, vectores, numero_sitios, interaccion_atractiva, temperatura_usada, energia_corte_usada)
        delta_siguiente = (1.0 - mezcla_usada) * perfil_delta + mezcla_usada * delta_calculada
        delta_siguiente = fijar_fase_global_delta(delta_siguiente)
        error = np.max(np.abs(delta_siguiente - perfil_delta))
        historial_error.append(float(error))
        historial_delta_media.append(float(np.mean(np.abs(delta_siguiente))))
        perfil_delta = delta_siguiente
        if error < tolerancia_usada:
            convergido = True
            break

    hamiltoniano = construir_hamiltoniano_bdg(numero_sitios, potencial_quimico, energia_zeeman, perfil_delta, potencial, hopping_usado, hopping_spin_orbita_usado)
    energias, vectores = diagonalizar_hamiltoniano(hamiltoniano)
    resultado = {"energias": energias, "vectores": vectores, "hamiltoniano": hamiltoniano, "delta": perfil_delta, "delta_media": float(np.mean(np.abs(perfil_delta))), "historial_error": np.array(historial_error), "historial_delta_media": np.array(historial_delta_media), "convergido": convergido, "iteraciones": iteracion + 1, "error_hermiticidad": error_hermiticidad(hamiltoniano), "error_particula_hueco": error_simetria_particula_hueco(hamiltoniano, numero_sitios)}
    return resultado


def resolver_autoconsistente_estandar(numero_sitios, potencial_quimico, energia_zeeman, interaccion_atractiva, delta_inicial): #Llama al solver usando siempre todos los parametros fisicos y numericos.
    potencial = potencial_electrostatico_nulo(numero_sitios)
    return diagonalizar_autoconsistente(numero_sitios, potencial_quimico, energia_zeeman, interaccion_atractiva, delta_inicial, potencial, hopping, hopping_spin_orbita, temperatura, energia_corte, max_iteraciones, tolerancia_delta, mezcla_delta)


###############################################################################
############################# Observables fisicos #############################
###############################################################################

def separar_componentes_bdg(vector_estado, numero_sitios): #Separa un autovector en u_up, u_down, v_down y -v_up.
    psi = vector_estado.reshape(numero_sitios, 4)
    u_arriba = psi[:, 0].copy()
    u_abajo = psi[:, 1].copy()
    v_abajo = psi[:, 2].copy()
    menos_v_arriba = psi[:, 3].copy()
    return u_arriba, u_abajo, v_abajo, menos_v_arriba


def densidad_estado(vector_estado, numero_sitios): #Calcula |u_up|^2 + |u_down|^2 + |v_up|^2 + |v_down|^2.
    u_arriba, u_abajo, v_abajo, menos_v_arriba = separar_componentes_bdg(vector_estado, numero_sitios)
    densidad = np.abs(u_arriba)**2 + np.abs(u_abajo)**2 + np.abs(v_abajo)**2 + np.abs(menos_v_arriba)**2
    densidad = densidad / max(float(np.sum(densidad)), 10**(-30))
    return densidad


def densidad_electron_y_hueco(vector_estado, numero_sitios): #Devuelve por separado la densidad electronica y la densidad de hueco.
    u_arriba, u_abajo, v_abajo, menos_v_arriba = separar_componentes_bdg(vector_estado, numero_sitios)
    densidad_electron = np.abs(u_arriba)**2 + np.abs(u_abajo)**2
    densidad_hueco = np.abs(v_abajo)**2 + np.abs(menos_v_arriba)**2
    norma = max(float(np.sum(densidad_electron + densidad_hueco)), 10**(-30))
    return densidad_electron / norma, densidad_hueco / norma


def indice_estado_mas_cercano_a_cero(energias): #Encuentra el autovalor mas cercano a E = 0.
    return int(np.argmin(np.abs(energias)))


def energia_minima_absoluta(energias): #Devuelve min |E|.
    return float(np.min(np.abs(energias)))


def transformar_particula_hueco(vector_estado, numero_sitios): #Aplica C = tau_y sigma_y K.
    psi = vector_estado.reshape(numero_sitios, 4)
    psi_ph = np.zeros_like(psi)
    for i in range(numero_sitios):
        psi_ph[i, :] = operador_particula_hueco_local @ np.conj(psi[i, :])
    return psi_ph.reshape(4 * numero_sitios)


def alinear_fase(estado_referencia, estado_objetivo): #Elimina la fase global arbitraria entre dos autovectores.
    solapamiento = np.vdot(estado_referencia, estado_objetivo)
    if abs(solapamiento) < 10**(-30):
        return estado_objetivo.copy()
    return estado_objetivo * np.exp(-1.0j * np.angle(solapamiento))


def construir_majoranas_desde_pareja(energias, vectores, numero_sitios): #Construye dos combinaciones tipo Majorana desde la pareja +-E mas cercana a cero.
    indices_positivos = np.where(energias >= 0.0)[0]
    indices_negativos = np.where(energias < 0.0)[0]
    indice_positivo = indices_positivos[np.argmin(np.abs(energias[indices_positivos]))]
    indice_negativo = indices_negativos[np.argmin(np.abs(energias[indices_negativos] + energias[indice_positivo]))]
    estado_positivo = vectores[:, indice_positivo]
    estado_negativo = vectores[:, indice_negativo]
    estado_phs = transformar_particula_hueco(estado_positivo, numero_sitios)
    estado_negativo = alinear_fase(estado_phs, estado_negativo)
    majorana_1 = (estado_positivo + estado_negativo) / np.sqrt(2.0)
    majorana_2 = -1.0j * (estado_positivo - estado_negativo) / np.sqrt(2.0)
    densidad_1 = densidad_estado(majorana_1, numero_sitios)
    densidad_2 = densidad_estado(majorana_2, numero_sitios)
    mitad = numero_sitios // 2
    if np.sum(densidad_1[:mitad]) >= np.sum(densidad_2[:mitad]):
        densidad_izquierda = densidad_1
        densidad_derecha = densidad_2
    else:
        densidad_izquierda = densidad_2
        densidad_derecha = densidad_1
    return densidad_izquierda, densidad_derecha


def lorentziana(energia, centro, anchura): #Funcion lorentziana para suavizar la LDOS.
    return anchura / (np.pi * ((energia - centro)**2 + anchura**2))


def calcular_ldos(numero_sitios, energias, vectores, energias_ldos, anchura): #Calcula la LDOS electronica local.
    ldos = np.zeros((len(energias_ldos), numero_sitios), dtype=float)
    for n, energia_n in enumerate(energias):
        if energia_n <= 1*10**(-12):
            continue
        psi = vectores[:, n].reshape(numero_sitios, 4)
        peso_electron = np.abs(psi[:, 0])**2 + np.abs(psi[:, 1])**2
        peso_hueco = np.abs(psi[:, 2])**2 + np.abs(psi[:, 3])**2
        for i, energia in enumerate(energias_ldos):
            ldos[i, :] += peso_electron * lorentziana(energia, energia_n, anchura) + peso_hueco * lorentziana(energia, -energia_n, anchura)
    return ldos


###############################################################################
######################### Bandas bulk del sistema infinito ####################
###############################################################################

def bandas_bulk_por_diagonalizacion(potencial_quimico, energia_zeeman, delta_uniforme, numero_k, hopping_usado, hopping_spin_orbita_usado): #Bandas bulk diagonalizando H(k) punto a punto.
    k = np.linspace(-np.pi, np.pi, numero_k)
    bandas = np.zeros((numero_k, 4), dtype=float)
    for i, k_i in enumerate(k):
        xi_k = 2.0 * hopping_usado - potencial_quimico - 2.0 * hopping_usado * np.cos(k_i)
        rashba_k = -2.0 * hopping_spin_orbita_usado * np.sin(k_i)
        hamiltoniano_k = xi_k * tau_z + rashba_k * sigma_y_tau_z + energia_zeeman * sigma_x + np.real(delta_uniforme) * tau_x - np.imag(delta_uniforme) * tau_y
        bandas[i, :] = la.eigvalsh(hamiltoniano_k, check_finite=False)
    return k, bandas


def bandas_bulk_analiticas(potencial_quimico, energia_zeeman, delta_uniforme, numero_k, hopping_usado, hopping_spin_orbita_usado): #Solucion analitica discreta para comparar con la diagonalizacion de H(k).
    k = np.linspace(-np.pi, np.pi, numero_k)
    xi_k = 2.0 * hopping_usado - potencial_quimico - 2.0 * hopping_usado * np.cos(k)
    rashba_k = -2.0 * hopping_spin_orbita_usado * np.sin(k)
    delta_abs = np.abs(delta_uniforme)
    termino_base = xi_k**2 + rashba_k**2 + energia_zeeman**2 + delta_abs**2
    termino_raiz = np.sqrt(energia_zeeman**2 * delta_abs**2 + xi_k**2 * (energia_zeeman**2 + rashba_k**2))
    energia_menor = np.sqrt(np.maximum(termino_base - 2.0 * termino_raiz, 0.0))
    energia_mayor = np.sqrt(np.maximum(termino_base + 2.0 * termino_raiz, 0.0))
    bandas = np.column_stack((-energia_mayor, -energia_menor, energia_menor, energia_mayor))
    return k, bandas


def gap_bulk_minimo(potencial_quimico, energia_zeeman, delta_uniforme, numero_k, hopping_usado, hopping_spin_orbita_usado): #Calcula min_k E_+(k).
    k, bandas = bandas_bulk_por_diagonalizacion(potencial_quimico, energia_zeeman, delta_uniforme, numero_k, hopping_usado, hopping_spin_orbita_usado)
    bandas_positivas = bandas[bandas > 10**(-12)]
    if bandas_positivas.size == 0:
        return 0.0
    return float(np.min(bandas_positivas))


def es_topologico_ideal(potencial_quimico, energia_zeeman, delta_media): #Criterio orientativo E_Z > sqrt(mu^2 + |Delta|^2).
    umbral = np.sqrt(potencial_quimico**2 + delta_media**2)
    return bool(energia_zeeman > umbral and delta_media > delta_minima_superconductora)


###############################################################################
############################## Barridos fisicos ###############################
###############################################################################

def barrer_zeeman(numero_sitios, potencial_quimico, valores_zeeman, interaccion_atractiva): #Barre E_Z y calcula el espectro finito autoconsistente.
    todos_los_espectros = []
    energias_minimas = []
    gaps_medios = []
    deltas_finales = []
    delta_inicial = perfil_delta_inicial_uniforme(numero_sitios, delta_semilla)
    for energia_zeeman in valores_zeeman:
        resultado = resolver_autoconsistente_estandar(numero_sitios, potencial_quimico, energia_zeeman, interaccion_atractiva, delta_inicial)
        todos_los_espectros.append(resultado["energias"])
        energias_minimas.append(energia_minima_absoluta(resultado["energias"]))
        gaps_medios.append(resultado["delta_media"])
        deltas_finales.append(resultado["delta"])
        delta_inicial = resultado["delta"]
    return np.array(todos_los_espectros), np.array(energias_minimas), np.array(gaps_medios), deltas_finales


def barrer_gap_bulk_zeeman(valores_zeeman, numero_sitios_barrido): #Barre E_Z y calcula el gap bulk usando Delta autoconsistente del centro del hilo.
    gaps_bulk = []
    gaps_medios = []
    delta_inicial = perfil_delta_inicial_uniforme(numero_sitios_barrido, delta_semilla)
    for energia_zeeman in valores_zeeman:
        resultado = resolver_autoconsistente_estandar(numero_sitios_barrido, potencial_quimico, energia_zeeman, interaccion_atractiva, delta_inicial)
        zona_central = resultado["delta"][numero_sitios_barrido // 4 : 3 * numero_sitios_barrido // 4]
        delta_bulk = np.mean(zona_central)
        gaps_bulk.append(gap_bulk_minimo(potencial_quimico, energia_zeeman, delta_bulk, numero_k_bulk, hopping, hopping_spin_orbita))
        gaps_medios.append(resultado["delta_media"])
        delta_inicial = resultado["delta"]
    return np.array(gaps_bulk), np.array(gaps_medios)


def barrer_longitud_cadena(valores_longitud, potencial_quimico, energia_zeeman, interaccion_atractiva): #Barre la longitud N para estudiar el splitting de los modos de borde.
    energias_minimas = []
    for N in valores_longitud:
        N = int(N)
        delta_inicial = perfil_delta_inicial_uniforme(N, delta_semilla)
        resultado = resolver_autoconsistente_estandar(N, potencial_quimico, energia_zeeman, interaccion_atractiva, delta_inicial)
        energias_minimas.append(energia_minima_absoluta(resultado["energias"]))
    return np.array(energias_minimas)


def calcular_mapa_mu_zeeman_minimo(valores_mu, valores_zeeman, numero_sitios_mapa): #Calcula log10(min|E|) en el plano (mu,E_Z) para la cadena finita.
    mapa_minimo = np.zeros((len(valores_zeeman), len(valores_mu)), dtype=float)
    for i, energia_zeeman in enumerate(valores_zeeman):
        delta_inicial = perfil_delta_inicial_uniforme(numero_sitios_mapa, delta_semilla)
        for j, mu in enumerate(valores_mu):
            resultado = resolver_autoconsistente_estandar(numero_sitios_mapa, mu, energia_zeeman, interaccion_atractiva, delta_inicial)
            mapa_minimo[i, j] = energia_minima_absoluta(resultado["energias"])
            delta_inicial = resultado["delta"]
    mapa_log = np.log10(np.maximum(mapa_minimo, 10**(-6) * meV))
    return mapa_log


def calcular_resultados_regimenes(numero_sitios): #Calcula una vez los tres regimenes principales para no repetir diagonalizaciones caras.
    resultados = {}
    casos = [("trivial", zeeman_trivial), ("intermedio", zeeman_intermedio), ("topologico", zeeman_topologico)]
    for nombre, energia_zeeman in casos:
        delta_inicial = perfil_delta_inicial_uniforme(numero_sitios, delta_semilla)
        resultados[nombre] = resolver_autoconsistente_estandar(numero_sitios, potencial_quimico, energia_zeeman, interaccion_atractiva, delta_inicial)
    return resultados


def calcular_resultados_interaccion(numero_sitios): #Calcula los perfiles al variar V_eff.
    resultados = []
    for valor_interaccion in [interaccion_debil, interaccion_media, interaccion_fuerte]:
        delta_inicial = perfil_delta_inicial_uniforme(numero_sitios, delta_semilla)
        resultado = resolver_autoconsistente_estandar(numero_sitios, potencial_quimico, zeeman_topologico, valor_interaccion, delta_inicial)
        resultados.append((valor_interaccion, resultado))
    return resultados


###############################################################################
################################## Graficar ###################################
###############################################################################

def preparar_estilo_graficas(): #Estilo sencillo para que las figuras sean claras en el TFG.
    plt.style.use("default")
    plt.rcParams.update({"figure.figsize": (9, 5.5), "axes.grid": True, "grid.alpha": 0.25, "font.size": 11, "axes.titlesize": 13, "axes.labelsize": 12, "legend.fontsize": 10, "savefig.dpi": 300, "savefig.bbox": "tight"})


def guardar_figura(figura, nombre_archivo): #Guarda la figura en la carpeta de resultados.
    os.makedirs(carpeta_figuras, exist_ok=True)
    figura.savefig(os.path.join(carpeta_figuras, nombre_archivo))


def graficar_convergencia_gap(resultado_topologico): #Muestra como converge Delta durante el proceso autoconsistente.
    figura, ejes = plt.subplots(1, 2, figsize=(12, 4))
    ejes[0].plot(resultado_topologico["historial_delta_media"], linewidth=2.2)
    ejes[0].set_title(r"Convergencia de $\langle|\Delta_i|\rangle$")
    ejes[0].set_xlabel("Iteracion")
    ejes[0].set_ylabel("Gap medio (meV)")
    ejes[1].semilogy(resultado_topologico["historial_error"], linewidth=2.2)
    ejes[1].set_title("Error maximo entre iteraciones")
    ejes[1].set_xlabel("Iteracion")
    ejes[1].set_ylabel("Error (meV)")
    figura.suptitle("Bucle autoconsistente BdG")
    figura.tight_layout()
    guardar_figura(figura, "convergencia_gap_autoconsistente.png")
    return figura


def graficar_perfil_gap_autoconsistente(resultados_regimenes): #Dibuja Delta_i en los tres regimenes.
    figura, eje = plt.subplots(figsize=(9, 5))
    casos = [("trivial", zeeman_trivial), ("intermedio", zeeman_intermedio), ("topologico", zeeman_topologico)]
    for nombre, energia_zeeman in casos:
        resultado = resultados_regimenes[nombre]
        eje.plot(np.arange(numero_sitios), np.abs(resultado["delta"]), linewidth=2.0, label=f"{nombre}: $E_Z={energia_zeeman:.2f}$ meV")
    eje.set_title(r"Perfil espacial del gap autoconsistente $|\Delta_i|$")
    eje.set_xlabel("Sitio i")
    eje.set_ylabel(r"$|\Delta_i|$ (meV)")
    eje.legend()
    figura.tight_layout()
    guardar_figura(figura, "perfil_gap_autoconsistente.png")
    return figura


def graficar_gap_bulk_vs_zeeman(valores_zeeman, gaps_bulk): #Grafica el gap bulk minimo frente a E_Z.
    figura, eje = plt.subplots(figsize=(9, 5))
    eje.plot(valores_zeeman, gaps_bulk, marker="o", markersize=4, linewidth=2.0)
    eje.set_title(r"Gap bulk minimo frente a $E_Z$")
    eje.set_xlabel(r"Energia de Zeeman $E_Z$ (meV)")
    eje.set_ylabel(r"$E_{\mathrm{gap}}^{\mathrm{bulk}}$ (meV)")
    figura.tight_layout()
    guardar_figura(figura, "gap_bulk_vs_zeeman_autoconsistente.png")
    return figura


def graficar_bandas_bulk_tres_regimenes(resultados_regimenes): #Grafica bandas teoricas y bandas obtenidas al diagonalizar H(k).
    casos = [("trivial", zeeman_trivial), ("intermedio", zeeman_intermedio), ("topologico", zeeman_topologico)]
    figura, ejes = plt.subplots(1, 3, figsize=(14, 4.2), sharey=True)
    colores = ["tab:blue", "tab:orange", "tab:green", "tab:red"]
    for eje, (nombre, energia_zeeman) in zip(ejes, casos):
        resultado = resultados_regimenes[nombre]
        zona_central = resultado["delta"][numero_sitios // 4 : 3 * numero_sitios // 4]
        delta_bulk = np.mean(zona_central)
        k, bandas_teoricas = bandas_bulk_analiticas(potencial_quimico, energia_zeeman, delta_bulk, numero_k_bulk, hopping, hopping_spin_orbita)
        k, bandas_numericas = bandas_bulk_por_diagonalizacion(potencial_quimico, energia_zeeman, delta_bulk, numero_k_bulk, hopping, hopping_spin_orbita)
        for banda in range(4):
            etiqueta_teorica = "formula analitica" if banda == 0 and nombre == "trivial" else None
            etiqueta_numerica = "diagonalizacion de $H(k)$" if banda == 0 and nombre == "trivial" else None
            eje.plot(k / np.pi, bandas_teoricas[:, banda], linewidth=1.8, color=colores[banda], label=etiqueta_teorica)
            eje.plot(k / np.pi, bandas_numericas[:, banda], linestyle="--", linewidth=1.2, color="black", alpha=0.75, label=etiqueta_numerica)
        eje.axhline(0.0, color="black", linewidth=1.0, alpha=0.45)
        eje.set_title(rf"{nombre}: $E_Z={energia_zeeman:.2f}$ meV, $\langle\Delta\rangle={resultado['delta_media']:.2f}$ meV")
        eje.set_xlabel(r"$ka/\pi$")
    ejes[0].set_ylabel("Energia (meV)")
    ejes[0].legend(loc="lower left", framealpha=0.9)
    figura.suptitle("Bandas bulk del nanohilo BdG autoconsistente")
    figura.tight_layout()
    guardar_figura(figura, "bandas_bulk_autoconsistente.png")
    return figura


def graficar_espectro_finito_vs_zeeman(valores_zeeman, espectros, gaps_medios): #Grafica el espectro finito frente a E_Z.
    figura, eje = plt.subplots(figsize=(9, 5.2))
    for i, energia_zeeman in enumerate(valores_zeeman):
        energias_cerca_cero = espectros[i][np.abs(espectros[i]) < 1.5 * meV]
        eje.scatter(np.full_like(energias_cerca_cero, energia_zeeman), energias_cerca_cero, s=8, color="black", alpha=0.65)
    delta_referencia = max(gaps_medios[0], 10**(-12))
    zeeman_guia = np.sqrt(potencial_quimico**2 + delta_referencia**2)
    eje.axvline(zeeman_guia, linestyle="--", linewidth=2.0, label=r"$E_Z\simeq\sqrt{\mu^2+\Delta^2}$")
    eje.axhline(0.0, color="black", linewidth=1.0, alpha=0.5)
    eje.set_title("Espectro finito autoconsistente frente al campo de Zeeman")
    eje.set_xlabel(r"Energia de Zeeman $E_Z$ (meV)")
    eje.set_ylabel("Energia (meV)")
    eje.legend()
    figura.tight_layout()
    guardar_figura(figura, "espectro_finito_vs_zeeman_autoconsistente.png")
    return figura


def graficar_ldos_tres_regimenes(resultados_regimenes): #Grafica LDOS(E,i) con la barra de color fuera de los paneles.
    casos = [("trivial", zeeman_trivial), ("intermedio", zeeman_intermedio), ("topologico", zeeman_topologico)]
    energias_ldos = np.linspace(-1.5 * meV, 1.5 * meV, numero_energias)
    figura = plt.figure(figsize=(14, 4.4))
    rejilla = figura.add_gridspec(1, 4, width_ratios=[1, 1, 1, 0.045], wspace=0.12)
    ejes = [figura.add_subplot(rejilla[0, i]) for i in range(3)]
    eje_barra = figura.add_subplot(rejilla[0, 3])
    imagen = None
    for eje, (nombre, energia_zeeman) in zip(ejes, casos):
        resultado = resultados_regimenes[nombre]
        ldos = calcular_ldos(numero_sitios, resultado["energias"], resultado["vectores"], energias_ldos, ensanchamiento_ldos)
        imagen = eje.imshow(ldos, aspect="auto", origin="lower", extent=[0, numero_sitios - 1, energias_ldos[0], energias_ldos[-1]], cmap="magma")
        eje.axhline(0.0, color="white", linewidth=1.0, alpha=0.85)
        eje.set_title(f"{nombre}: $E_Z={energia_zeeman:.2f}$ meV")
        eje.set_xlabel("Sitio i")
    ejes[0].set_ylabel("Energia (meV)")
    for eje in ejes[1:]:
        eje.tick_params(labelleft=False)
    barra = figura.colorbar(imagen, cax=eje_barra)
    barra.set_label("LDOS")
    figura.suptitle(r"LDOS$(E,i)$ del nanohilo autoconsistente")
    guardar_figura(figura, "ldos_tres_regimenes_autoconsistente.png")
    return figura


def graficar_modo_cero_topologico_y_trivial(resultados_regimenes): #Compara la densidad del estado mas cercano a cero en un caso trivial y uno topologico.
    casos = [("trivial", zeeman_trivial), ("topologico", zeeman_topologico)]
    figura, ejes = plt.subplots(1, 2, figsize=(12, 4.2), sharey=True)
    for eje, (nombre, energia_zeeman) in zip(ejes, casos):
        resultado = resultados_regimenes[nombre]
        indice = indice_estado_mas_cercano_a_cero(resultado["energias"])
        densidad = densidad_estado(resultado["vectores"][:, indice], numero_sitios)
        eje.plot(np.arange(numero_sitios), densidad, linewidth=2.2)
        eje.set_title(f"{nombre}: $E={resultado['energias'][indice]:.3e}$ meV")
        eje.set_xlabel("Sitio i")
    ejes[0].set_ylabel("Densidad normalizada")
    figura.suptitle("Estado BdG mas cercano a energia cero")
    figura.tight_layout()
    guardar_figura(figura, "modo_cero_topologico_y_trivial_autoconsistente.png")
    return figura


def graficar_particula_hueco_modo_cero(resultado_topologico): #Grafica la parte electronica y de hueco del estado mas cercano a cero.
    indice = indice_estado_mas_cercano_a_cero(resultado_topologico["energias"])
    estado = resultado_topologico["vectores"][:, indice]
    densidad_electron, densidad_hueco = densidad_electron_y_hueco(estado, numero_sitios)
    figura, eje = plt.subplots(figsize=(9, 5))
    eje.plot(np.arange(numero_sitios), densidad_electron, linewidth=2.2, label="parte electronica")
    eje.plot(np.arange(numero_sitios), densidad_hueco, linewidth=2.2, linestyle="--", label="parte de hueco")
    eje.set_title("Contenido electron-hueco del estado mas cercano a cero")
    eje.set_xlabel("Sitio i")
    eje.set_ylabel("Densidad normalizada")
    eje.legend()
    figura.tight_layout()
    guardar_figura(figura, "particula_hueco_modo_cero_autoconsistente.png")
    return figura


def graficar_majoranas_separadas(resultado_topologico): #Grafica las dos combinaciones tipo Majorana localizadas a izquierda y derecha.
    densidad_izquierda, densidad_derecha = construir_majoranas_desde_pareja(resultado_topologico["energias"], resultado_topologico["vectores"], numero_sitios)
    figura, eje = plt.subplots(figsize=(9, 5))
    eje.plot(np.arange(numero_sitios), densidad_izquierda, linewidth=2.3, label="Majorana izquierda")
    eje.plot(np.arange(numero_sitios), densidad_derecha, linewidth=2.3, label="Majorana derecha")
    eje.set_title("Descomposicion del modo casi cero en dos Majoranas")
    eje.set_xlabel("Sitio i")
    eje.set_ylabel("Densidad normalizada")
    eje.legend()
    figura.tight_layout()
    guardar_figura(figura, "majoranas_separadas_autoconsistente.png")
    return figura


def graficar_splitting_vs_longitud(): #Grafica la energia minima frente a la longitud N.
    valores_longitud = np.arange(40, 181, 10)
    energias_minimas = barrer_longitud_cadena(valores_longitud, potencial_quimico, zeeman_topologico, interaccion_atractiva)
    figura, eje = plt.subplots(figsize=(9, 5))
    eje.plot(valores_longitud, energias_minimas, marker="o", markersize=4, linewidth=1.7)
    eje.set_yscale("log")
    eje.set_title("Splitting de los modos de borde frente a la longitud")
    eje.set_xlabel("Numero de sitios N")
    eje.set_ylabel(r"$\min |E|$ (meV)")
    figura.tight_layout()
    guardar_figura(figura, "splitting_longitud_autoconsistente.png")
    return figura


def graficar_comparacion_interaccion_localizacion(resultados_interaccion): #Compara como cambia la localizacion al variar V_eff.
    figura, eje = plt.subplots(figsize=(9, 5))
    for valor_interaccion, resultado in resultados_interaccion:
        indice = indice_estado_mas_cercano_a_cero(resultado["energias"])
        densidad = densidad_estado(resultado["vectores"][:, indice], numero_sitios)
        eje.plot(np.arange(numero_sitios), densidad, linewidth=2.0, label=rf"$V_{{eff}}={valor_interaccion:.1f}$ meV, $\langle\Delta\rangle={resultado['delta_media']:.2f}$ meV")
    eje.set_title("Efecto de la interaccion atractiva en el modo de borde")
    eje.set_xlabel("Sitio i")
    eje.set_ylabel("Densidad normalizada")
    eje.legend()
    figura.tight_layout()
    guardar_figura(figura, "comparacion_interaccion_localizacion_autoconsistente.png")
    return figura


def graficar_perfil_gap_vs_interaccion(resultados_interaccion): #Compara el perfil |Delta_i| al variar V_eff.
    figura, eje = plt.subplots(figsize=(9, 5))
    for valor_interaccion, resultado in resultados_interaccion:
        eje.plot(np.arange(numero_sitios), np.abs(resultado["delta"]), linewidth=2.0, label=rf"$V_{{eff}}={valor_interaccion:.1f}$ meV")
    eje.set_title(r"Perfil autoconsistente $|\Delta_i|$ al variar $V_{eff}$")
    eje.set_xlabel("Sitio i")
    eje.set_ylabel(r"$|\Delta_i|$ (meV)")
    eje.legend()
    figura.tight_layout()
    guardar_figura(figura, "perfil_gap_vs_interaccion_autoconsistente.png")
    return figura


def graficar_mapa_mu_zeeman_minimo(mapa_log, valores_mu, valores_zeeman): #Grafica log10(min|E|) en el plano (mu,E_Z).
    figura, eje = plt.subplots(figsize=(8.5, 5.5))
    imagen = eje.imshow(mapa_log, origin="lower", aspect="auto", extent=[valores_mu[0], valores_mu[-1], valores_zeeman[0], valores_zeeman[-1]], cmap="magma")
    eje.set_title(r"Mapa finito de $\log_{10}(\min |E|)$")
    eje.set_xlabel(r"Potencial quimico $\mu$ (meV)")
    eje.set_ylabel(r"Energia de Zeeman $E_Z$ (meV)")
    barra = figura.colorbar(imagen, ax=eje, pad=0.03)
    barra.set_label(r"$\log_{10}(E/\mathrm{meV})$")
    figura.tight_layout()
    guardar_figura(figura, "mapa_mu_zeeman_minimo_autoconsistente.png")
    return figura


###############################################################################
################################## Ejecucion ##################################
###############################################################################

if __name__ == "__main__":
    preparar_estilo_graficas()

    potencial_prueba = potencial_electrostatico_nulo(numero_sitios)
    hamiltoniano_prueba = construir_hamiltoniano_bdg(numero_sitios, potencial_quimico, zeeman_topologico, delta_semilla, potencial_prueba, hopping, hopping_spin_orbita)
    print("Error hermiticidad:", error_hermiticidad(hamiltoniano_prueba))
    print("Error simetria particula-hueco:", error_simetria_particula_hueco(hamiltoniano_prueba, numero_sitios))

    resultados_regimenes = calcular_resultados_regimenes(numero_sitios)
    resultado_topologico = resultados_regimenes["topologico"]
    resultados_interaccion = calcular_resultados_interaccion(numero_sitios)

    valores_zeeman = np.linspace(0.0, 0.75 * meV, numero_puntos_barrido)
    espectros, energias_minimas, gaps_medios, deltas_finales = barrer_zeeman(numero_sitios, potencial_quimico, valores_zeeman, interaccion_atractiva)

    valores_zeeman_bulk = np.linspace(0.0, 0.75 * meV, 19)
    gaps_bulk, gaps_medios_bulk = barrer_gap_bulk_zeeman(valores_zeeman_bulk, max(60, numero_sitios // 2))

    valores_mu_mapa = np.linspace(-0.8 * meV, 0.8 * meV, 13)
    valores_zeeman_mapa = np.linspace(0.0, 0.8 * meV, 13)
    mapa_log = calcular_mapa_mu_zeeman_minimo(valores_mu_mapa, valores_zeeman_mapa, numero_sitios_mapa)

    graficar_convergencia_gap(resultado_topologico)
    graficar_perfil_gap_autoconsistente(resultados_regimenes)
    graficar_gap_bulk_vs_zeeman(valores_zeeman_bulk, gaps_bulk)
    graficar_bandas_bulk_tres_regimenes(resultados_regimenes)
    graficar_espectro_finito_vs_zeeman(valores_zeeman, espectros, gaps_medios)
    graficar_ldos_tres_regimenes(resultados_regimenes)
    graficar_modo_cero_topologico_y_trivial(resultados_regimenes)
    graficar_particula_hueco_modo_cero(resultado_topologico)
    graficar_majoranas_separadas(resultado_topologico)
    graficar_splitting_vs_longitud()
    graficar_comparacion_interaccion_localizacion(resultados_interaccion)
    graficar_perfil_gap_vs_interaccion(resultados_interaccion)
    graficar_mapa_mu_zeeman_minimo(mapa_log, valores_mu_mapa, valores_zeeman_mapa)

    plt.show()

