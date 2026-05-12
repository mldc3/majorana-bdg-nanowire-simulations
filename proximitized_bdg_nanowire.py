"""
proximitized_bdg_nanowire.py

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
from matplotlib.lines import Line2D


###############################################################################
################ Constantes fisicas y parametros del nanohilo #################
###############################################################################

# Trabajamos en unidades practicas: energia en meV, longitud en nm y temperatura en K.
meV = 1.0
nm = 1.0
kelvin = 1.0

constante_boltzmann = 0.08617333262 * meV / kelvin       # k_B en meV/K.
hbar2_sobre_2m_e = 38.0998212 * meV * nm**2              # hbar^2/(2m_e) en meV nm^2.
magneton_bohr = 0.05788381806 * meV                      # mu_B en meV/T.

# Parametros efectivos de un nanohilo tipo InSb/InAs proximitizado.
masa_efectiva = 0.015                                    # m*/m_e.
factor_g = 40.0                                          # Factor g efectivo.
alpha_rashba = 80.0 * meV * nm                           # Acoplo Rashba alpha_R.

# Malla espacial. Para probar usa numero_sitios_prueba; para figuras finales cambia numero_sitios a numero_sitios_final.
numero_sitios_prueba = 180                                # Valor razonable para comprobar que todo funciona.
numero_sitios_mapa_prueba = 80                           # Tamaño reducido para mapas 2D, porque requieren muchas diagonalizaciones.
numero_sitios = numero_sitios_prueba                     # Cambia esta linea a numero_sitios_final para las figuras finales.
numero_sitios_mapa = numero_sitios_mapa_prueba           # Cambia esta linea a numero_sitios_mapa_final para mapas finales.
paso_espacial = 20.0 * nm
longitud_hilo = (numero_sitios - 1) * paso_espacial

# Escalas discretas que salen de diferencias finitas.
hopping = hbar2_sobre_2m_e / (masa_efectiva * paso_espacial**2)
hopping_spin_orbita = alpha_rashba / (2.0 * paso_espacial)

# Parametros fisicos principales.
potencial_quimico = 0.15 * meV
delta_inducido = 0.25 * meV
temperatura = 0.02 * kelvin

# Perfil inducido por el superconductor: una pequeña zona normal y entrada suave al tramo superconductivo.
sitio_inicio_superconductor = 4
suavidad_borde_superconductor = 2.0

# Barrera electrostatica suave cerca del extremo izquierdo, que representa una compuerta local.
altura_barrera = 0.20 * meV
anchura_barrera = 5.0
centro_barrera = 3.0

# Casos representativos. El criterio ideal homogéneo es E_Z > sqrt(mu^2 + Delta^2).
zeeman_trivial = 0.0 * meV
zeeman_transicion = np.sqrt(potencial_quimico**2 + delta_inducido**2)
zeeman_topologico = 0.55 * meV

# Parametros de barrido y visualizacion.
numero_energias = 220
numero_k_bandas = 320
numero_k_marcadores = 45
ensanchamiento_ldos = 0.025 * meV
numero_puntos_barrido = 17
numero_puntos_mapa = 11

# Valores alternativos para comparar la localizacion al variar el gap inducido y el Rashba.
delta_debil = 0.15 * meV
delta_medio = 0.25 * meV
delta_fuerte = 0.40 * meV
alpha_debil = 30.0 * meV * nm
alpha_medio = 80.0 * meV * nm
alpha_fuerte = 140.0 * meV * nm

# Carpeta de salida para las figuras del TFG.
directorio_figuras = "figurasnanohilo"


def calcular_hoppings_discretos(alpha_rashba_usado, paso_espacial_usado): #Calcula t y t_so para un valor concreto de alpha_R y del paso a.
    hopping_calculado = hbar2_sobre_2m_e / (masa_efectiva * paso_espacial_usado**2)
    hopping_so_calculado = alpha_rashba_usado / (2.0 * paso_espacial_usado)
    return hopping_calculado, hopping_so_calculado


def energia_zeeman_desde_campo(campo_magnetico): #Convierte campo magnetico B en energia de Zeeman E_Z = g mu_B B / 2.
    return 0.5 * factor_g * magneton_bohr * campo_magnetico


def campo_desde_energia_zeeman(energia_zeeman): #Convierte una energia de Zeeman en el campo magnetico equivalente.
    return energia_zeeman / (0.5 * factor_g * magneton_bohr)


###################### MATRICES DE PAULI ######################

matriz_identidad_2 = np.eye(2, dtype=complex)
matriz_sigma_x_2 = np.array([[0.0, 1.0], [1.0, 0.0]], dtype=complex)
matriz_sigma_y_2 = np.array([[0.0, -1.0j], [1.0j, 0.0]], dtype=complex)
matriz_sigma_z_2 = np.array([[1.0, 0.0], [0.0, -1.0]], dtype=complex)
matriz_tau_x_2 = np.array([[0.0, 1.0], [1.0, 0.0]], dtype=complex)
matriz_tau_y_2 = np.array([[0.0, -1.0j], [1.0j, 0.0]], dtype=complex)
matriz_tau_z_2 = np.array([[1.0, 0.0], [0.0, -1.0]], dtype=complex)

# En la base compacta (u_up, u_down, v_down, -v_up) usamos productos tau \otimes sigma.
sigma_x = np.kron(matriz_identidad_2, matriz_sigma_x_2)
sigma_y = np.kron(matriz_identidad_2, matriz_sigma_y_2)
sigma_z = np.kron(matriz_identidad_2, matriz_sigma_z_2)
tau_x = np.kron(matriz_tau_x_2, matriz_identidad_2)
tau_y = np.kron(matriz_tau_y_2, matriz_identidad_2)
tau_z = np.kron(matriz_tau_z_2, matriz_identidad_2)
sigma_y_tau_z = sigma_y @ tau_z
operador_particula_hueco_local = tau_y @ sigma_y


###################### PERFILES ESPACIALES ######################

def posiciones_sitios(numero_sitios_usado, paso_espacial_usado): #Devuelve la posicion x_i de cada sitio de la cadena.
    return np.arange(numero_sitios_usado, dtype=float) * paso_espacial_usado


def convertir_a_perfil_delta(delta, numero_sitios_usado): #Convierte un numero complejo o un array en un perfil Delta_i de longitud N.
    if np.ndim(delta) == 0:
        perfil_delta = np.full(numero_sitios_usado, complex(delta), dtype=complex)
    else:
        perfil_delta = np.asarray(delta, dtype=complex).copy()
    if len(perfil_delta) != numero_sitios_usado:
        raise ValueError("El perfil Delta_i no tiene el mismo numero de sitios que el hilo.")
    return perfil_delta


def perfil_gap_inducido(numero_sitios_usado, delta=delta_inducido, sitio_inicio=sitio_inicio_superconductor, suavidad=suavidad_borde_superconductor, fase=0.0): #Construye un perfil suave para el gap inducido por proximidad.
    indices = np.arange(numero_sitios_usado, dtype=float)
    if sitio_inicio <= 0:
        peso_superconductor = np.ones(numero_sitios_usado, dtype=float)
    else:
        peso_superconductor = 0.5 * (1.0 + np.tanh((indices - sitio_inicio) / suavidad))
    perfil_delta = delta * peso_superconductor * np.exp(1.0j * fase)
    return perfil_delta


def potencial_electrostatico_nulo(numero_sitios_usado): #Devuelve V_i = 0 en todo el hilo.
    return np.zeros(numero_sitios_usado, dtype=float)


def potencial_barrera_gaussiana(numero_sitios_usado, altura=altura_barrera, anchura=anchura_barrera, centro=centro_barrera): #Potencial gaussiano que representa una barrera electrostatica suave.
    indices = np.arange(numero_sitios_usado, dtype=float)
    potencial = altura * np.exp(-0.5 * ((indices - centro) / anchura)**2)
    return potencial


def potencial_dot_suave(numero_sitios_usado, profundidad=0.0 * meV, longitud_dot=6.0, suavidad=2.0): #Potencial suave en el extremo izquierdo, util para simular una region tipo dot.
    indices = np.arange(numero_sitios_usado, dtype=float)
    potencial = 0.5 * profundidad * (np.tanh((indices - longitud_dot) / suavidad) - 1.0)
    return potencial


def potencial_nanohilo_por_defecto(numero_sitios_usado): #Perfil electrostatico usado en las figuras principales del nanohilo.
    return potencial_barrera_gaussiana(numero_sitios_usado, altura=altura_barrera, anchura=anchura_barrera, centro=centro_barrera)


def perfiles_nanohilo_base(numero_sitios_usado): #Devuelve explicitamente el perfil de gap y el potencial usados como caso base.
    perfil_delta = perfil_gap_inducido(numero_sitios_usado, delta=delta_inducido, sitio_inicio=sitio_inicio_superconductor, suavidad=suavidad_borde_superconductor)
    potencial = potencial_nanohilo_por_defecto(numero_sitios_usado)
    return perfil_delta, potencial


###############################################################################
############################## Hamiltoniano BdG ###############################
###############################################################################

def construir_hamiltoniano_bdg(numero_sitios_usado, potencial_quimico_usado, energia_zeeman_usada, perfil_delta, potencial, hopping_usado, hopping_spin_orbita_usado): #Construye la matriz BdG 4N x 4N del nanohilo.
    potencial = np.asarray(potencial, dtype=float).copy()
    perfil_delta = convertir_a_perfil_delta(perfil_delta, numero_sitios_usado)
    if len(potencial) != numero_sitios_usado:
        raise ValueError("El perfil V_i no tiene el mismo numero de sitios que el hilo.")
    hamiltoniano = np.zeros((4 * numero_sitios_usado, 4 * numero_sitios_usado), dtype=complex)
    bloque_hopping = -hopping_usado * tau_z + 1.0j * hopping_spin_orbita_usado * sigma_y_tau_z
    for i in range(numero_sitios_usado):
        bloque_i = slice(4 * i, 4 * (i + 1))
        delta_i = perfil_delta[i]
        bloque_local = ((2.0 * hopping_usado - potencial_quimico_usado + potencial[i]) * tau_z + energia_zeeman_usada * sigma_x + delta_i.real * tau_x - delta_i.imag * tau_y)
        hamiltoniano[bloque_i, bloque_i] = bloque_local
        if i < numero_sitios_usado - 1:
            bloque_j = slice(4 * (i + 1), 4 * (i + 2))
            hamiltoniano[bloque_i, bloque_j] = bloque_hopping
            hamiltoniano[bloque_j, bloque_i] = bloque_hopping.conj().T
    hamiltoniano = 0.5 * (hamiltoniano + hamiltoniano.conj().T)
    return hamiltoniano


def diagonalizar_nanohilo(numero_sitios_usado, potencial_quimico_usado, energia_zeeman_usada, perfil_delta, potencial, hopping_usado, hopping_spin_orbita_usado): #Diagonaliza el Hamiltoniano BdG del nanohilo finito.
    hamiltoniano = construir_hamiltoniano_bdg(numero_sitios_usado, potencial_quimico_usado, energia_zeeman_usada, perfil_delta, potencial, hopping_usado, hopping_spin_orbita_usado)
    energias, vectores = la.eigh(hamiltoniano)
    orden = np.argsort(energias)
    energias = energias[orden]
    vectores = vectores[:, orden]
    return energias, vectores, hamiltoniano




def diagonalizar_energias_nanohilo(numero_sitios_usado, potencial_quimico_usado, energia_zeeman_usada, perfil_delta, potencial, hopping_usado, hopping_spin_orbita_usado): #Calcula solo autovalores cuando no hacen falta autovectores; acelera mucho los barridos.
    hamiltoniano = construir_hamiltoniano_bdg(numero_sitios_usado, potencial_quimico_usado, energia_zeeman_usada, perfil_delta, potencial, hopping_usado, hopping_spin_orbita_usado)
    energias = la.eigvalsh(hamiltoniano)
    energias = np.sort(energias)
    return energias

def error_hermiticidad(hamiltoniano): #Mide numericamente si H = H^dagger.
    numerador = la.norm(hamiltoniano - hamiltoniano.conj().T)
    denominador = max(la.norm(hamiltoniano), 10**(-30))
    return float(numerador / denominador)


def error_simetria_particula_hueco(hamiltoniano, numero_sitios_usado): #Comprueba la simetria BdG C H^* C^dagger = -H.
    operador_global = np.kron(np.eye(numero_sitios_usado, dtype=complex), operador_particula_hueco_local)
    transformado = operador_global @ hamiltoniano.conj() @ operador_global.conj().T
    numerador = la.norm(transformado + hamiltoniano)
    denominador = max(la.norm(hamiltoniano), 10**(-30))
    return float(numerador / denominador)


###############################################################################
############################# Observables fisicos #############################
###############################################################################

def separar_componentes_bdg(vector_estado, numero_sitios_usado): #Separa un autovector en la base compacta (u_up, u_down, v_down, -v_up).
    psi = vector_estado.reshape(numero_sitios_usado, 4)
    u_arriba = psi[:, 0].copy()
    u_abajo = psi[:, 1].copy()
    v_abajo = psi[:, 2].copy()
    menos_v_arriba = psi[:, 3].copy()
    return u_arriba, u_abajo, v_abajo, menos_v_arriba


def densidad_estado(vector_estado, numero_sitios_usado): #Calcula la densidad total |u_up|^2 + |u_down|^2 + |v_up|^2 + |v_down|^2.
    u_arriba, u_abajo, v_abajo, menos_v_arriba = separar_componentes_bdg(vector_estado, numero_sitios_usado)
    densidad = np.abs(u_arriba)**2 + np.abs(u_abajo)**2 + np.abs(v_abajo)**2 + np.abs(menos_v_arriba)**2
    densidad = densidad / max(float(np.sum(densidad)), 10**(-30))
    return densidad


def densidad_electron_y_hueco(vector_estado, numero_sitios_usado): #Devuelve por separado la densidad electronica y la densidad de hueco.
    u_arriba, u_abajo, v_abajo, menos_v_arriba = separar_componentes_bdg(vector_estado, numero_sitios_usado)
    densidad_electron = np.abs(u_arriba)**2 + np.abs(u_abajo)**2
    densidad_hueco = np.abs(v_abajo)**2 + np.abs(menos_v_arriba)**2
    norma = max(float(np.sum(densidad_electron + densidad_hueco)), 10**(-30))
    return densidad_electron / norma, densidad_hueco / norma


def peso_en_bordes(vector_estado, numero_sitios_usado, fraccion_borde=0.20): #Calcula qué fraccion de la densidad vive cerca de los dos extremos.
    densidad = densidad_estado(vector_estado, numero_sitios_usado)
    ancho = max(1, int(fraccion_borde * numero_sitios_usado))
    peso = float(np.sum(densidad[:ancho]) + np.sum(densidad[-ancho:]))
    return peso


def indice_estado_mas_cercano_a_cero(energias): #Encuentra el indice del autovalor mas cercano a E = 0.
    return int(np.argmin(np.abs(energias)))


def energia_minima_absoluta(energias): #Devuelve min |E|, util para detectar estados casi cero.
    return float(np.min(np.abs(energias)))


def transformar_particula_hueco(vector_estado, numero_sitios_usado): #Aplica la simetria particula-hueco al estado BdG.
    psi = vector_estado.reshape(numero_sitios_usado, 4)
    psi_ph = np.zeros_like(psi)
    for i in range(numero_sitios_usado):
        psi_ph[i, :] = operador_particula_hueco_local @ np.conj(psi[i, :])
    return psi_ph.reshape(4 * numero_sitios_usado)


def alinear_fase(estado_referencia, estado_objetivo): #Elimina la fase global arbitraria entre dos autovectores.
    solapamiento = np.vdot(estado_referencia, estado_objetivo)
    if abs(solapamiento) < 10**(-30):
        estado_alineado = estado_objetivo.copy()
    else:
        estado_alineado = estado_objetivo * np.exp(-1.0j * np.angle(solapamiento))
    return estado_alineado


def construir_majoranas_desde_pareja(energias, vectores, numero_sitios_usado): #Construye dos combinaciones tipo Majorana desde la pareja pm E mas cercana a cero.
    indices_positivos = np.where(energias >= 0.0)[0]
    indices_negativos = np.where(energias < 0.0)[0]
    indice_positivo = indices_positivos[np.argmin(np.abs(energias[indices_positivos]))]
    indice_negativo = indices_negativos[np.argmin(np.abs(energias[indices_negativos] + energias[indice_positivo]))]
    estado_positivo = vectores[:, indice_positivo]
    estado_negativo = vectores[:, indice_negativo]
    estado_phs = transformar_particula_hueco(estado_positivo, numero_sitios_usado)
    estado_negativo = alinear_fase(estado_phs, estado_negativo)
    majorana_1 = (estado_positivo + estado_negativo) / np.sqrt(2.0)
    majorana_2 = -1.0j * (estado_positivo - estado_negativo) / np.sqrt(2.0)
    densidad_1 = densidad_estado(majorana_1, numero_sitios_usado)
    densidad_2 = densidad_estado(majorana_2, numero_sitios_usado)
    mitad = numero_sitios_usado // 2
    peso_izquierdo_1 = np.sum(densidad_1[:mitad])
    peso_izquierdo_2 = np.sum(densidad_2[:mitad])
    if peso_izquierdo_1 >= peso_izquierdo_2:
        densidad_izquierda = densidad_1
        densidad_derecha = densidad_2
    else:
        densidad_izquierda = densidad_2
        densidad_derecha = densidad_1
    return densidad_izquierda, densidad_derecha


def lorentziana(energia, centro, anchura): #Funcion lorentziana usada para suavizar la densidad local de estados.
    return anchura / (np.pi * ((energia - centro)**2 + anchura**2))


def calcular_ldos(numero_sitios_usado, energias, vectores, energias_ldos, anchura): #Calcula la LDOS electronica: |u|^2 L(E-En) + |v|^2 L(E+En) para En > 0.
    ldos = np.zeros((len(energias_ldos), numero_sitios_usado), dtype=float)
    for n, energia_n in enumerate(energias):
        if energia_n <= 0.0:
            continue
        densidad_electron, densidad_hueco = densidad_electron_y_hueco(vectores[:, n], numero_sitios_usado)
        for i, energia in enumerate(energias_ldos):
            ldos[i, :] += densidad_electron * lorentziana(energia, energia_n, anchura) + densidad_hueco * lorentziana(energia, -energia_n, anchura)
    return ldos


def es_topologico_ideal(potencial_quimico_usado, energia_zeeman_usada, delta_usada, alpha_usado): #Criterio homogéneo orientativo del nanohilo.
    condicion_zeeman = energia_zeeman_usada > np.sqrt(potencial_quimico_usado**2 + np.abs(delta_usada)**2)
    condicion_rashba = abs(alpha_usado) > 10**(-12)
    return bool(condicion_zeeman and condicion_rashba)


###############################################################################
########################## Bandas bulk y barridos #############################
###############################################################################

def matriz_bulk_k(punto_k, potencial_quimico_usado, energia_zeeman_usada, delta_uniforme, hopping_usado, hopping_spin_orbita_usado): #Hamiltoniano 4x4 del nanohilo homogeneo en espacio k.
    xi_k = 2.0 * hopping_usado - potencial_quimico_usado - 2.0 * hopping_usado * np.cos(punto_k)
    rashba_k = -2.0 * hopping_spin_orbita_usado * np.sin(punto_k)
    hamiltoniano_k = xi_k * tau_z + rashba_k * sigma_y_tau_z + energia_zeeman_usada * sigma_x + np.real(delta_uniforme) * tau_x - np.imag(delta_uniforme) * tau_y
    return hamiltoniano_k


def bandas_bulk_diagonalizando(potencial_quimico_usado, energia_zeeman_usada, delta_uniforme, puntos_k, hopping_usado, hopping_spin_orbita_usado): #Bandas obtenidas diagonalizando H(k) punto a punto.
    bandas = np.zeros((len(puntos_k), 4), dtype=float)
    for i, punto_k in enumerate(puntos_k):
        bandas[i, :] = la.eigvalsh(matriz_bulk_k(punto_k, potencial_quimico_usado, energia_zeeman_usada, delta_uniforme, hopping_usado, hopping_spin_orbita_usado))
    return bandas


def bandas_bulk_teoricas(potencial_quimico_usado, energia_zeeman_usada, delta_uniforme, puntos_k, hopping_usado, hopping_spin_orbita_usado): #Bandas analiticas del nanohilo homogeneo en red.
    delta_abs = np.abs(delta_uniforme)
    bandas = np.zeros((len(puntos_k), 4), dtype=float)
    for i, punto_k in enumerate(puntos_k):
        xi_k = 2.0 * hopping_usado - potencial_quimico_usado - 2.0 * hopping_usado * np.cos(punto_k)
        rashba_k = -2.0 * hopping_spin_orbita_usado * np.sin(punto_k)
        A = xi_k**2 + rashba_k**2 + energia_zeeman_usada**2 + delta_abs**2
        B = 2.0 * np.sqrt(energia_zeeman_usada**2 * (xi_k**2 + delta_abs**2) + xi_k**2 * rashba_k**2)
        energia_menor = np.sqrt(max(A - B, 0.0))
        energia_mayor = np.sqrt(max(A + B, 0.0))
        bandas[i, :] = np.array([-energia_mayor, -energia_menor, energia_menor, energia_mayor])
    return bandas


def gap_bulk_minimo(potencial_quimico_usado, energia_zeeman_usada, delta_uniforme, hopping_usado, hopping_spin_orbita_usado): #Calcula el minimo de las bandas positivas del sistema bulk.
    puntos_k = np.linspace(-np.pi, np.pi, numero_k_bandas)
    bandas = bandas_bulk_diagonalizando(potencial_quimico_usado, energia_zeeman_usada, delta_uniforme, puntos_k, hopping_usado, hopping_spin_orbita_usado)
    bandas_positivas = bandas[bandas > 0.0]
    return float(np.min(bandas_positivas))


def barrer_zeeman(numero_sitios_usado, potencial_quimico_usado, valores_zeeman, perfil_delta, potencial, hopping_usado, hopping_spin_orbita_usado): #Barre E_Z y calcula el espectro finito y min|E|.
    todos_los_espectros = []
    energias_minimas = []
    for energia_zeeman in valores_zeeman:
        energias = diagonalizar_energias_nanohilo(numero_sitios_usado, potencial_quimico_usado, energia_zeeman, perfil_delta, potencial, hopping_usado, hopping_spin_orbita_usado)
        todos_los_espectros.append(energias)
        energias_minimas.append(energia_minima_absoluta(energias))
    return np.array(todos_los_espectros), np.array(energias_minimas)


def barrer_potencial_quimico(numero_sitios_usado, valores_mu, energia_zeeman_usada, perfil_delta, potencial, hopping_usado, hopping_spin_orbita_usado): #Barre el potencial quimico y calcula el espectro del nanohilo finito.
    todos_los_espectros = []
    energias_minimas = []
    for mu in valores_mu:
        energias = diagonalizar_energias_nanohilo(numero_sitios_usado, mu, energia_zeeman_usada, perfil_delta, potencial, hopping_usado, hopping_spin_orbita_usado)
        todos_los_espectros.append(energias)
        energias_minimas.append(energia_minima_absoluta(energias))
    return np.array(todos_los_espectros), np.array(energias_minimas)


def calcular_diagrama_fase_mu_zeeman(valores_mu, valores_zeeman, delta_usada, alpha_usado): #Diagrama ideal del nanohilo homogeneo en el plano (mu, E_Z).
    diagrama = np.zeros((len(valores_zeeman), len(valores_mu)), dtype=float)
    for i, energia_zeeman in enumerate(valores_zeeman):
        for j, mu in enumerate(valores_mu):
            diagrama[i, j] = 1.0 if es_topologico_ideal(mu, energia_zeeman, delta_usada, alpha_usado) else 0.0
    return diagrama


def barrer_rashba(valores_alpha, numero_sitios_usado, perfil_delta, potencial): #Barre alpha_R y calcula min|E| y peso de borde.
    energias_minimas = []
    pesos_bordes = []
    for alpha in valores_alpha:
        hopping_usado, hopping_spin_orbita_usado = calcular_hoppings_discretos(alpha, paso_espacial)
        energias, vectores, hamiltoniano = diagonalizar_nanohilo(numero_sitios_usado, potencial_quimico, zeeman_topologico, perfil_delta, potencial, hopping_usado, hopping_spin_orbita_usado)
        indice = indice_estado_mas_cercano_a_cero(energias)
        energias_minimas.append(energia_minima_absoluta(energias))
        pesos_bordes.append(peso_en_bordes(vectores[:, indice], numero_sitios_usado))
    return np.array(energias_minimas), np.array(pesos_bordes)


def calcular_diagrama_alpha_zeeman(valores_alpha, valores_zeeman, numero_sitios_usado, perfil_delta, potencial): #Mapa de min|E| y criterio ideal en el plano (alpha_R, E_Z).
    mapa_log = np.zeros((len(valores_zeeman), len(valores_alpha)), dtype=float)
    mapa_ideal = np.zeros_like(mapa_log)
    for i, energia_zeeman in enumerate(valores_zeeman):
        for j, alpha in enumerate(valores_alpha):
            hopping_usado, hopping_spin_orbita_usado = calcular_hoppings_discretos(alpha, paso_espacial)
            energias = diagonalizar_energias_nanohilo(numero_sitios_usado, potencial_quimico, energia_zeeman, perfil_delta, potencial, hopping_usado, hopping_spin_orbita_usado)
            mapa_log[i, j] = np.log10(max(energia_minima_absoluta(energias) / meV, 10**(-16)))
            mapa_ideal[i, j] = 1.0 if es_topologico_ideal(potencial_quimico, energia_zeeman, delta_inducido, alpha) else 0.0
    return mapa_log, mapa_ideal


def barrer_paso_malla(valores_paso, longitud_objetivo, energia_zeeman_usada): #Comprueba cómo cambia min|E| al modificar el paso espacial manteniendo la longitud fisica aproximada.
    energias_minimas = []
    numeros_sitios = []
    for paso in valores_paso:
        N = int(round(longitud_objetivo / paso)) + 1
        hopping_usado, hopping_spin_orbita_usado = calcular_hoppings_discretos(alpha_rashba, paso)
        perfil_delta = perfil_gap_inducido(N, delta=delta_inducido, sitio_inicio=sitio_inicio_superconductor, suavidad=suavidad_borde_superconductor)
        potencial = potencial_nanohilo_por_defecto(N)
        energias = diagonalizar_energias_nanohilo(N, potencial_quimico, energia_zeeman_usada, perfil_delta, potencial, hopping_usado, hopping_spin_orbita_usado)
        energias_minimas.append(energia_minima_absoluta(energias))
        numeros_sitios.append(N)
    return np.array(energias_minimas), np.array(numeros_sitios)


###############################################################################
################################## Graficas ###################################
###############################################################################

def preparar_estilo_graficas(): #Estilo sencillo para que las figuras sean claras en el TFG.
    plt.style.use("default")
    plt.rcParams.update({"figure.figsize": (9, 5.5), "axes.grid": True, "grid.alpha": 0.23, "font.size": 11, "axes.titlesize": 13, "axes.labelsize": 12, "legend.fontsize": 9, "savefig.dpi": 300})


def guardar_figura(figura, nombre_archivo): #Guarda la figura en la carpeta de salida con margenes ajustados.
    os.makedirs(directorio_figuras, exist_ok=True)
    ruta = os.path.join(directorio_figuras, nombre_archivo)
    figura.savefig(ruta, dpi=300, bbox_inches="tight")
    return ruta


def graficar_bandas_bulk_tres_regimenes(): #Grafica bandas teoricas y bandas obtenidas diagonalizando H(k).
    casos = [("trivial", zeeman_trivial), ("transición", zeeman_transicion), ("topológico", zeeman_topologico)]
    puntos_k = np.linspace(-np.pi, np.pi, numero_k_bandas)
    puntos_k_marcadores = np.linspace(-np.pi, np.pi, numero_k_marcadores)
    figura, ejes = plt.subplots(1, 3, figsize=(13.2, 4.0), sharey=True)
    for eje, (nombre, energia_zeeman) in zip(ejes, casos):
        bandas_teoricas = bandas_bulk_teoricas(potencial_quimico, energia_zeeman, delta_inducido, puntos_k, hopping, hopping_spin_orbita)
        bandas_diag = bandas_bulk_diagonalizando(potencial_quimico, energia_zeeman, delta_inducido, puntos_k_marcadores, hopping, hopping_spin_orbita)
        for indice_banda in range(4):
            eje.plot(puntos_k / np.pi, bandas_teoricas[:, indice_banda] / meV, linewidth=1.9)
            eje.scatter(puntos_k_marcadores / np.pi, bandas_diag[:, indice_banda] / meV, s=10, marker="x", alpha=0.75)
        eje.axhline(0.0, color="black", linewidth=1.0, alpha=0.55)
        eje.set_title(f"{nombre}: $E_Z={energia_zeeman / meV:.2f}$ meV")
        eje.set_xlabel(r"$ka/\pi$")
        leyenda_bandas = [Line2D([0], [0], linewidth=1.9, label="expresión analítica"), Line2D([0], [0], marker="x", linestyle="", markersize=6, label="diag. de $H(k)$")]
        eje.legend(handles=leyenda_bandas, loc="upper center")
    ejes[0].set_ylabel("Energía (meV)")
    figura.suptitle("Bandas bulk del nanohilo BdG")
    figura.tight_layout()
    guardar_figura(figura, "nanohilo_realista_bandas_bulk.png")
    return figura


def graficar_espectro_finito_vs_zeeman(datos_barrido): #Grafica el espectro finito frente a E_Z.
    valores_zeeman, espectros, energias_minimas = datos_barrido
    figura, eje = plt.subplots(figsize=(9.6, 5.3))
    for i, energia_zeeman in enumerate(valores_zeeman):
        energias_cerca_cero = espectros[i][np.abs(espectros[i]) < 1.2 * meV]
        eje.scatter(np.full_like(energias_cerca_cero, energia_zeeman / meV), energias_cerca_cero / meV, s=9, color="black", alpha=0.62)
    eje.axvline(zeeman_transicion / meV, linestyle="--", linewidth=2.0, label="criterio ideal")
    eje.axhline(0.0, color="black", linewidth=1.0, alpha=0.50)
    eje.set_title("Espectro finito del nanohilo frente a Zeeman")
    eje.set_xlabel(r"Energía de Zeeman $E_Z$ (meV)")
    eje.set_ylabel("Energía BdG (meV)")
    eje.legend()
    figura.tight_layout()
    guardar_figura(figura, "nanohilo_realista_espectro_zeeman.png")
    return figura


def graficar_espectro_finito_vs_mu(datos_barrido): #Grafica el espectro finito frente a mu, con E_Z fijo en el regimen topologico.
    valores_mu, espectros, energias_minimas = datos_barrido
    figura, eje = plt.subplots(figsize=(9.6, 5.3))
    for i, mu in enumerate(valores_mu):
        energias_cerca_cero = espectros[i][np.abs(espectros[i]) < 1.2 * meV]
        eje.scatter(np.full_like(energias_cerca_cero, mu / meV), energias_cerca_cero / meV, s=9, color="black", alpha=0.62)
    mu_critico = np.sqrt(max(zeeman_topologico**2 - delta_inducido**2, 0.0))
    eje.axvline(-mu_critico / meV, linestyle="--", linewidth=2.0, label="frontera ideal")
    eje.axvline(+mu_critico / meV, linestyle="--", linewidth=2.0)
    eje.axhline(0.0, color="black", linewidth=1.0, alpha=0.50)
    eje.set_title("Espectro finito frente al potencial químico")
    eje.set_xlabel(r"Potencial químico $\mu$ (meV)")
    eje.set_ylabel("Energía BdG (meV)")
    eje.legend()
    figura.tight_layout()
    guardar_figura(figura, "nanohilo_realista_espectro_mu.png")
    return figura


def graficar_diagrama_fase_mu_zeeman(): #Grafica el diagrama ideal de fase en el plano (mu, E_Z).
    valores_mu = np.linspace(-0.8 * meV, 0.8 * meV, 220)
    valores_zeeman = np.linspace(0.0, 0.9 * meV, 220)
    diagrama = calcular_diagrama_fase_mu_zeeman(valores_mu, valores_zeeman, delta_inducido, alpha_rashba)
    figura, eje = plt.subplots(figsize=(7.7, 5.3))
    imagen = eje.imshow(diagrama, origin="lower", aspect="auto", extent=[valores_mu[0] / meV, valores_mu[-1] / meV, valores_zeeman[0] / meV, valores_zeeman[-1] / meV], cmap="viridis")
    eje.set_title(r"Diagrama ideal del nanohilo: $(\mu,E_Z)$")
    eje.set_xlabel(r"Potencial químico $\mu$ (meV)")
    eje.set_ylabel(r"Energía de Zeeman $E_Z$ (meV)")
    barra = figura.colorbar(imagen, ax=eje, pad=0.025)
    barra.set_label("0 = trivial, 1 = topológico")
    figura.tight_layout()
    guardar_figura(figura, "nanohilo_realista_diagrama_mu_zeeman.png")
    return figura


def graficar_ldos_tres_regimenes(perfil_delta, potencial): #Grafica LDOS(E,i) para regimen trivial, transicion y topologico.
    casos = [("trivial", zeeman_trivial), ("transición", zeeman_transicion), ("topológico", zeeman_topologico)]
    energias_ldos = np.linspace(-1.0 * meV, 1.0 * meV, numero_energias)
    figura, ejes = plt.subplots(1, 3, figsize=(13.6, 4.15), sharey=True)
    for eje, (nombre, energia_zeeman) in zip(ejes, casos):
        energias, vectores, hamiltoniano = diagonalizar_nanohilo(numero_sitios, potencial_quimico, energia_zeeman, perfil_delta, potencial, hopping, hopping_spin_orbita)
        ldos = calcular_ldos(numero_sitios, energias, vectores, energias_ldos, ensanchamiento_ldos)
        imagen = eje.imshow(ldos, aspect="auto", origin="lower", extent=[0, numero_sitios - 1, energias_ldos[0] / meV, energias_ldos[-1] / meV], cmap="magma")
        eje.axhline(0.0, color="white", linewidth=1.0, alpha=0.75)
        eje.set_title(f"{nombre}: $E_Z={energia_zeeman / meV:.2f}$ meV")
        eje.set_xlabel("Sitio i")
    ejes[0].set_ylabel("Energía (meV)")
    figura.suptitle(r"LDOS$(E,i)$ del nanohilo")
    figura.subplots_adjust(left=0.065, right=0.875, bottom=0.16, top=0.80, wspace=0.08)
    eje_barra = figura.add_axes([0.905, 0.16, 0.017, 0.64])
    figura.colorbar(imagen, cax=eje_barra, label="LDOS")
    guardar_figura(figura, "nanohilo_realista_ldos.png")
    return figura


def graficar_modo_cero_topologico_y_trivial(perfil_delta, potencial): #Compara el estado mas cercano a cero en un caso trivial y uno topologico.
    casos = [("trivial", zeeman_trivial), ("topológico", zeeman_topologico)]
    figura, ejes = plt.subplots(1, 2, figsize=(12.0, 4.1), sharey=True)
    for eje, (nombre, energia_zeeman) in zip(ejes, casos):
        energias, vectores, hamiltoniano = diagonalizar_nanohilo(numero_sitios, potencial_quimico, energia_zeeman, perfil_delta, potencial, hopping, hopping_spin_orbita)
        indice = indice_estado_mas_cercano_a_cero(energias)
        densidad = densidad_estado(vectores[:, indice], numero_sitios)
        eje.plot(np.arange(numero_sitios), densidad, linewidth=2.25)
        eje.set_title(f"{nombre}: E = {energias[indice] / meV:.3e} meV")
        eje.set_xlabel("Sitio i")
    ejes[0].set_ylabel("Densidad normalizada")
    figura.suptitle("Estado BdG más cercano a energía cero")
    figura.tight_layout()
    guardar_figura(figura, "nanohilo_realista_estado_cero.png")
    return figura


def graficar_majoranas_separadas(perfil_delta, potencial): #Grafica las dos combinaciones tipo Majorana del modo casi cero.
    energias, vectores, hamiltoniano = diagonalizar_nanohilo(numero_sitios, potencial_quimico, zeeman_topologico, perfil_delta, potencial, hopping, hopping_spin_orbita)
    densidad_izquierda, densidad_derecha = construir_majoranas_desde_pareja(energias, vectores, numero_sitios)
    figura, eje = plt.subplots(figsize=(9.0, 5.2))
    eje.plot(np.arange(numero_sitios), densidad_izquierda, linewidth=2.35, label="Majorana izquierda")
    eje.plot(np.arange(numero_sitios), densidad_derecha, linewidth=2.35, label="Majorana derecha")
    eje.set_title("Descomposición del modo casi cero en dos Majoranas")
    eje.set_xlabel("Sitio i")
    eje.set_ylabel("Densidad normalizada")
    eje.legend()
    figura.tight_layout()
    guardar_figura(figura, "nanohilo_realista_majoranas.png")
    return figura


def graficar_particula_hueco_modo_cero(perfil_delta, potencial): #Grafica contenido electronico y de hueco del estado mas cercano a cero.
    energias, vectores, hamiltoniano = diagonalizar_nanohilo(numero_sitios, potencial_quimico, zeeman_topologico, perfil_delta, potencial, hopping, hopping_spin_orbita)
    indice = indice_estado_mas_cercano_a_cero(energias)
    densidad_electron, densidad_hueco = densidad_electron_y_hueco(vectores[:, indice], numero_sitios)
    figura, eje = plt.subplots(figsize=(9.0, 5.2))
    eje.plot(np.arange(numero_sitios), densidad_electron, linewidth=2.25, label="parte electrónica")
    eje.plot(np.arange(numero_sitios), densidad_hueco, linewidth=2.25, linestyle="--", label="parte de hueco")
    eje.set_title("Contenido partícula-hueco del modo cercano a cero")
    eje.set_xlabel("Sitio i")
    eje.set_ylabel("Densidad normalizada")
    eje.legend()
    figura.tight_layout()
    guardar_figura(figura, "nanohilo_realista_electron_hueco.png")
    return figura


def graficar_barrido_rashba(): #Grafica min|E| y peso en bordes al variar alpha_R.
    valores_alpha = np.linspace(0.0, 160.0 * meV * nm, 21)
    perfil_delta_mapa, potencial_mapa = perfiles_nanohilo_base(numero_sitios_mapa)
    energias_minimas, pesos_bordes = barrer_rashba(valores_alpha, numero_sitios_mapa, perfil_delta_mapa, potencial_mapa)
    figura, eje_1 = plt.subplots(figsize=(9.0, 5.4))
    eje_2 = eje_1.twinx()
    eje_1.plot(valores_alpha / (meV * nm), energias_minimas / meV, linewidth=2.2, label=r"$\min |E|$")
    eje_2.plot(valores_alpha / (meV * nm), pesos_bordes, linewidth=2.2, color="tab:orange", label="peso en bordes")
    eje_1.set_title("Barrido del acoplo Rashba")
    eje_1.set_xlabel(r"$\alpha_R$ (meV nm)")
    eje_1.set_ylabel(r"$\min |E|$ (meV)")
    eje_2.set_ylabel("Peso en bordes", color="tab:orange")
    lineas = eje_1.get_lines() + eje_2.get_lines()
    etiquetas = [linea.get_label() for linea in lineas]
    eje_1.legend(lineas, etiquetas, loc="upper right")
    figura.tight_layout()
    guardar_figura(figura, "nanohilo_realista_barrido_rashba.png")
    return figura


def graficar_diagrama_alpha_zeeman(): #Grafica un mapa en el plano (alpha_R, E_Z).
    valores_alpha = np.linspace(0.0, 160.0 * meV * nm, numero_puntos_mapa)
    valores_zeeman = np.linspace(0.0, 0.85 * meV, numero_puntos_mapa)
    perfil_delta_mapa, potencial_mapa = perfiles_nanohilo_base(numero_sitios_mapa)
    mapa_log, mapa_ideal = calcular_diagrama_alpha_zeeman(valores_alpha, valores_zeeman, numero_sitios_mapa, perfil_delta_mapa, potencial_mapa)
    figura, ejes = plt.subplots(1, 2, figsize=(12.0, 4.6), sharey=True)
    imagen_1 = ejes[0].imshow(mapa_log, origin="lower", aspect="auto", extent=[valores_alpha[0] / (meV * nm), valores_alpha[-1] / (meV * nm), valores_zeeman[0] / meV, valores_zeeman[-1] / meV], cmap="magma")
    imagen_2 = ejes[1].imshow(mapa_ideal, origin="lower", aspect="auto", extent=[valores_alpha[0] / (meV * nm), valores_alpha[-1] / (meV * nm), valores_zeeman[0] / meV, valores_zeeman[-1] / meV], cmap="viridis", vmin=0.0, vmax=1.0)
    ejes[0].set_title(r"$\log_{10}(\min |E|)$")
    ejes[1].set_title("Criterio ideal")
    for eje in ejes:
        eje.set_xlabel(r"$\alpha_R$ (meV nm)")
    ejes[0].set_ylabel(r"$E_Z$ (meV)")
    figura.suptitle(r"Diagrama en el plano $(\alpha_R,E_Z)$")
    figura.subplots_adjust(left=0.07, right=0.88, bottom=0.15, top=0.82, wspace=0.10)
    eje_barra_1 = figura.add_axes([0.47, 0.15, 0.017, 0.64])
    eje_barra_2 = figura.add_axes([0.91, 0.15, 0.017, 0.64])
    figura.colorbar(imagen_1, cax=eje_barra_1, label=r"$\log_{10}(E/\mathrm{meV})$")
    figura.colorbar(imagen_2, cax=eje_barra_2, label="0 = trivial, 1 = topológico")
    guardar_figura(figura, "nanohilo_realista_diagrama_alpha_zeeman.png")
    return figura


def graficar_convergencia_paso_malla(): #Grafica una comprobacion sencilla de convergencia con el paso espacial.
    valores_paso = np.array([30.0, 25.0, 20.0, 15.0]) * nm
    longitud_objetivo = (numero_sitios - 1) * paso_espacial
    energias_minimas, numeros_sitios = barrer_paso_malla(valores_paso, longitud_objetivo, zeeman_topologico)
    figura, eje = plt.subplots(figsize=(8.2, 5.1))
    eje.plot(valores_paso / nm, energias_minimas / meV, marker="o", linewidth=2.25)
    for paso, energia, N in zip(valores_paso, energias_minimas, numeros_sitios):
        eje.annotate(f"N={N}", (paso / nm, energia / meV), textcoords="offset points", xytext=(6, 6))
    eje.invert_xaxis()
    eje.set_title("Convergencia al refinar la malla espacial")
    eje.set_xlabel("Paso espacial a (nm)")
    eje.set_ylabel(r"$\min |E|$ (meV)")
    figura.tight_layout()
    guardar_figura(figura, "nanohilo_realista_convergencia_malla.png")
    return figura


def graficar_comparacion_gap_localizacion(perfil_delta_base, potencial_base): #Compara cómo cambia la localizacion al variar el gap inducido.
    valores_delta = [("gap débil", delta_debil), ("gap medio", delta_medio), ("gap fuerte", delta_fuerte)]
    figura, eje = plt.subplots(figsize=(9.0, 5.2))
    for nombre, delta in valores_delta:
        perfil_delta = perfil_gap_inducido(numero_sitios, delta=delta, sitio_inicio=sitio_inicio_superconductor, suavidad=suavidad_borde_superconductor)
        energias, vectores, hamiltoniano = diagonalizar_nanohilo(numero_sitios, potencial_quimico, zeeman_topologico, perfil_delta, potencial_base, hopping, hopping_spin_orbita)
        indice = indice_estado_mas_cercano_a_cero(energias)
        densidad = densidad_estado(vectores[:, indice], numero_sitios)
        eje.plot(np.arange(numero_sitios), densidad, linewidth=2.0, label=f"{nombre}: Delta={delta / meV:.2f} meV")
    eje.set_title("Efecto del gap inducido en la localización")
    eje.set_xlabel("Sitio i")
    eje.set_ylabel("Densidad normalizada")
    eje.legend()
    figura.tight_layout()
    guardar_figura(figura, "nanohilo_realista_gap_localizacion.png")
    return figura


def generar_todas_las_graficas(): #Ejecuta las figuras que aparecen en los apuntes de la seccion del nanohilo realista.
    preparar_estilo_graficas()
    perfil_delta, potencial = perfiles_nanohilo_base(numero_sitios)
    hamiltoniano_prueba = construir_hamiltoniano_bdg(numero_sitios, potencial_quimico, zeeman_topologico, perfil_delta, potencial, hopping, hopping_spin_orbita)
    print("N usado en las figuras principales:", numero_sitios)
    print("Dimensión de la matriz BdG:", 4 * numero_sitios)
    print("Para resultado final cambia numero_sitios = numero_sitios_final y, si quieres, numero_sitios_mapa = numero_sitios_mapa_final.")
    print("Error hermiticidad:", error_hermiticidad(hamiltoniano_prueba))
    print("Error simetria particula-hueco:", error_simetria_particula_hueco(hamiltoniano_prueba, numero_sitios))
    print("Campo equivalente al caso topologico [T]:", campo_desde_energia_zeeman(zeeman_topologico))
    valores_zeeman = np.linspace(0.0, 0.75 * meV, numero_puntos_barrido)
    espectros_zeeman, energias_minimas_zeeman = barrer_zeeman(numero_sitios, potencial_quimico, valores_zeeman, perfil_delta, potencial, hopping, hopping_spin_orbita)
    datos_barrido_zeeman = (valores_zeeman, espectros_zeeman, energias_minimas_zeeman)
    valores_mu = np.linspace(-0.8 * meV, 0.8 * meV, numero_puntos_barrido)
    espectros_mu, energias_minimas_mu = barrer_potencial_quimico(numero_sitios, valores_mu, zeeman_topologico, perfil_delta, potencial, hopping, hopping_spin_orbita)
    datos_barrido_mu = (valores_mu, espectros_mu, energias_minimas_mu)
    figuras = {}
    figuras["bandas_bulk"] = graficar_bandas_bulk_tres_regimenes()
    figuras["espectro_zeeman"] = graficar_espectro_finito_vs_zeeman(datos_barrido_zeeman)
    figuras["espectro_mu"] = graficar_espectro_finito_vs_mu(datos_barrido_mu)
    figuras["diagrama_mu_zeeman"] = graficar_diagrama_fase_mu_zeeman()
    figuras["ldos"] = graficar_ldos_tres_regimenes(perfil_delta, potencial)
    figuras["estado_cero"] = graficar_modo_cero_topologico_y_trivial(perfil_delta, potencial)
    figuras["majoranas"] = graficar_majoranas_separadas(perfil_delta, potencial)
    figuras["electron_hueco"] = graficar_particula_hueco_modo_cero(perfil_delta, potencial)
    figuras["barrido_rashba"] = graficar_barrido_rashba()
    figuras["diagrama_alpha_zeeman"] = graficar_diagrama_alpha_zeeman()
    figuras["convergencia_malla"] = graficar_convergencia_paso_malla()
    figuras["gap_localizacion"] = graficar_comparacion_gap_localizacion(perfil_delta, potencial)
    return figuras


###############################################################################
############################### Ejecucion final ###############################
###############################################################################

if __name__ == "__main__":
    figuras = generar_todas_las_graficas()
    plt.show()






"""
import numpy as np
import scipy.linalg as la
import matplotlib.pyplot as plt

###############################################################################
################ Constantes fisicas y parametros del nanohilo #################
###############################################################################

# Trabajamos en unidades practicas: energia en meV y longitud en nm
meV = 1.0
nm = 1.0
kelvin = 1.0

constante_boltzmann = 0.08617333262 * meV / kelvin       # k_B en meV/K.
hbar2_sobre_2m_e = 38.0998212 * meV * nm**2              # hbar^2/(2m_e) en meV nm^2.
magneton_bohr = 0.05788381806 * meV                      # mu_B en meV/T.

# Parámetros efectivos de un nanohilo tipo InSb/InAs proximitizado.
masa_efectiva = 0.015                                    # m*/m_e.
factor_g = 40.0                                          # Factor g efectivo.
alpha_rashba = 80.0 * meV * nm                           # Acoplo Rashba alpha_R.

# Malla espacial. El código sigue siendo denso para que se lea claramente.
numero_sitios = 60
paso_espacial = 20.0 * nm
longitud_hilo = (numero_sitios - 1) * paso_espacial

# Escalas discretas que salen de diferencias finitas.
hopping = hbar2_sobre_2m_e / (masa_efectiva * paso_espacial**2)
hopping_spin_orbita = alpha_rashba / (2.0 * paso_espacial)

# Parámetros físicos principales.
potencial_quimico = 0.15 * meV
delta_inducido = 0.25 * meV
temperatura = 0.02 * kelvin

# En el nanohilo el superconductor puede cubrir todo el hilo o empezar después
# de una pequeña zona normal. Para la mayoría de figuras usamos una cobertura
# casi completa, pero dejamos estos parámetros explícitos.
sitio_inicio_superconductor = 4
suavidad_borde_superconductor = 2.0

# Barrera suave cerca del extremo izquierdo, parecida a una compuerta local.
altura_barrera = 0.20 * meV
anchura_barrera = 5.0
centro_barrera = 3.0

# Casos representativos. El criterio ideal homogéneo es
# E_Z > sqrt(mu^2 + Delta^2), así que el valor intermedio marca la transición.
zeeman_trivial = 0.0 * meV
zeeman_transicion = np.sqrt(potencial_quimico**2 + delta_inducido**2)
zeeman_topologico = 0.55 * meV

# Parámetros para la LDOS y los barridos.
numero_energias = 320
ensanchamiento_ldos = 0.025 * meV
numero_puntos_barrido = 35

# Valores alternativos para comparar física.
delta_debil = 0.15 * meV
delta_medio = 0.25 * meV
delta_fuerte = 0.40 * meV

alpha_debil = 30.0 * meV * nm
alpha_medio = 80.0 * meV * nm
alpha_fuerte = 140.0 * meV * nm

def calcular_hoppings_discretos(alpha_rashba_usado, paso_espacial_usado): #Calcula t y t_so para un valor concreto de alpha_R y del paso a
    return hbar2_sobre_2m_e / (masa_efectiva * paso_espacial_usado**2), alpha_rashba_usado / (2.0 * paso_espacial_usado)

def energia_zeeman_desde_campo(campo_magnetico): #Convierte campo magnético B en energía de Zeeman E_Z = g mu_B B / 2
    return 0.5 * factor_g * magneton_bohr * campo_magnetico

def campo_desde_energia_zeeman(energia_zeeman): #Convierte una energía de Zeeman en el campo magnético equivalente
    return energia_zeeman / (0.5 * factor_g * magneton_bohr)



###################### MATRICES DE PAULI ######################

matriz_identidad_2 = np.eye(2, dtype=complex)
matriz_sigma_x_2 = np.array([[0.0, 1.0],[1.0, 0.0],],dtype=complex,)
matriz_sigma_y_2 = np.array([[0.0, -1.0j], [1.0j, 0.0]], dtype=complex)
matriz_sigma_z_2 = np.array([[1.0, 0.0],[0.0, -1.0],],dtype=complex,)

matriz_tau_x_2 = np.array([[0.0, 1.0], [1.0, 0.0]], dtype=complex)
matriz_tau_y_2 = np.array([[0.0, -1.0j], [1.0j, 0.0]], dtype=complex)
matriz_tau_z_2 = np.array([[1.0, 0.0], [0.0, -1.0]], dtype=complex)

# En la base (Nambu, espin) usamos productos de Kronecker tau \otimes sigma
sigma_x = np.kron(matriz_identidad_2, matriz_sigma_x_2)
sigma_y = np.kron(matriz_identidad_2, matriz_sigma_y_2)
sigma_z = np.kron(matriz_identidad_2, matriz_sigma_z_2)

tau_x = np.kron(matriz_tau_x_2, matriz_identidad_2)
tau_y = np.kron(matriz_tau_y_2, matriz_identidad_2)
tau_z = np.kron(matriz_tau_z_2, matriz_identidad_2)

sigma_y_tau_z = sigma_y @ tau_z
operador_particula_hueco_local = tau_y @ sigma_y



###################### PERFILES ESPACIALES ######################

def posiciones_sitios(numero_sitios, paso_espacial_usado=None): #Devuelve la posición x_i de cada sitio de la cadena
    if paso_espacial_usado is None:
        paso_espacial_usado = paso_espacial
    posiciones = np.arange(numero_sitios, dtype=float) * paso_espacial_usado

    return posiciones


def convertir_a_perfil_delta(delta, numero_sitios): #Convierte un número complejo o un array en un perfil Delta_i de longitud N
    if np.ndim(delta) == 0:
        perfil_delta = np.full(numero_sitios, complex(delta), dtype=complex)
    else:
        perfil_delta = np.asarray(delta, dtype=complex).copy()

    if len(perfil_delta) != numero_sitios:
        raise ValueError("El perfil Delta_i no tiene el mismo número de sitios que el hilo.")

    return perfil_delta


def perfil_gap_inducido(numero_sitios,delta=delta_inducido,sitio_inicio=sitio_inicio_superconductor,suavidad=suavidad_borde_superconductor,fase=0.0,):#Construye un perfil suave para el gap inducido por el superconductor
    indices = np.arange(numero_sitios, dtype=float)
    #Si sitio_inicio = 0: todo el hilo está cubierto
    if sitio_inicio <= 0:
        peso_superconductor = np.ones(numero_sitios, dtype=float)
    #Si sitio_inicio > 0: queda una pequeña región normal en el extremo izquierdo
    else:
        peso_superconductor = 0.5 * (1.0 + np.tanh((indices - sitio_inicio) / suavidad))
    perfil_delta = delta * peso_superconductor * np.exp(1.0j * fase)
    
    return perfil_delta


def potencial_electrostatico_nulo(numero_sitios): #Devuelve V_i = 0 en todo el hilo
    return np.zeros(numero_sitios, dtype=float)


def potencial_barrera_gaussiana(numero_sitios,altura=altura_barrera,anchura=anchura_barrera,centro=centro_barrera,): #Potencial gaussiano que representa una barrera electrostática suave
    indices = np.arange(numero_sitios, dtype=float)
    return altura * np.exp(-0.5 * ((indices - centro) / anchura)**2)


def potencial_dot_suave(numero_sitios, profundidad=0.0 * meV, longitud_dot=6.0, suavidad=2.0): #Potencial suave en el extremo izquierdo, útil para simular una región tipo dot
    indices = np.arange(numero_sitios, dtype=float)
    return 0.5 * profundidad * (np.tanh((indices - longitud_dot) / suavidad) - 1.0)


def potencial_nanohilo_por_defecto(numero_sitios): #Perfil electrostático usado por defecto en las figuras del nanohilo
    return potencial_barrera_gaussiana(numero_sitios)


###############################################################################
############################## Hamiltoniano BdG ###############################
###############################################################################

def construir_hamiltoniano_bdg(numero_sitios, potencial_quimico, energia_zeeman, perfil_delta, potencial=None, hopping_usado=None, hopping_spin_orbita_usado=None,):


    if potencial is None:
        potencial = potencial_electrostatico_nulo(numero_sitios)
    else:
        potencial = np.asarray(potencial, dtype=float).copy()

    if len(potencial) != numero_sitios:
        raise ValueError("El perfil V_i no tiene el mismo número de sitios que el hilo.")

    perfil_delta = convertir_a_perfil_delta(perfil_delta, numero_sitios)

    if hopping_usado is None:
        hopping_usado = hopping

    if hopping_spin_orbita_usado is None:
        hopping_spin_orbita_usado = hopping_spin_orbita

    hamiltoniano = np.zeros((4 * numero_sitios, 4 * numero_sitios), dtype=complex)
    bloque_hopping = -hopping_usado * tau_z + 1.0j * hopping_spin_orbita_usado * sigma_y_tau_z

    for i in range(numero_sitios):
        bloque_i = slice(4 * i, 4 * (i + 1))
        delta_i = perfil_delta[i]
        bloque_local = ((2.0 * hopping_usado - potencial_quimico + potencial[i]) * tau_z + energia_zeeman * sigma_x + delta_i.real * tau_x - delta_i.imag * tau_y)
        hamiltoniano[bloque_i, bloque_i] = bloque_local
        if i < numero_sitios - 1:
            bloque_j = slice(4 * (i + 1), 4 * (i + 2))
            hamiltoniano[bloque_i, bloque_j] = bloque_hopping
            hamiltoniano[bloque_j, bloque_i] = bloque_hopping.conj().T

    hamiltoniano = 0.5 * (hamiltoniano + hamiltoniano.conj().T)
    return hamiltoniano


def diagonalizar_nanohilo(numero_sitios, potencial_quimico, energia_zeeman, perfil_delta, potencial, hopping_usado, hopping_spin_orbita_usado,):
    #Diagonaliza el Hamiltoniano BdG del nanohilo finito.
    hamiltoniano = construir_hamiltoniano_bdg(numero_sitios,potencial_quimico,energia_zeeman, perfil_delta, potencial=potencial, hopping_usado=hopping_usado, hopping_spin_orbita_usado=hopping_spin_orbita_usado,)
    energias, vectores = la.eigh(hamiltoniano)
    orden = np.argsort(energias)
    energias = energias[orden]
    vectores = vectores[:, orden]
    return energias, vectores, hamiltoniano


def error_hermiticidad(hamiltoniano): #Mide numéricamente si H = H^\dagger
    numerador = la.norm(hamiltoniano - hamiltoniano.conj().T)
    denominador = max(la.norm(hamiltoniano), 10**(-30))
    return float(numerador / denominador)


def error_simetria_particula_hueco(hamiltoniano, numero_sitios): #Comprueba la simetría partícula-hueco BdG: C H^* C^\dagger = -H
    operador_global = np.kron(np.eye(numero_sitios, dtype=complex), operador_particula_hueco_local)
    hamiltoniano_transformado = operador_global @ hamiltoniano.conj() @ operador_global.conj().T
    numerador = la.norm(hamiltoniano_transformado + hamiltoniano)
    denominador = max(la.norm(hamiltoniano), 10**(-30))
    return float(numerador / denominador)


###############################################################################
############################# Observables fisicos #############################
###############################################################################

def separar_componentes_bdg(vector_estado, numero_sitios): #Separa un autovector en las cuatro componentes locales de la base compacta
    estado = np.asarray(vector_estado, dtype=complex).reshape(numero_sitios, 4)
    u_arriba = estado[:, 0]
    u_abajo = estado[:, 1]
    v_abajo = estado[:, 2]
    menos_v_arriba = estado[:, 3]
    return u_arriba, u_abajo, v_abajo, menos_v_arriba


def densidad_estado(vector_estado, numero_sitios): #Calcula la densidad espacial |u|^2 + |v|^2
    u_arriba, u_abajo, v_abajo, menos_v_arriba = separar_componentes_bdg(vector_estado,numero_sitios,)
    densidad = (np.abs(u_arriba)**2 + np.abs(u_abajo)**2 + np.abs(v_abajo)**2 + np.abs(menos_v_arriba)**2)
    return densidad / max(np.sum(densidad), 10**(-30))


def densidad_electron_y_hueco(vector_estado, numero_sitios): #Devuelve por separado la densidad electrónica y la densidad de hueco
    u_arriba, u_abajo, v_abajo, menos_v_arriba = separar_componentes_bdg(vector_estado, numero_sitios,)
    densidad_electron = np.abs(u_arriba)**2 + np.abs(u_abajo)**2
    densidad_hueco = np.abs(v_abajo)**2 + np.abs(menos_v_arriba)**2
    norma = max(np.sum(densidad_electron + densidad_hueco), 10**(-30))
    return densidad_electron / norma, densidad_hueco / norma


def carga_bdg_estado(vector_estado, numero_sitios): #Calcula la carga BdG total: peso electrónico menos peso de hueco
    densidad_electron, densidad_hueco = densidad_electron_y_hueco(vector_estado, numero_sitios)
    carga = np.sum(densidad_electron - densidad_hueco)
    return float(np.real(carga))


def peso_en_bordes(vector_estado, numero_sitios, fraccion_borde=0.20): #Calcula qué fracción de la densidad vive cerca de los dos extremos
    densidad = densidad_estado(vector_estado, numero_sitios)
    numero_borde = max(1, int(fraccion_borde * numero_sitios))
    peso = np.sum(densidad[:numero_borde]) + np.sum(densidad[-numero_borde:])
    return float(peso)


def espin_estado(vector_estado, numero_sitios): #Calcula el valor medio de sigma_x, sigma_y y sigma_z para un estado
    estado = np.asarray(vector_estado, dtype=complex).reshape(numero_sitios, 4)

    espin_x_local = np.zeros(numero_sitios, dtype=float)
    espin_y_local = np.zeros(numero_sitios, dtype=float)
    espin_z_local = np.zeros(numero_sitios, dtype=float)
    for i in range(numero_sitios):
        psi_i = estado[i, :]
        espin_x_local[i] = np.real(np.vdot(psi_i, sigma_x @ psi_i))
        espin_y_local[i] = np.real(np.vdot(psi_i, sigma_y @ psi_i))
        espin_z_local[i] = np.real(np.vdot(psi_i, sigma_z @ psi_i))

    espin_x_total = np.sum(espin_x_local)
    espin_y_total = np.sum(espin_y_local)
    espin_z_total = np.sum(espin_z_local)
    return (float(espin_x_total),float(espin_y_total),float(espin_z_total),espin_x_local,espin_y_local,espin_z_local,)


def longitud_localizacion_desde_densidad(densidad): #Estima la longitud de localización ajustando log(rho) cerca del borde
    densidad = np.asarray(densidad, dtype=float)
    mitad = len(densidad) // 2
    densidad_izquierda = densidad[:mitad]
    posiciones = np.arange(mitad, dtype=float)
    maximo = max(float(np.max(densidad_izquierda)), 10**(-30))
    mascara = densidad_izquierda > maximo * 10**(-3)
    if np.sum(mascara) < 3:
        return np.nan

    pendiente, ordenada = np.polyfit(posiciones[mascara], np.log(densidad_izquierda[mascara]), 1)
    if pendiente >= 0.0:
        return np.nan
    longitud = -2.0 / pendiente

    return float(longitud)


def indice_estado_mas_cercano_a_cero(energias): #Encuentra el índice del autovalor más cercano a E = 0
    return int(np.argmin(np.abs(energias)))


def energia_minima_absoluta(energias): #Devuelve min |E|, útil para detectar estados casi cero
    return float(np.min(np.abs(energias)))


def transformar_particula_hueco(vector_estado, numero_sitios): #Aplica la simetría partícula-hueco al estado BdG
    estado = np.asarray(vector_estado, dtype=complex).reshape(numero_sitios, 4)
    estado_ph = np.zeros_like(estado)
    for i in range(numero_sitios):
        estado_ph[i, :] = operador_particula_hueco_local @ estado[i, :].conj()

    return estado_ph.reshape(4 * numero_sitios)


def alinear_fase(estado_referencia, estado_objetivo): #Elimina la fase global arbitraria entre dos autovectores
    solapamiento = np.vdot(estado_referencia, estado_objetivo)
    if abs(solapamiento) < 10**(-30):
        return estado_objetivo.copy()
    estado_alineado = estado_objetivo * np.exp(-1.0j * np.angle(solapamiento))

    return estado_alineado


def construir_majoranas_desde_pareja(energias, vectores, numero_sitios): #Construye dos combinaciones tipo Majorana desde la pareja pm E más cercana a cero
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


def calcular_ldos(numero_sitios, energias, vectores, energias_ldos, anchura): #Calcula la LDOS(E,i) sumando contribuciones lorentzianas
    ldos = np.zeros((len(energias_ldos), numero_sitios), dtype=float)
    for n, energia_n in enumerate(energias):
        densidad_n = densidad_estado(vectores[:, n], numero_sitios)
        for i, energia in enumerate(energias_ldos):
            ldos[i, :] += densidad_n * lorentziana(energia, energia_n, anchura)

    return ldos


def es_topologico_ideal(potencial_quimico, energia_zeeman, delta_media, alpha_usado=None): #Criterio homogéneo ideal del nanohilo: E_Z > sqrt(mu^2 + Delta^2). Si alpha_R = 0, no lo marcamos como topológico aunque se cierre el gap
    if alpha_usado is None:
        alpha_usado = alpha_rashba
    condicion_zeeman = energia_zeeman > np.sqrt(potencial_quimico**2 + np.abs(delta_media)**2)
    condicion_rashba = abs(alpha_usado) > 10**(-12)
    return bool(condicion_zeeman and condicion_rashba)


def bandas_bulk_nanohilo(potencial_quimico,energia_zeeman,delta_uniforme, numero_k=500, hopping_usado=None, hopping_spin_orbita_usado=None,): #Calcula las bandas bulk de la versión infinita del Hamiltoniano discretizado
    k = np.linspace(-np.pi, np.pi, numero_k)
    bandas = np.zeros((numero_k, 4), dtype=float)
    if hopping_usado is None:
        hopping_usado = hopping

    if hopping_spin_orbita_usado is None:
        hopping_spin_orbita_usado = hopping_spin_orbita

    for i, k_i in enumerate(k):
        xi_k = 2.0 * hopping_usado - potencial_quimico - 2.0 * hopping_usado * np.cos(k_i)
        rashba_k = -2.0 * hopping_spin_orbita_usado * np.sin(k_i)
        hamiltoniano_k = (xi_k * tau_z + rashba_k * sigma_y_tau_z + energia_zeeman * sigma_x + np.real(delta_uniforme) * tau_x - np.imag(delta_uniforme) * tau_y)
        bandas[i, :] = la.eigvalsh(hamiltoniano_k)

    return k, bandas


def gap_bulk_minimo(potencial_quimico, energia_zeeman, delta_uniforme, hopping_usado=None, hopping_spin_orbita_usado=None,): #Calcula el mínimo de las bandas positivas del sistema bulk
    k, bandas = bandas_bulk_nanohilo(potencial_quimico, energia_zeeman, delta_uniforme, numero_k=350, hopping_usado=hopping_usado, hopping_spin_orbita_usado=hopping_spin_orbita_usado,)
    bandas_positivas = bandas[bandas > 10**(-10)]
    if len(bandas_positivas) == 0:
        return 0.0

    return float(np.min(bandas_positivas))


###############################################################################
############################## Barridos fisicos ###############################
###############################################################################

def barrer_zeeman(numero_sitios,potencial_quimico,valores_zeeman,perfil_delta,potencial,): #Barre E_Z y calcula el espectro finito, min|E| y el gap bulk ideal
    todos_los_espectros = []
    energias_minimas = []
    gaps_bulk = []
    delta_media = float(np.mean(np.abs(perfil_delta)))
    for energia_zeeman in valores_zeeman:
        energias, vectores, hamiltoniano = diagonalizar_nanohilo(numero_sitios,potencial_quimico,energia_zeeman,perfil_delta=perfil_delta,potencial=potencial,)
        todos_los_espectros.append(energias)
        energias_minimas.append(energia_minima_absoluta(energias))
        gaps_bulk.append(gap_bulk_minimo(potencial_quimico, energia_zeeman, delta_media))

    return np.array(todos_los_espectros), np.array(energias_minimas), np.array(gaps_bulk)


def barrer_potencial_quimico(numero_sitios,valores_mu,energia_zeeman,perfil_delta=None,potencial=None,): #Barre el potencial químico y calcula el espectro del nanohilo finito
    todos_los_espectros = []
    energias_minimas = []
    gaps_bulk = []
    delta_media = float(np.mean(np.abs(perfil_delta)))
    for mu in valores_mu:
        energias, vectores, hamiltoniano = diagonalizar_nanohilo(numero_sitios, mu, energia_zeeman, perfil_delta=perfil_delta, potencial=potencial,)
        todos_los_espectros.append(energias)
        energias_minimas.append(energia_minima_absoluta(energias))
        gaps_bulk.append(gap_bulk_minimo(mu, energia_zeeman, delta_media))

    return np.array(todos_los_espectros), np.array(energias_minimas), np.array(gaps_bulk)


def barrer_longitud_cadena(valores_longitud, potencial_quimico, energia_zeeman, delta=delta_inducido): #Barre la longitud del nanohilo para estudiar el splitting de Majoranas
    energias_minimas = []
    for N in valores_longitud:
        N = int(N)
        perfil_delta = perfil_gap_inducido(N, delta=delta, sitio_inicio=0)
        potencial = potencial_electrostatico_nulo(N)
        energias, vectores, hamiltoniano = diagonalizar_nanohilo(N,potencial_quimico,energia_zeeman,perfil_delta=perfil_delta,potencial=potencial,)
        energias_minimas.append(energia_minima_absoluta(energias))

    return np.array(energias_minimas)


def calcular_diagrama_fase_mu_zeeman(valores_mu, valores_zeeman, delta=delta_inducido): #Diagrama ideal del nanohilo homogéneo en el plano (mu, E_Z)
    diagrama = np.zeros((len(valores_zeeman), len(valores_mu)))
    for i, energia_zeeman in enumerate(valores_zeeman):
        for j, mu in enumerate(valores_mu):
            diagrama[i, j] = 1.0 if es_topologico_ideal(mu, energia_zeeman, delta) else 0.0

    return diagrama


def barrer_rashba(valores_alpha, numero_sitios_barrido=44): #Barre alpha_R y calcula min|E| para ver cómo ayuda el acoplo espín-órbita
    energias_minimas = []
    pesos_bordes = []
    perfil_delta = perfil_gap_inducido(numero_sitios_barrido, sitio_inicio=0)
    potencial = potencial_electrostatico_nulo(numero_sitios_barrido)
    for alpha in valores_alpha:
        hopping_usado, hopping_so_usado = calcular_hoppings_discretos(alpha, paso_espacial)
        energias, vectores, hamiltoniano = diagonalizar_nanohilo(numero_sitios_barrido, potencial_quimico, zeeman_topologico, perfil_delta=perfil_delta, potencial=potencial, hopping_usado=hopping_usado,hopping_spin_orbita_usado=hopping_so_usado,)
        indice = indice_estado_mas_cercano_a_cero(energias)
        energias_minimas.append(energia_minima_absoluta(energias))
        pesos_bordes.append(peso_en_bordes(vectores[:, indice], numero_sitios_barrido))

    return np.array(energias_minimas), np.array(pesos_bordes)


def calcular_diagrama_alpha_zeeman(valores_alpha, valores_zeeman): # Calcula un mapa de min|E| en el plano (alpha_R, E_Z)
    numero_sitios_mapa = 34
    mapa_minimo = np.zeros((len(valores_zeeman), len(valores_alpha)))
    mapa_ideal = np.zeros((len(valores_zeeman), len(valores_alpha)))
    perfil_delta = perfil_gap_inducido(numero_sitios_mapa, sitio_inicio=0)
    potencial = potencial_electrostatico_nulo(numero_sitios_mapa)
    for i, energia_zeeman in enumerate(valores_zeeman):
        for j, alpha in enumerate(valores_alpha):
            hopping_usado, hopping_so_usado = calcular_hoppings_discretos(alpha, paso_espacial)
            energias, vectores, hamiltoniano = diagonalizar_nanohilo(numero_sitios_mapa,potencial_quimico,energia_zeeman,perfil_delta=perfil_delta,potencial=potencial,hopping_usado=hopping_usado,hopping_spin_orbita_usado=hopping_so_usado,)
            mapa_minimo[i, j] = np.log10(energia_minima_absoluta(energias) / meV + 10**(-8))
            mapa_ideal[i, j] = 1.0 if es_topologico_ideal(potencial_quimico,energia_zeeman,delta_inducido,alpha_usado=alpha,) else 0.0

    return mapa_minimo, mapa_ideal


def barrer_gap_bulk_zeeman(valores_zeeman): #Barre E_Z y compara el gap bulk con la energía mínima del sistema finito
    numero_sitios_barrido = 44
    gaps_bulk = []
    energias_minimas = []
    perfil_delta = perfil_gap_inducido(numero_sitios_barrido, sitio_inicio=0)
    potencial = potencial_electrostatico_nulo(numero_sitios_barrido)
    for energia_zeeman in valores_zeeman:
        gaps_bulk.append(gap_bulk_minimo(potencial_quimico, energia_zeeman, delta_inducido))
        energias, vectores, hamiltoniano = diagonalizar_nanohilo( numero_sitios_barrido, potencial_quimico, energia_zeeman, perfil_delta=perfil_delta, potencial=potencial,)
        energias_minimas.append(energia_minima_absoluta(energias))

    return np.array(gaps_bulk), np.array(energias_minimas)


def barrer_temperatura(valores_temperatura): # Barre temperatura y calcula la LDOS de borde en E = 0 con ensanchamiento térmico
    numero_sitios_barrido = 44
    ldos_borde_cero = []
    anchuras_termicas = []
    perfil_delta = perfil_gap_inducido(numero_sitios_barrido, sitio_inicio=0)
    potencial = potencial_electrostatico_nulo(numero_sitios_barrido)
    energias, vectores, hamiltoniano = diagonalizar_nanohilo(numero_sitios_barrido,potencial_quimico,zeeman_topologico,perfil_delta=perfil_delta,potencial=potencial,)
    for T in valores_temperatura:
        anchura = max(ensanchamiento_ldos, constante_boltzmann * T)
        ldos = calcular_ldos(numero_sitios_barrido, energias, vectores, np.array([0.0]), anchura)
        ldos_borde_cero.append(ldos[0, 0] + ldos[0, -1])
        anchuras_termicas.append(anchura)

    return np.array(ldos_borde_cero), np.array(anchuras_termicas)


def barrer_paso_malla(valores_paso): #Comprueba cómo cambia min|E| al modificar el paso espacial
    energias_minimas = []
    numeros_sitios = []
    longitud_fisica = 900.0 * nm
    for paso in valores_paso:
        N = int(longitud_fisica / paso) + 1
        hopping_usado, hopping_so_usado = calcular_hoppings_discretos(alpha_rashba, paso)
        perfil_delta = perfil_gap_inducido(N, sitio_inicio=0)
        potencial = potencial_electrostatico_nulo(N)
        energias, vectores, hamiltoniano = diagonalizar_nanohilo(N,potencial_quimico,zeeman_topologico, perfil_delta=perfil_delta, potencial=potencial, hopping_usado=hopping_usado, hopping_spin_orbita_usado=hopping_so_usado, )
        energias_minimas.append(energia_minima_absoluta(energias))
        numeros_sitios.append(N)

    return np.array(energias_minimas), np.array(numeros_sitios)


def barrer_barrera_gaussiana(valores_altura, numero_sitios_barrido=44): #Barre la altura de la barrera y mide energía mínima y peso en bordes
    energias_minimas = []
    pesos_bordes = []
    densidades = []
    perfil_delta = perfil_gap_inducido(numero_sitios_barrido)
    for altura in valores_altura:
        potencial = potencial_barrera_gaussiana(numero_sitios_barrido, altura=altura)
        energias, vectores, hamiltoniano = diagonalizar_nanohilo(numero_sitios_barrido,potencial_quimico,zeeman_topologico,perfil_delta=perfil_delta,potencial=potencial,)
        indice = indice_estado_mas_cercano_a_cero(energias)
        densidad = densidad_estado(vectores[:, indice], numero_sitios_barrido)
        energias_minimas.append(energia_minima_absoluta(energias))
        pesos_bordes.append(peso_en_bordes(vectores[:, indice], numero_sitios_barrido))
        densidades.append(densidad)

    return np.array(energias_minimas), np.array(pesos_bordes), densidades


def calcular_mapa_mu_zeeman_minimo(valores_mu, valores_zeeman): #Calcula un mapa finito de log10(min|E|) en el plano (mu, E_Z)
    numero_sitios_mapa = 34
    mapa_log = np.zeros((len(valores_zeeman), len(valores_mu)))
    perfil_delta = perfil_gap_inducido(numero_sitios_mapa, sitio_inicio=0)
    potencial = potencial_electrostatico_nulo(numero_sitios_mapa)
    for i, energia_zeeman in enumerate(valores_zeeman):
        for j, mu in enumerate(valores_mu):
            energias, vectores, hamiltoniano = diagonalizar_nanohilo(numero_sitios_mapa,mu,energia_zeeman,perfil_delta=perfil_delta,potencial=potencial,)
            mapa_log[i, j] = np.log10(energia_minima_absoluta(energias) / meV + 10**(-8))

    return mapa_log


def barrer_diagnosticos_modo_zeeman(valores_zeeman, numero_sitios_barrido=44): #Barre E_Z y mide peso de borde, carga, espín y longitud de localización
    pesos_bordes = []
    cargas_bdg = []
    espines_x = []
    espines_y = []
    espines_z = []
    longitudes = []
    perfil_delta = perfil_gap_inducido(numero_sitios_barrido, sitio_inicio=0)
    potencial = potencial_electrostatico_nulo(numero_sitios_barrido)
    for energia_zeeman in valores_zeeman:
        energias, vectores, hamiltoniano = diagonalizar_nanohilo(numero_sitios_barrido,potencial_quimico,energia_zeeman,perfil_delta=perfil_delta,potencial=potencial,)
        indice = indice_estado_mas_cercano_a_cero(energias)
        estado = vectores[:, indice]
        densidad = densidad_estado(estado, numero_sitios_barrido)
        spin_x, spin_y, spin_z, spin_x_local, spin_y_local, spin_z_local = espin_estado(estado,numero_sitios_barrido,)
        pesos_bordes.append(peso_en_bordes(estado, numero_sitios_barrido))
        cargas_bdg.append(carga_bdg_estado(estado, numero_sitios_barrido))
        espines_x.append(spin_x)
        espines_y.append(spin_y)
        espines_z.append(spin_z)
        longitudes.append(longitud_localizacion_desde_densidad(densidad))

    return (np.array(pesos_bordes),np.array(cargas_bdg),np.array(espines_x),np.array(espines_y),np.array(espines_z),np.array(longitudes),)


def barrer_gap_homogeneo_vs_perfil_inducido(valores_zeeman, numero_sitios_barrido=44): #Compara el nanohilo cubierto completamente con el perfil inducido suave
    energias_homogeneas = []
    energias_perfil = []
    perfil_delta_homogeneo = perfil_gap_inducido(numero_sitios_barrido, sitio_inicio=0)
    perfil_delta_suave = perfil_gap_inducido(numero_sitios_barrido)
    potencial_nulo = potencial_electrostatico_nulo(numero_sitios_barrido)
    potencial_suave = potencial_nanohilo_por_defecto(numero_sitios_barrido)
    for energia_zeeman in valores_zeeman:
        energias_h, vectores_h, hamiltoniano_h = diagonalizar_nanohilo(numero_sitios_barrido,potencial_quimico,energia_zeeman,perfil_delta=perfil_delta_homogeneo,potencial=potencial_nulo,)
        energias_p, vectores_p, hamiltoniano_p = diagonalizar_nanohilo(numero_sitios_barrido,potencial_quimico,energia_zeeman,perfil_delta=perfil_delta_suave,potencial=potencial_suave,)
        energias_homogeneas.append(energia_minima_absoluta(energias_h))
        energias_perfil.append(energia_minima_absoluta(energias_p))

    return np.array(energias_homogeneas), np.array(energias_perfil)


###############################################################################
################################## Graficar ###################################
###############################################################################

def preparar_estilo_graficas(): #Estilo sencillo para que las figuras sean claras en un informe
    plt.style.use("default")
    plt.rcParams.update({"figure.figsize": (9, 5.5),"axes.grid": True,"grid.alpha": 0.25,"font.size": 11,"axes.titlesize": 13,"axes.labelsize": 12,"legend.fontsize": 10,})


def graficar_bandas_bulk_tres_regimenes(): # Grafica las bandas bulk en régimen trivial, transición y topológico con diferente campo Zeeman
    casos = [("trivial", zeeman_trivial),("transición", zeeman_transicion),("topológico", zeeman_topologico),]

    figura, ejes = plt.subplots(1, 3, figsize=(13, 4), sharey=True)
    for eje, (nombre, energia_zeeman) in zip(ejes, casos):
        k, bandas = bandas_bulk_nanohilo(potencial_quimico,energia_zeeman,delta_inducido,)
        for indice_banda in range(4):
            eje.plot(k / np.pi, bandas[:, indice_banda], linewidth=1.8)
        eje.axhline(0.0, color="black", linewidth=1.0, alpha=0.5)
        eje.set_title(f"{nombre}: $E_Z={energia_zeeman / meV:.2f}$ meV")
        eje.set_xlabel(r"$ka/\pi$")
    ejes[0].set_ylabel("Energía (meV)")
    figura.suptitle("Bandas bulk del nanohilo BdG")
    figura.tight_layout()
    return figura



def graficar_espectro_finito_vs_zeeman(datos_barrido=None): #Grafica el espectro finito frente a E_Z.
    valores_zeeman = np.linspace(0.0, 0.75 * meV, numero_puntos_barrido)
    espectros, energias_minimas, gaps_medios = barrer_zeeman(numero_sitios,potencial_quimico,valores_zeeman,)

    figura, eje = plt.subplots()
    for i, energia_zeeman in enumerate(valores_zeeman):
        energias_cerca_cero = espectros[i][np.abs(espectros[i]) < 1.2 * meV]
        eje.scatter(np.full_like(energias_cerca_cero, energia_zeeman), energias_cerca_cero, s=9, color="black", alpha=0.65,)

    delta_referencia = max(gaps_medios[0], 10**(-12))                       # Marcamos una guia orientativa usando el gap medio de E_Z = 0
    zeeman_guia = np.sqrt(potencial_quimico**2 + delta_referencia**2)
    eje.axvline(zeeman_guia, linestyle="--", linewidth=2.0, label=r"$E_Z\simeq\sqrt{\mu^2+\Delta^2}$")
    eje.axhline(0.0, color="black", linewidth=1.0, alpha=0.5)
    eje.set_title("Espectro finito autoconsistente frente al campo de Zeeman")
    eje.set_xlabel(r"Energia de Zeeman $E_Z$ (meV)")
    eje.set_ylabel("Energia (meV)")
    eje.legend()
    figura.tight_layout()

    return figura



def graficar_energia_minima_vs_zeeman(datos_barrido=None): #Grafica min|E| frente a E_Z 
    valores_zeeman = np.linspace(0.0, 0.75 * meV, numero_puntos_barrido)
    espectros, energias_minimas, gaps_medios = barrer_zeeman(numero_sitios,potencial_quimico,valores_zeeman,)

    figura, eje = plt.subplots()
    eje.plot(valores_zeeman / meV, energias_minimas / meV, linewidth=2.2, label="sistema finito")
    eje.axvline(zeeman_transicion / meV, linestyle="--", linewidth=2.0, label="criterio ideal")
    eje.set_title(r"Energía más cercana a cero: $\min |E|$")
    eje.set_xlabel(r"Energía de Zeeman $E_Z$ (meV)")
    eje.set_ylabel(r"$\min |E|$ (meV)")
    eje.legend()

    figura.tight_layout()

    return figura

def graficar_espectro_finito_vs_mu(datos_barrido=None): #Grafica el espectro finito frente a mu, con E_Z fijo en el regimen topologico
    valores_mu = np.linspace(-0.8 * meV, 0.8 * meV, numero_puntos_barrido)
    espectros, energias_minimas, gaps_medios = barrer_potencial_quimico(numero_sitios,valores_mu,zeeman_topologico,)

    figura, eje = plt.subplots()
    for i, mu in enumerate(valores_mu):
        energias_cerca_cero = espectros[i][np.abs(espectros[i]) < 1.2 * meV]
        eje.scatter(np.full_like(energias_cerca_cero, mu),energias_cerca_cero,s=9,color="black",alpha=0.65,)

    mu_critico = np.sqrt(max(zeeman_topologico**2 - delta_inducido**2, 0.0))
    eje.axvline(-mu_critico / meV, linestyle="--", linewidth=2.0, label="frontera ideal")
    eje.axvline(+mu_critico / meV, linestyle="--", linewidth=2.0)
    eje.axhline(0.0, color="black", linewidth=1.0, alpha=0.5)

    eje.set_title("Espectro finito frente al potencial químico")
    eje.set_xlabel(r"Potencial químico $\mu$ (meV)")
    eje.set_ylabel("Energía BdG (meV)")
    eje.legend()

    figura.tight_layout()

    return figura


def graficar_energia_minima_vs_mu(datos_barrido=None): #Grafica min|E| frente al potencial químico.
    valores_mu = np.linspace(-0.8 * meV, 0.8 * meV, numero_puntos_barrido)
    espectros, energias_minimas, gaps_bulk = barrer_potencial_quimico(numero_sitios,valores_mu,zeeman_topologico,)
    
    figura, eje = plt.subplots()
    mu_critico = np.sqrt(max(zeeman_topologico**2 - delta_inducido**2, 0.0))
    eje.plot(valores_mu / meV, energias_minimas / meV, linewidth=2.2)
    eje.axvline(-mu_critico / meV, linestyle="--", linewidth=2.0, label="frontera ideal")
    eje.axvline(+mu_critico / meV, linestyle="--", linewidth=2.0)
    eje.set_title(r"Energía más cercana a cero frente a $\mu$")
    eje.set_xlabel(r"Potencial químico $\mu$ (meV)")
    eje.set_ylabel(r"$\min |E|$ (meV)")
    eje.legend()
    figura.tight_layout()

    return figura



def graficar_ldos_tres_regimenes(): #Grafica LDOS(E,j) para regimen trivial, intermedio y topologico
    casos = [("trivial", zeeman_trivial),("transición", zeeman_transicion),("topologico", zeeman_topologico),]
    energias_ldos = np.linspace(-1.0 * meV, 1.0 * meV, numero_energias)
    figura = plt.figure(figsize=(14.5, 4)) #Creamos una figura algo más ancha para dejar una columna propia a la barra de color
    rejilla = figura.add_gridspec(1, 4, width_ratios=[1.0, 1.0, 1.0, 0.045], wspace=0.12) #Reservamos tres columnas para LDOS y una columna estrecha para la barra
    ejes = [figura.add_subplot(rejilla[0, i]) for i in range(3)] #Creamos los tres ejes principales de la LDOS
    ejes[1].sharey(ejes[0]) #Compartimos eje y para comparar los tres regímenes
    ejes[2].sharey(ejes[0]) #Compartimos eje y para comparar los tres regímenes
    eje_colorbar = figura.add_subplot(rejilla[0, 3]) #Creamos un eje independiente para la barra de color

    perfil_delta = perfil_gap_inducido(numero_sitios)
    potencial = potencial_nanohilo_por_defecto(numero_sitios)

    for eje, (nombre, energia_zeeman) in zip(ejes, casos):
        energias, vectores, hamiltoniano = diagonalizar_nanohilo(numero_sitios,potencial_quimico,energia_zeeman,perfil_delta=perfil_delta,potencial=potencial,)
        ldos = calcular_ldos(numero_sitios, energias, vectores, energias_ldos, ensanchamiento_ldos)
        imagen = eje.imshow(ldos,aspect="auto",origin="lower",extent=[0, numero_sitios - 1, energias_ldos[0] / meV, energias_ldos[-1] / meV],cmap="magma",)

        eje.axhline(0.0, color="white", linewidth=1.0, alpha=0.7)
        eje.set_title(f"{nombre}: $E_Z={energia_zeeman / meV:.2f}$ meV")
        eje.set_xlabel("Sitio i")

    ejes[0].set_ylabel("Energía (meV)")
    figura.suptitle(r"LDOS$(E,i)$ del nanohilo")
    barra = figura.colorbar(imagen, cax=eje_colorbar) #Colocamos la barra de color en su propia columna para que no invada ningún panel
    barra.set_label("LDOS") #Etiqueta de la barra de color
    figura.suptitle(r"LDOS$(E,i)$ del nanohilo autoconsistente", y=1.03) #Subimos un poco el título para que no choque con los paneles

    return figura


def graficar_modo_cero_topologico_y_trivial(): #Compara el estado más cercano a cero en un caso trivial y uno topológico
    casos = [("trivial", zeeman_trivial),("topológico", zeeman_topologico),]

    figura, ejes = plt.subplots(1, 2, figsize=(12, 4), sharey=True)
    for eje, (nombre, energia_zeeman) in zip(ejes, casos):
        energias, vectores, hamiltoniano = diagonalizar_nanohilo(numero_sitios,potencial_quimico,energia_zeeman,)
        indice = indice_estado_mas_cercano_a_cero(energias)
        densidad = densidad_estado(vectores[:, indice], numero_sitios)
        eje.plot(np.arange(numero_sitios), densidad, linewidth=2.2)
        eje.set_title(f"{nombre}: E = {energias[indice] / meV:.3e} meV")
        eje.set_xlabel("Sitio i")

    ejes[0].set_ylabel("Densidad normalizada")
    figura.suptitle("Estado BdG más cercano a energía cero")
    figura.tight_layout()

    return figura


def graficar_majoranas_separadas(): #Grafica las dos combinaciones tipo Majorana del modo casi cero
    energias, vectores, hamiltoniano = diagonalizar_nanohilo(numero_sitios,potencial_quimico,zeeman_topologico,)
    densidad_izquierda, densidad_derecha = construir_majoranas_desde_pareja(energias,vectores,numero_sitios,)

    figura, eje = plt.subplots()
    eje.plot(densidad_izquierda, linewidth=2.4, label="Majorana izquierda")
    eje.plot(densidad_derecha, linewidth=2.4, label="Majorana derecha")
    eje.set_title("Descomposición del modo casi cero en dos Majoranas")
    eje.set_xlabel("Sitio i")
    eje.set_ylabel("Densidad normalizada")
    eje.legend()
    figura.tight_layout()

    return figura


def graficar_particula_hueco_modo_cero(): #Grafica contenido electrónico y de hueco del estado más cercano a cero
    energias, vectores, hamiltoniano = diagonalizar_nanohilo(numero_sitios, potencial_quimico, zeeman_topologico,)
    indice = indice_estado_mas_cercano_a_cero(energias)
    densidad_electron, densidad_hueco = densidad_electron_y_hueco(vectores[:, indice], numero_sitios)

    figura, eje = plt.subplots()
    eje.plot(densidad_electron, linewidth=2.2, label="parte electrónica")
    eje.plot(densidad_hueco, linewidth=2.2, linestyle="--", label="parte de hueco")
    eje.set_title("Contenido partícula-hueco del modo cercano a cero")
    eje.set_xlabel("Sitio i")
    eje.set_ylabel("Densidad normalizada")
    eje.legend()

    figura.tight_layout()

    return figura


def graficar_perfil_gap_y_potencial(): #Grafica el perfil de gap inducido y la barrera electrostática
    perfil_delta = perfil_gap_inducido(numero_sitios)
    potencial = potencial_nanohilo_por_defecto(numero_sitios)
    posiciones = posiciones_sitios(numero_sitios) / nm

    figura, eje_1 = plt.subplots()
    eje_2 = eje_1.twinx()

    eje_1.plot(posiciones, np.abs(perfil_delta) / meV, linewidth=2.4, color="tab:blue", label=r"$|\Delta_i|$")
    eje_2.plot(posiciones, potencial / meV, linewidth=2.4, color="tab:red", label=r"$V_i$")

    eje_1.set_title("Perfiles espaciales del nanohilo")
    eje_1.set_xlabel("Posición x (nm)")
    eje_1.set_ylabel(r"Gap inducido $|\Delta_i|$ (meV)", color="tab:blue")
    eje_2.set_ylabel(r"Potencial $V_i$ (meV)", color="tab:red")

    lineas = eje_1.get_lines() + eje_2.get_lines()
    etiquetas = [linea.get_label() for linea in lineas]
    eje_1.legend(lineas, etiquetas, loc="upper right")
    
    figura.tight_layout()

    return figura


def graficar_splitting_vs_longitud(): #Grafica el splitting del modo de borde frente a la longitud del hilo
    valores_longitud = np.arange(28, 101, 6)
    energias_minimas = barrer_longitud_cadena(valores_longitud,potencial_quimico,zeeman_topologico,)

    figura, eje = plt.subplots()
    eje.semilogy(valores_longitud, energias_minimas / meV, marker="o", linewidth=2.2)
    eje.set_title("Splitting de Majoranas frente a la longitud")
    eje.set_xlabel("Número de sitios N")
    eje.set_ylabel(r"$\min |E|$ (meV)")
    figura.tight_layout()

    return figura


def graficar_comparacion_gap_localizacion(): #Compara cómo cambia la localización de la densidad al variar el gap inducido
    casos = [("gap débil", delta_debil),("gap medio", delta_medio),("gap fuerte", delta_fuerte),]

    figura, eje = plt.subplots()
    for nombre, delta in casos:
        perfil_delta = perfil_gap_inducido(numero_sitios, delta=delta, sitio_inicio=0)
        potencial = potencial_electrostatico_nulo(numero_sitios)
        energias, vectores, hamiltoniano = diagonalizar_nanohilo(numero_sitios,potencial_quimico,zeeman_topologico,perfil_delta=perfil_delta,potencial=potencial,)
        indice = indice_estado_mas_cercano_a_cero(energias)
        densidad = densidad_estado(vectores[:, indice], numero_sitios)
        eje.plot(densidad, linewidth=2.0, label=f"{nombre}: Delta={delta / meV:.2f} meV")
    eje.set_title("Efecto del gap inducido en la localización")
    eje.set_xlabel("Sitio i")
    eje.set_ylabel("Densidad normalizada")
    eje.legend()
    figura.tight_layout()

    return figura


def graficar_diagrama_fase_mu_zeeman(): #Grafica el diagrama ideal de fase en el plano (mu, E_Z)
    valores_mu = np.linspace(-0.8 * meV, 0.8 * meV, 220)
    valores_zeeman = np.linspace(0.0, 0.9 * meV, 220)
    diagrama = calcular_diagrama_fase_mu_zeeman(valores_mu, valores_zeeman)
    
    figura, eje = plt.subplots()
    imagen = eje.imshow(diagrama,aspect="auto",origin="lower",extent=[valores_mu[0] / meV, valores_mu[-1] / meV, valores_zeeman[0] / meV, valores_zeeman[-1] / meV],cmap="viridis",)
    eje.set_title(r"Diagrama ideal del nanohilo: $(\mu,E_Z)$")
    eje.set_xlabel(r"Potencial químico $\mu$ (meV)")
    eje.set_ylabel(r"Energía de Zeeman $E_Z$ (meV)")
    barra = figura.colorbar(imagen, ax=eje)
    barra.set_label("0 = trivial, 1 = topológico")
    figura.tight_layout()

    return figura


def graficar_barrido_rashba(): #Grafica min|E| y peso en bordes al variar alpha_R
    valores_alpha = np.linspace(0.0, 160.0 * meV * nm, 27)
    energias_minimas, pesos_bordes = barrer_rashba(valores_alpha)

    figura, eje_1 = plt.subplots()
    eje_2 = eje_1.twinx()
    eje_1.plot(valores_alpha / (meV * nm), energias_minimas / meV, linewidth=2.2, color="tab:blue", label=r"$\min |E|$")
    eje_2.plot(valores_alpha / (meV * nm), pesos_bordes, linewidth=2.2, color="tab:orange", label="peso en bordes")
    eje_1.set_title("Barrido del acoplo Rashba")
    eje_1.set_xlabel(r"$\alpha_R$ (meV nm)")
    eje_1.set_ylabel(r"$\min |E|$ (meV)", color="tab:blue")
    eje_2.set_ylabel("Peso en bordes", color="tab:orange")
    lineas = eje_1.get_lines() + eje_2.get_lines()
    etiquetas = [linea.get_label() for linea in lineas]
    eje_1.legend(lineas, etiquetas, loc="upper right")
    figura.tight_layout()

    return figura


def graficar_diagrama_alpha_zeeman(): #Grafica un mapa en el plano (alpha_R, E_Z)
    valores_alpha = np.linspace(0.0, 160.0 * meV * nm, 21)
    valores_zeeman = np.linspace(0.0, 0.85 * meV, 21)
    mapa_minimo, mapa_ideal = calcular_diagrama_alpha_zeeman(valores_alpha, valores_zeeman)

    figura, ejes = plt.subplots(1, 2, figsize=(12, 4), sharey=True)
    imagen_1 = ejes[0].imshow(mapa_minimo,aspect="auto",origin="lower",extent=[valores_alpha[0] / (meV * nm), valores_alpha[-1] / (meV * nm), valores_zeeman[0] / meV, valores_zeeman[-1] / meV], cmap="magma",)
    ejes[0].set_title(r"$\log_{10}(\min |E|)$")
    ejes[0].set_xlabel(r"$\alpha_R$ (meV nm)")
    ejes[0].set_ylabel(r"$E_Z$ (meV)")
    figura.colorbar(imagen_1, ax=ejes[0], label=r"$\log_{10}(E/\mathrm{meV})$")
    imagen_2 = ejes[1].imshow(mapa_ideal,aspect="auto", origin="lower", extent=[valores_alpha[0] / (meV * nm), valores_alpha[-1] / (meV * nm), valores_zeeman[0] / meV, valores_zeeman[-1] / meV], cmap="viridis",)
    ejes[1].set_title("Criterio ideal")
    ejes[1].set_xlabel(r"$\alpha_R$ (meV nm)")
    figura.colorbar(imagen_2, ax=ejes[1], label="0 = trivial, 1 = topológico")
    figura.suptitle(r"Diagrama en el plano $(\alpha_R,E_Z)$")
    figura.tight_layout()

    return figura


def graficar_gap_bulk_vs_zeeman(datos_bulk=None): #Grafica el gap bulk mínimo frente a E_Z
    valores_zeeman = np.linspace(0.0, 0.75 * meV, 25)
    gaps_bulk, energias_minimas = barrer_gap_bulk_zeeman(valores_zeeman)

    figura, eje = plt.subplots()
    eje.plot(valores_zeeman / meV, gaps_bulk / meV, linewidth=2.2)
    eje.axvline(zeeman_transicion / meV, linestyle="--", linewidth=2.0, label="criterio ideal")
    eje.set_title("Gap bulk del nanohilo homogéneo")
    eje.set_xlabel(r"Energía de Zeeman $E_Z$ (meV)")
    eje.set_ylabel("Gap bulk mínimo (meV)")
    eje.legend()
    figura.tight_layout()

    return figura


def graficar_comparacion_gap_bulk_borde(datos_bulk=None): #Compara el gap bulk con la energía más baja de un hilo finito
    valores_zeeman = np.linspace(0.0, 0.75 * meV, 25)
    gaps_bulk, energias_minimas = barrer_gap_bulk_zeeman(valores_zeeman)
    
    figura, eje = plt.subplots()
    eje.plot(valores_zeeman / meV, gaps_bulk / meV, linewidth=2.2, label="gap bulk")
    eje.plot(valores_zeeman / meV, energias_minimas / meV, linewidth=2.2, label="min|E| finito")
    eje.axvline(zeeman_transicion / meV, linestyle="--", linewidth=2.0, label="criterio ideal")
    eje.set_title("Comparación bulk-borde en el nanohilo")
    eje.set_xlabel(r"Energía de Zeeman $E_Z$ (meV)")
    eje.set_ylabel("Energía (meV)")
    eje.legend()
    figura.tight_layout()

    return figura


def graficar_perfil_gap_vs_region_superconductora(): #Compara perfiles Delta_i al mover el inicio de la cobertura superconductora
    sitios_inicio = [0, 4, 10, 16]
    posiciones = posiciones_sitios(numero_sitios) / nm

    figura, eje = plt.subplots()
    for inicio in sitios_inicio:
        perfil_delta = perfil_gap_inducido(numero_sitios, sitio_inicio=inicio)
        eje.plot(posiciones, np.abs(perfil_delta) / meV, linewidth=2.0, label=f"inicio SC = {inicio}")
    eje.set_title("Perfil del gap inducido para distintas coberturas")
    eje.set_xlabel("Posición x (nm)")
    eje.set_ylabel(r"$|\Delta_i|$ (meV)")
    eje.legend()
    figura.tight_layout()

    return figura


def graficar_barrido_temperatura(): #Grafica cómo la temperatura ensancha la señal LDOS de borde en E = 0
    valores_temperatura = np.linspace(0.01 * kelvin, 0.40 * kelvin, 24)
    ldos_borde_cero, anchuras_termicas = barrer_temperatura(valores_temperatura)

    figura, eje_1 = plt.subplots()
    eje_2 = eje_1.twinx()
    eje_1.plot(valores_temperatura / kelvin, ldos_borde_cero, linewidth=2.2, color="tab:blue", label="LDOS borde en E=0")
    eje_2.plot(valores_temperatura / kelvin, anchuras_termicas / meV, linewidth=2.2, color="tab:red", label="anchura térmica")
    eje_1.set_title("Efecto de la temperatura en la LDOS de borde")
    eje_1.set_xlabel("Temperatura T (K)")
    eje_1.set_ylabel("LDOS de borde", color="tab:blue")
    eje_2.set_ylabel("Anchura efectiva (meV)", color="tab:red")
    lineas = eje_1.get_lines() + eje_2.get_lines()
    etiquetas = [linea.get_label() for linea in lineas]
    eje_1.legend(lineas, etiquetas, loc="upper right")
    figura.tight_layout()

    return figura


def graficar_convergencia_paso_malla(): #Grafica una comprobación sencilla de convergencia con el paso espacial
    valores_paso = np.array([30.0, 25.0, 20.0, 15.0]) * nm
    energias_minimas, numeros_sitios = barrer_paso_malla(valores_paso)

    figura, eje = plt.subplots()
    eje.plot(valores_paso / nm, energias_minimas / meV, marker="o", linewidth=2.2)
    for paso, energia, N in zip(valores_paso, energias_minimas, numeros_sitios):
        eje.annotate(f"N={N}", (paso / nm, energia / meV), textcoords="offset points", xytext=(6, 6))
    eje.invert_xaxis()
    eje.set_title("Convergencia al refinar la malla espacial")
    eje.set_xlabel("Paso espacial a (nm)")
    eje.set_ylabel(r"$\min |E|$ (meV)")
    figura.tight_layout()

    return figura


def graficar_barrido_barrera_gaussiana(): #Grafica el efecto de la altura de la barrera suave
    valores_altura = np.linspace(0.0, 0.6 * meV, 24)
    energias_minimas, pesos_bordes, densidades = barrer_barrera_gaussiana(valores_altura)

    figura, ejes = plt.subplots(1, 2, figsize=(12, 4))
    ejes[0].plot(valores_altura / meV, energias_minimas / meV, linewidth=2.2)
    ejes[0].set_title("Energía mínima")
    ejes[0].set_xlabel("Altura de barrera (meV)")
    ejes[0].set_ylabel(r"$\min |E|$ (meV)")
    ejes[1].plot(valores_altura / meV, pesos_bordes, linewidth=2.2, color="tab:orange")
    ejes[1].set_title("Peso en los bordes")
    ejes[1].set_xlabel("Altura de barrera (meV)")
    ejes[1].set_ylabel("Peso de borde")
    figura.tight_layout()

    return figura


def graficar_mapa_mu_zeeman_minimo(): #Grafica un mapa finito de log10(min|E|) en el plano (mu, E_Z)
    valores_mu = np.linspace(-0.7 * meV, 0.7 * meV, 25)
    valores_zeeman = np.linspace(0.0, 0.8 * meV, 25)
    mapa_log = calcular_mapa_mu_zeeman_minimo(valores_mu, valores_zeeman)

    figura, eje = plt.subplots()
    imagen = eje.imshow(mapa_log,aspect="auto",origin="lower",extent=[valores_mu[0] / meV, valores_mu[-1] / meV, valores_zeeman[0] / meV, valores_zeeman[-1] / meV],cmap="magma",)
    eje.set_title(r"Mapa finito de $\log_{10}(\min |E|)$")
    eje.set_xlabel(r"Potencial químico $\mu$ (meV)")
    eje.set_ylabel(r"Energía de Zeeman $E_Z$ (meV)")
    barra = figura.colorbar(imagen, ax=eje)
    barra.set_label(r"$\log_{10}(E/\mathrm{meV})$")
    figura.tight_layout()

    return figura


def graficar_peso_bordes_vs_zeeman(datos_diagnosticos=None): #Grafica el peso de borde del modo más cercano a cero
    valores_zeeman = np.linspace(0.0, 0.75 * meV, 25)
    datos_diagnosticos = barrer_diagnosticos_modo_zeeman(valores_zeeman)
    pesos_bordes, cargas_bdg, espines_x, espines_y, espines_z, longitudes = datos_diagnosticos

    figura, eje = plt.subplots()
    eje.plot(valores_zeeman / meV, pesos_bordes, linewidth=2.2)
    eje.axvline(zeeman_transicion / meV, linestyle="--", linewidth=2.0, label="criterio ideal")
    eje.set_title("Peso de borde del modo cercano a cero")
    eje.set_xlabel(r"Energía de Zeeman $E_Z$ (meV)")
    eje.set_ylabel("Peso en bordes")
    eje.legend()
    figura.tight_layout()

    return figura



def graficar_espin_modo_vs_zeeman(datos_diagnosticos=None): #Grafica el espín total del modo cercano a cero frente a E_Z
    valores_zeeman = np.linspace(0.0, 0.75 * meV, 25)
    datos_diagnosticos = barrer_diagnosticos_modo_zeeman(valores_zeeman)
    pesos_bordes, cargas_bdg, espines_x, espines_y, espines_z, longitudes = datos_diagnosticos

    figura, eje = plt.subplots()
    eje.plot(valores_zeeman / meV, espines_x, linewidth=2.0, label=r"$\langle\sigma_x\rangle$")
    eje.plot(valores_zeeman / meV, espines_y, linewidth=2.0, label=r"$\langle\sigma_y\rangle$")
    eje.plot(valores_zeeman / meV, espines_z, linewidth=2.0, label=r"$\langle\sigma_z\rangle$")
    eje.axvline(zeeman_transicion / meV, linestyle="--", linewidth=2.0, label="criterio ideal")
    eje.set_title("Espín del modo cercano a cero")
    eje.set_xlabel(r"Energía de Zeeman $E_Z$ (meV)")
    eje.set_ylabel("Valor medio de espín")
    eje.legend()
    figura.tight_layout()

    return figura


def graficar_longitud_localizacion_vs_zeeman(datos_diagnosticos=None): #Grafica una estimación de la longitud de localización del modo
    valores_zeeman = np.linspace(0.0, 0.75 * meV, 25)
    datos_diagnosticos = barrer_diagnosticos_modo_zeeman(valores_zeeman)
    pesos_bordes, cargas_bdg, espines_x, espines_y, espines_z, longitudes = datos_diagnosticos

    figura, eje = plt.subplots()
    eje.plot(valores_zeeman / meV, longitudes * paso_espacial / nm, linewidth=2.2)
    eje.axvline(zeeman_transicion / meV, linestyle="--", linewidth=2.0, label="criterio ideal")
    eje.set_title("Longitud de localización estimada")
    eje.set_xlabel(r"Energía de Zeeman $E_Z$ (meV)")
    eje.set_ylabel("Longitud de localización (nm)")
    eje.legend()
    figura.tight_layout()

    return figura


def graficar_comparacion_gap_homogeneo_inducido(): #Compara gap homogéneo con perfil inducido suave e inhomogéneo
    valores_zeeman = np.linspace(0.0, 0.75 * meV, 25)
    energias_homogeneas, energias_perfil = barrer_gap_homogeneo_vs_perfil_inducido(valores_zeeman)

    figura, eje = plt.subplots()
    eje.plot(valores_zeeman / meV, energias_homogeneas / meV, linewidth=2.2, label="gap homogéneo")
    eje.plot(valores_zeeman / meV, energias_perfil / meV, linewidth=2.2, label="perfil inducido + barrera")
    eje.axvline(zeeman_transicion / meV, linestyle="--", linewidth=2.0, label="criterio ideal")
    eje.set_title("Efecto de la inhomogeneidad del nanohilo")
    eje.set_xlabel(r"Energía de Zeeman $E_Z$ (meV)")
    eje.set_ylabel(r"$\min |E|$ (meV)")
    eje.legend()
    figura.tight_layout()

    return figura


def bandas_bulk_desde_bloques(potencial_quimico, energia_zeeman, delta_uniforme, numero_k=80, hopping_usado=None, hopping_spin_orbita_usado=None,):
    k = np.linspace(-np.pi, np.pi, numero_k)
    bandas = np.zeros((numero_k, 4), dtype=float)

    bloque_local_bulk = ( (2.0 * hopping_usado - potencial_quimico) * tau_z + energia_zeeman * sigma_x + np.real(delta_uniforme) * tau_x - np.imag(delta_uniforme) * tau_y)
    bloque_hopping = ( -hopping_usado * tau_z + 1.0j * hopping_spin_orbita_usado * sigma_y_tau_z)

    for i, k_i in enumerate(k):
        hamiltoniano_k = ( bloque_local_bulk + np.exp(1.0j * k_i) * bloque_hopping + np.exp(-1.0j * k_i) * bloque_hopping.conj().T )
        bandas[i, :] = la.eigvalsh(hamiltoniano_k)

    return k, bandas

def graficar_bandas_bulk_tres_regimenes():
    casos = [("trivial", zeeman_trivial), ("transición", zeeman_transicion), ("topológico", zeeman_topologico),]

    figura, ejes = plt.subplots(1, 3, figsize=(13, 4), sharey=True)
    for eje, (nombre, energia_zeeman) in zip(ejes, casos):
        # Bandas usando la expresión cerrada en k.
        k, bandas_teoricas = bandas_bulk_nanohilo(potencial_quimico, energia_zeeman, delta_inducido, numero_k=500,)
        # Bandas reconstruidas directamente desde los bloques h_i y T.
        k_bloques, bandas_bloques = bandas_bulk_desde_bloques(potencial_quimico, energia_zeeman, delta_inducido, numero_k=70,)

        for indice_banda in range(4):
            eje.plot(k / np.pi, bandas_teoricas[:, indice_banda] / meV, linewidth=1.8,)
            eje.scatter(k_bloques / np.pi, bandas_bloques[:, indice_banda] / meV, s=12, facecolors="none", edgecolors="black", alpha=0.75,)

        eje.axhline(0.0, color="black", linewidth=1.0, alpha=0.5)
        eje.set_title(f"{nombre}: $E_Z={energia_zeeman / meV:.2f}$ meV")
        eje.set_xlabel(r"$ka/\pi$")

    ejes[0].set_ylabel("Energía (meV)")
    figura.suptitle("Bandas bulk del nanohilo BdG")

    # Leyenda manual para no repetirla en cada panel.
    ejes[0].plot([], [], color="black", linewidth=1.8, label="expresión en $k$")
    ejes[0].scatter([], [], s=12, facecolors="none", edgecolors="black", label="bloques tight-binding")
    ejes[0].legend(loc="best")

    figura.tight_layout()

    return figura


###############################################################################
################################## Gráficas ##################################
###############################################################################

preparar_estilo_graficas()

hamiltoniano_prueba = construir_hamiltoniano_bdg(numero_sitios,potencial_quimico,zeeman_topologico,perfil_gap_inducido(numero_sitios),potencial=potencial_nanohilo_por_defecto(numero_sitios),)

print("Error hermiticidad:", error_hermiticidad(hamiltoniano_prueba))
print("Error simetría partícula-hueco:", error_simetria_particula_hueco(hamiltoniano_prueba, numero_sitios))
print("Campo equivalente al caso topológico [T]:", campo_desde_energia_zeeman(zeeman_topologico))

valores_zeeman = np.linspace(0.0, 0.75 * meV, numero_puntos_barrido)
barrido_zeeman = (valores_zeeman,*barrer_zeeman(numero_sitios,potencial_quimico,valores_zeeman,),)
valores_mu = np.linspace(-0.8 * meV, 0.8 * meV, numero_puntos_barrido)
barrido_mu = (valores_mu,*barrer_potencial_quimico(numero_sitios, valores_mu, zeeman_topologico,),)
valores_zeeman_compacto = np.linspace(0.0, 0.75 * meV, 21)
datos_bulk = (valores_zeeman_compacto,*barrer_gap_bulk_zeeman(valores_zeeman_compacto),)
datos_diagnosticos = (valores_zeeman_compacto, barrer_diagnosticos_modo_zeeman(valores_zeeman_compacto),)
graficar_bandas_bulk_tres_regimenes()
graficar_espectro_finito_vs_zeeman(datos_barrido=barrido_zeeman,)
graficar_energia_minima_vs_zeeman(datos_barrido=barrido_zeeman,)
graficar_espectro_finito_vs_mu(datos_barrido=barrido_mu,)
graficar_energia_minima_vs_mu(datos_barrido=barrido_mu,)
graficar_ldos_tres_regimenes()
graficar_modo_cero_topologico_y_trivial()
graficar_majoranas_separadas()
graficar_particula_hueco_modo_cero()
graficar_perfil_gap_y_potencial()
graficar_splitting_vs_longitud()
graficar_comparacion_gap_localizacion()
graficar_diagrama_fase_mu_zeeman()
graficar_barrido_rashba()
graficar_diagrama_alpha_zeeman()
graficar_gap_bulk_vs_zeeman(datos_bulk=datos_bulk,)
graficar_comparacion_gap_bulk_borde(datos_bulk=datos_bulk,)
graficar_perfil_gap_vs_region_superconductora()
graficar_barrido_temperatura()
graficar_convergencia_paso_malla()
graficar_barrido_barrera_gaussiana()
graficar_mapa_mu_zeeman_minimo()
graficar_peso_bordes_vs_zeeman(datos_diagnosticos=datos_diagnosticos,)
graficar_espin_modo_vs_zeeman(datos_diagnosticos=datos_diagnosticos,)
graficar_longitud_localizacion_vs_zeeman(datos_diagnosticos=datos_diagnosticos,)
graficar_comparacion_gap_homogeneo_inducido()
graficar_bandas_bulk_tres_regimenes()
plt.show()
"""