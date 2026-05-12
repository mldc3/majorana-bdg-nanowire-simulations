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
