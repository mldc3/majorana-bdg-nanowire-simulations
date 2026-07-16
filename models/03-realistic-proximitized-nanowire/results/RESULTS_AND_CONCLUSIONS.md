# Numerical Implementation

*Architecture, algorithms, observables and unchanged Python source*

> The source file reproduced in this volume is copied byte for byte from the uploaded file `nanohilosimetricofinal(1).py`. It is documented as supplied. No function, parameter, comment, plotting instruction or numerical choice has been edited.

**Source length:** 1675 lines  
**SHA-256:** `631cd5ba0cc6d8b73a5ac4b4275dfbe9d4900616c480f8e9fa59210bc8948825`

## Contents

1. [Purpose and computational organization](#chapter-1)
2. [Units and parameter layers](#chapter-2)
   - [Internal units](#internal-units)
   - [Global defaults and configuration dictionary](#global-defaults-and-configuration-dictionary)
   - [Multiple mesh sizes](#multiple-mesh-sizes)
3. [Mesh, energy scales and helper functions](#chapter-3)
4. [Matrix representation](#chapter-4)
5. [Profiles and electrostatics](#chapter-5)
6. [Hamiltonian construction](#chapter-6)
7. [Diagonalization strategies](#chapter-7)
8. [Numerical symmetry checks](#chapter-8)
9. [Electrostatic self-consistency in the supplied source](#chapter-9)
10. [High-level device workflow](#chapter-10)
11. [Physical observables](#chapter-11)
12. [Parameter sweeps and bulk reference](#chapter-12)
13. [Plotting layer](#chapter-13)
14. [Built-in numerical tests](#chapter-14)
15. [Main execution sequence](#chapter-15)
16. [Function-to-equation map](#chapter-16)
17. [Unchanged source listing](#chapter-17)

# Chapter 1

## Purpose and computational organization

The supplied script is a self-contained numerical study of a one-dimensional proximitized semiconductor nanowire. It defines physical constants and model parameters, constructs smooth superconducting and electrostatic profiles, obtains a Hartree-corrected electrostatic base, builds dense or sparse BdG Hamiltonians, extracts low-energy eigenstates, evaluates spectral and spatial diagnostics, performs parameter sweeps and produces the figures used to interpret trivial, transition and topological regimes.

The code is divided by comment banners rather than by separate Python modules. Its logical layers are:

1. constants and parameter configuration;

2. mesh and unit conversion;

3. Pauli and Nambu matrices;

4. proximity and electrostatic profiles;

5. normal and BdG Hamiltonian assembly;

6. dense and sparse diagonalization;

7. Hartree iteration;

8. high-level base construction;

9. observables and Majorana localization;

10. parameter sweeps and bulk diagnostics;

11. plotting functions;

12. numerical tests and the main execution block.

# Chapter 2

## Units and parameter layers

### Internal units

The script sets 

$$
\texttt{meV}=\texttt{nm}=\texttt{kelvin}=1,
$$

<p align="right"><strong>(2.1)</strong></p>

so numerical arrays carry implicit units. Universal constants are converted once: 

$$
k_B=0.08617333262\ \mathrm{meV/K},
$$

<p align="right"><strong>(2.2)</strong></p>

$$
\frac{\hbar^2}{2m_e}=38.0998212\ \mathrm{meV\,nm^2},
$$

<p align="right"><strong>(2.3)</strong></p>

$$
\mu_B=0.05788381806\ \mathrm{meV/T}.
$$

<p align="right"><strong>(2.4)</strong></p>

This makes all energies numerically expressible in meV and all lengths in nm.

### Global defaults and configuration dictionary

The file contains global defaults and a dictionary named `configuracion_realista_suave`. The dictionary is copied by `obtener_configuracion_modelo`; downstream functions use `dict.get` so configuration-specific values override global defaults while unspecified quantities fall back to the globals.

The active dictionary contains a chemical potential, the beginning of superconducting coverage, a barrier height, the covered-region potential shift, a Hartree strength and a representative topological Zeeman energy. This design lets the code pass one compact object through the high-level workflow.

### Multiple mesh sizes

The source deliberately separates meshes by task:

- `numero_sitios_prueba`: inexpensive tests;

- `numero_sitios_hartree`: auxiliary Hartree calculation;

- `numero_sitios_barridos`: parameter sweeps;

- `numero_sitios_estados`: LDOS and near-zero states;

- `numero_sitios_perfiles`: spatial-profile figures.

The main state mesh has 801 sites. A refined electrostatic base is obtained by solving Hartree on a coarser mesh and interpolating the resulting potential to a finer BdG mesh.

# Chapter 3

## Mesh, energy scales and helper functions

### Mesh creation

`crear_malla_1d` returns 

$$
x_i=\operatorname{linspace}(0,L,N),\qquad a=\frac{L}{N-1}.
$$

<p align="right"><strong>(3.1)</strong></p>

The endpoint-inclusive mesh makes the physical wire length exactly $L$.

### Discrete hopping scales

`calcular_hoppings_discretos` evaluates 

$$
t=\frac{\hbar^2}{2m^*a^2},\qquad t_{\mathrm{SO}}=\frac{\alpha_R}{2a}.
$$

<p align="right"><strong>(3.2)</strong></p>

The function accepts a mesh-specific spacing, which is essential because the Hartree, sweep and state meshes use different $a$ values.

### Zeeman-field conversion

`gamma_zeeman_desde_B` and `campo_desde_gamma_zeeman` implement the inverse pair 

$$
\Gamma=\frac{g\mu_BB}{2},\qquad B=\frac{\Gamma}{g\mu_B/2}.
$$

<p align="right"><strong>(3.3)</strong></p>

### Stable Fermi function

`fermi` supports scalars and arrays. At $T\le0$ it returns a step function. At finite temperature it clips $E/(k_BT)$ to $[-700,700]$ before exponentiation, preventing overflow without affecting physically meaningful occupation values.

# Chapter 4

## Matrix representation

### Pauli matrices and Kronecker products

The script first defines $2\times2$ spin and Nambu matrices, then builds the $4\times4$ operators 

$$
\sigma_x=I_\tau\otimes\sigma_x^{(2)},
$$

<p align="right"><strong>(4.1)</strong></p>

$$
\sigma_y=I_\tau\otimes\sigma_y^{(2)},
$$

<p align="right"><strong>(4.2)</strong></p>

$$
\tau_j=\tau_j^{(2)}\otimes I_\sigma.
$$

<p align="right"><strong>(4.3)</strong></p>

The local particle–hole unitary part is 

$$
U_C=\tau_y\sigma_y,
$$

<p align="right"><strong>(4.4)</strong></p>

so the antiunitary transformation is $\mathcal C=U_CK$.

### Local basis

Every four-component site vector is ordered as 

$$
(u_\uparrow,u_\downarrow,v_\downarrow,-v_\uparrow).
$$

<p align="right"><strong>(4.5)</strong></p>

All component extraction, density formulas and symmetry transformations use this order.

# Chapter 5

## Profiles and electrostatics

### Superconducting coverage

`perfil_cobertura_superconductor` evaluates the two-tanh profile and clips the result to $[0,1]$. The right interface is generated as $L-x_L$, so the default geometry is symmetric.

### External potential

`construir_potencial_externo_1d` calculates the covered-region shift and two Gaussian barriers separately, then returns 

$$
(V_{\mathrm{ext}},V_{\mathrm{SC}},V_{\mathrm{barrier}},f_{\mathrm{SC}}).
$$

<p align="right"><strong>(5.1)</strong></p>

Returning the components rather than only the total allows the plotting functions to display and diagnose each contribution.

### Proximity parameters

`calcular_parametros_proximidad` returns 

$$
\gamma(x),\quad\Delta_{\mathrm{ind}}(x),\quad Z(x),\quad f_{\mathrm{SC}}(x).
$$

<p align="right"><strong>(5.2)</strong></p>

A small positive denominator regularizer avoids division by zero when a scale is set to zero.

### Screened interaction matrix

`construir_matriz_interaccion_apantallada` uses broadcasting to create the dense distance matrix $|x_i-x_j|$ and evaluates 

$$
\nu_{ij}=U_0\exp[-|x_i-x_j|/\lambda_{\mathrm{scr}}].
$$

<p align="right"><strong>(5.3)</strong></p>

This operation scales as $O(N^2)$ in memory and time.

### Hartree update

`actualizar_hartree` forms the charge relative to a supplied reference, optionally subtracts its mean and multiplies by the interaction matrix. When the interaction matrix is numerically zero, it returns a zero potential directly.

# Chapter 6

## Hamiltonian construction

### Normal-state Hamiltonian

`construir_H_normal_1d` builds a $2N\times2N$ spinful Hamiltonian with local blocks 

$$
h_{ii}=(2t-\mu+V_{\mathrm{eff},i})I_2+\Gamma\sigma_x
$$

<p align="right"><strong>(6.1)</strong></p>

and nearest-neighbour blocks 

$$
h_{i,i+1}=-tI_2+\mathrm{i} t_{\mathrm{SO}}\sigma_y.
$$

<p align="right"><strong>(6.2)</strong></p>

The final matrix is explicitly symmetrized as $(H+H^\dagger)/2$ to suppress round-off-level non-Hermiticity.

### Dense BdG Hamiltonian

`construir_hamiltoniano_bdg_realista` creates a dense $4N\times4N$ array. Its local block is 

$$
H_{ii}=Z_i(2t-\mu+V_{\mathrm{eff},i})\tau_z+Z_i\Gamma\sigma_x +\Re\Delta_i\tau_x-\Im\Delta_i\tau_y,
$$

<p align="right"><strong>(6.3)</strong></p>

and its bond block is 

$$
H_{i,i+1}=\sqrt{Z_iZ_{i+1}} (-t\tau_z+\mathrm{i} t_{\mathrm{SO}}\sigma_y\tau_z).
$$

<p align="right"><strong>(6.4)</strong></p>

The opposite block is assigned as the conjugate transpose.

### Sparse BdG Hamiltonian

`construir_hamiltoniano_bdg_realista_disperso` assembles the identical block structure in LIL format and converts it to CSR. The sparse representation stores only $O(N)$ nonzero blocks, which makes near-zero eigensolvers feasible for the high-resolution meshes.

# Chapter 7

## Diagonalization strategies

### Dense eigensolvers

`diagonalizar_hamiltoniano` uses `numpy.linalg.eigh` and returns eigenvalues and eigenvectors sorted in ascending energy. `diagonalizar_energias` uses `eigvalsh` when vectors are unnecessary.

### Shift-invert near zero

`diagonalizar_cerca_cero_bdg` switches to `scipy.sparse.linalg.eigsh` when the dimension exceeds 300. With 

$$
\texttt{sigma}=0,
$$

<p align="right"><strong>(7.1)</strong></p>

shift-invert returns eigenvalues closest to zero rather than those of largest magnitude. If exact zero causes a factorization failure, the function retries with a small nonzero shift. For small systems it diagonalizes densely and selects the requested number of states by $|E|$.

### Complexity

A full dense Hermitian diagonalization scales as $O((4N)^3)$ and stores $O((4N)^2)$ entries. Sparse shift-invert avoids the full spectrum and is therefore the central method for LDOS and end-state diagnostics. The Hartree matrix remains dense, so the two-mesh strategy reduces the most expensive electrostatic part.

# Chapter 8

## Numerical symmetry checks

### Hermiticity residual

`error_hermiticidad` calculates 

$$
\epsilon_H=\frac{\left\lVert H-H^\dagger\right\rVert}{\left\lVert H\right\rVert}.
$$

<p align="right"><strong>(8.1)</strong></p>

Sparse and dense matrices are handled separately.

### Matrix particle–hole residual

For matrices small enough to convert to dense form, `error_simetria_particula_hueco` builds $I_N\otimes U_C$ and evaluates 

$$
\epsilon_C=\frac{\left\lVert U_CH^*U_C^\dagger+H\right\rVert}{\left\lVert H\right\rVert}.
$$

<p align="right"><strong>(8.2)</strong></p>

Large matrices return `NaN` to avoid an expensive dense conversion.

### Spectral pairing residual

`error_espectro_particula_hueco` sorts the energies and compares opposite ends of the list: 

$$
\epsilon_E=\max_n|E_n+E_{-n}|.
$$

<p align="right"><strong>(8.3)</strong></p>

This check can be applied even when only a symmetric near-zero subset has been calculated.

# Chapter 9

## Electrostatic self-consistency in the supplied source

### Normal-state density

`densidad_normal_1d` sums occupied normal-state eigenvectors: 

$$
n_i=\sum_n f(E_n)\left(|\psi_{n\uparrow}(i)|^2+|\psi_{n\downarrow}(i)|^2\right).
$$

<p align="right"><strong>(9.1)</strong></p>

### Diagnostic BdG density

`densidad_bdg_para_hartree` separately computes, from the available positive-energy BdG states, 

$$
n_i^{\mathrm{BdG}}=\sum_{E_n>0}\left[\rho_{u,n}(i)f(E_n)+\rho_{v,n}(i)(1-f(E_n))\right].
$$

<p align="right"><strong>(9.2)</strong></p>

The returned metadata records whether the eigenvectors form only a subspace of the full $4N$ spectrum.

### Hartree loop

`resolver_Hartree_normal` initializes $U_H=0$ and uses a zero reference density. At each step it constructs the normal Hamiltonian, diagonalizes the full $2N\times2N$ matrix, computes the normal density, evaluates the screened Hartree potential and mixes it with the previous iterate. It stores both a relative norm error and the maximum local change.

The return dictionary contains the final spectrum, eigenvectors, Hamiltonian, density, $U_H$, $V_{\mathrm{eff}}$, convergence histories, a Boolean flag and the number of iterations.

# Chapter 10

## High-level device workflow

### Electrostatic base

`resolver_base_electrostatica` performs the following operations:

1. load the configuration;

2. create the mesh and hopping scales;

3. construct $V_{\mathrm{ext}}$;

4. solve the normal Hartree loop unless $U_0=0$;

5. construct $\gamma$, $\Delta_{\mathrm{ind}}$ and $Z$;

6. return all arrays and metadata in one dictionary.

### Refined base

`resolver_base_electrostatica_refinada` solves Hartree on a coarse mesh, constructs a fine base with Hartree disabled, interpolates $U_H$ and the density onto the fine mesh and rebuilds $V_{\mathrm{eff}}$. The proximity and external-potential profiles on the fine mesh are evaluated directly rather than interpolated.

### BdG solution from a base

`resolver_bdg_desde_base` calls the near-zero eigensolver with the converged electrostatic arrays and returns the low-energy spectrum, eigenvectors, Hamiltonian, diagnostic BdG density and Hartree convergence metadata.

### Three representative regimes

`calcular_resultados_regimenes` uses one base for all three fields. `calcular_resultados_regimenes_desde_bases` instead accepts separate electrostatic bases for the trivial, transition and topological fields. The final main block uses the second route.

# Chapter 11

## Physical observables

### Component extraction and normalized density

`separar_componentes_bdg` reshapes a vector into $N\times4$. `densidad_estado` sums all four absolute squares and divides by the total norm. `densidad_electron_hueco_total` returns normalized electron, hole and total densities separately.

### Selecting low-energy states

The functions `indice_estado_mas_cercano_a_cero`, `indice_estado_positivo_mas_cercano_a_cero` and `energia_minima_absoluta` provide consistent state selection. `pesos_bordes` integrates the normalized density within a default $150\,\mathrm{nm}$ interval at each end.

### Particle–hole transformation

`transformar_particula_hueco` applies $U_CK$ site by site. `alinear_fase` removes a global relative phase by making the overlap with a reference state real and positive.

### Localized Majorana combinations

`construir_majoranas_desde_pareja` selects the two states closest to zero and forms their $2\times2$ matrix of left-half weight. Diagonalizing this matrix finds the linear combinations with minimum and maximum left weight. They are labeled right and left, converted to normalized spatial densities and returned together with the two source energies.

This subspace-localization method is insensitive to an arbitrary unitary rotation of nearly degenerate eigenvectors and is therefore more stable than assigning one raw eigensolver vector to each edge.

### LDOS-like symmetric spectral map

`calcular_ldos` uses the normalized total BdG density of each positive-energy state and adds Lorentzians at both $+E_n$ and $-E_n$: 

$$
\rho(i,E)=\sum_{E_n>0}\rho_n(i)\left[L_\eta(E-E_n)+L_\eta(E+E_n)\right].
$$

<p align="right"><strong>(11.1)</strong></p>

The plotting function can display the logarithm of this quantity to reveal weak subgap weight.

# Chapter 12

## Parameter sweeps and bulk reference

### Zeeman sweep

`barrer_zeeman_desde_base` repeatedly solves the low-energy BdG problem, stores the full near-zero subset, the minimum $|E|$ and the left/right edge weights. An optional flag can recalculate the electrostatic base at every Zeeman value.

### Chemical-potential sweep

`barrer_mu_desde_base` reuses one electrostatic base and rebuilds only the BdG Hamiltonian for each chemical potential. It stores the near-zero spectra and minimum energies.

### Bulk Hamiltonian

`construir_hamiltoniano_bulk_k_realista` evaluates 

$$
H(k)=H_0+T\mathrm{e}^{\mathrm{i} ka}+T^\dagger\mathrm{e}^{-\mathrm{i} ka}
$$

<p align="right"><strong>(12.1)</strong></p>

from covered-region averages. `calcular_bandas_bulk_realistas` diagonalizes this $4\times4$ matrix on a momentum grid.

### Critical estimates

`parametros_bulk_desde_base` averages $V_{\mathrm{eff}}$, $\Delta_{\mathrm{ind}}$ and $Z$ where $f_{\mathrm{SC}}>0.9$. The critical functions implement 

$$
\Gamma_c=\sqrt{\mu_{\mathrm{eff}}^2+(\Delta_{\mathrm{ind}}/Z)^2}
$$

<p align="right"><strong>(12.2)</strong></p>

and 

$$
\mu_c^\pm=V_{\mathrm{eff}}\pm\sqrt{\Gamma^2-(\Delta_{\mathrm{ind}}/Z)^2}.
$$

<p align="right"><strong>(12.3)</strong></p>

These values define the dashed lines and shaded topological intervals in the figures.

# Chapter 13

## Plotting layer

The plotting functions are intentionally separate from the numerical solvers. They receive already calculated bases or result dictionaries and create:

- full and zoomed bulk-band comparisons;

- finite subgap levels in the three representative regimes;

- symmetric LDOS maps;

- left/right Majorana decompositions;

- electron–hole content;

- finite spectra versus $\Gamma$ and $\mu$;

- superconducting, potential and density profiles.

The helper `aplicar_cero_visual_a_modo_topologico` can place sufficiently small topological energies exactly at zero for visual display while the diagnostic functions retain the unmodified numerical value.

# Chapter 14

## Built-in numerical tests

`ejecutar_pruebas_basicas` constructs a small test system and checks:

1. Hermiticity of the BdG matrix;

2. matrix particle–hole symmetry;

3. $\pm E$ spectral pairing;

4. recovery of the normal spectrum when $\Delta=0$;

5. vanishing induced gap for zero interface coupling;

6. saturation toward the parent gap at strong coupling;

7. vanishing Hartree potential for $U_0=0$.

Each residual is compared with an internal threshold and printed as `OK` or `ADVERTENCIA`.

# Chapter 15

## Main execution sequence

When executed as a script, the file:

1. prepares the Matplotlib style and runs basic tests;

2. prints the active mesh sizes;

3. constructs a sweep base and a fine state base;

4. estimates critical Zeeman values;

5. recalculates separate Hartree bases for the three regimes;

6. solves the low-energy BdG problem for each regime;

7. performs Zeeman and chemical-potential sweeps;

8. generates the complete figure set;

9. prints spectral, LDOS and localization diagnostics;

10. displays all figures.

The workflow is deterministic apart from ordinary eigensolver phase choices and possible rotations inside numerically degenerate subspaces; physical densities are insensitive to a global phase.

# Chapter 16

## Function-to-equation map

| **Function**                                   | **Mathematical role**                                      |
|:-----------------------------------------------|:-----------------------------------------------------------|
| `crear_malla_1d`                               | $x_i=ia$, $a=L/(N-1)$.                                     |
| `calcular_hoppings_discretos`                  | $t=\hbar^2/(2m^*a^2)$ and $t_{\mathrm{SO}}=\alpha_R/(2a)$. |
| `perfil_cobertura_superconductor`              | Two-tanh superconducting coverage.                         |
| `calcular_parametros_proximidad`               | $\gamma$, $\Delta_{\mathrm{ind}}$, $Z$.                    |
| `construir_potencial_externo_1d`               | Global offset, covered-region shift and Gaussian barriers. |
| `construir_matriz_interaccion_apantallada`     | $\nu_{ij}=U_0e^{-|x_i-x_j|/\lambda}$.                      |
| `resolver_Hartree_normal`                      | Fixed-point iteration for $U_H$ from normal-state density. |
| `construir_hamiltoniano_bdg_realista`          | Dense block-tridiagonal BdG matrix.                        |
| `construir_hamiltoniano_bdg_realista_disperso` | Sparse CSR BdG matrix.                                     |
| `diagonalizar_cerca_cero_bdg`                  | Shift-invert low-energy eigensolver.                       |
| `calcular_ldos`                                | Symmetric Lorentzian BdG spectral map.                     |
| `construir_majoranas_desde_pareja`             | Left/right localization inside the near-zero subspace.     |
| `barrer_zeeman_desde_base`                     | Finite spectrum as a function of $\Gamma$.                 |
| `barrer_mu_desde_base`                         | Finite spectrum as a function of $\mu$.                    |
| `construir_hamiltoniano_bulk_k_realista`       | Homogeneous $4\times4$ lattice Hamiltonian.                |
| `estimar_gamma_critico_desde_base`             | Renormalized homogeneous critical estimate.                |


# Chapter 17

## Unchanged source listing

The following listing is the supplied Python file in full. Line numbers are added by LaTeX during typesetting; they are not part of the stored source file.

````python
import os
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["OMP_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"
import numpy as np
import numpy.linalg as la
import matplotlib.pyplot as plt
import scipy.sparse as sp
import scipy.sparse.linalg as spla
from matplotlib.lines import Line2D
scipy_sparse_disponible = True

###############################################################################
################ Constantes sacados de literatura o papers ####################
###############################################################################

# Definimos las unidades que vamos a utilizar en la simulacion
meV = 1.0
nm = 1.0
kelvin = 1.0

#Definición de constantes universales
constante_boltzmann = 0.08617333262 * meV / kelvin          # k_B en meV/K
hbar2_sobre_2m_e = 38.0998212 * meV * nm**2                 # hbar^2/(2m_e) en meV nm^2
magneton_bohr = 0.05788381806 * meV                         # mu_B en meV/T

# Parametros de InSb (usados por Woods, Stanescu y Das Sarma)
masa_efectiva = 0.014                                       # m*/m_e para InSb (Woods, Stanescu & Das Sarma, PRB 2018)
epsilon_relativa = 17.7                                     # Permitividad relativa de InSb (Woods, Stanescu & Das Sarma, PRB 2018)
factor_g = 40.0                                             # Factor g efectivo grande típico de InSb; valor fenomenologico del modelo 1 (Mourik et al., Science 2012)
alpha_rashba = 25.0 * meV * nm                              # 250 meV A = 25 meV nm, valor usado para hilos finitos (Woods, Stanescu & Das Sarma)

potencial_quimico = 0.18 * meV                              # mu efectivo respecto a la subbanda activa (este parametro se ajusta)
temperatura = 0.02 * kelvin                                 # Temperatura numerica
Delta0_padre = 0.30 * meV                                   # Gap del superconductor padre usado como escala en el modelo de proximidad (Woods, Stanescu & Das Sarma, PRB 2018)
gamma0_interfaz = 0.30 * meV                                # Acoplo SM-SC maximo; parametro efectivo elegido para obtener Delta_ind approx 0.15 meV cuando gamma0 = Delta0

# Electrostatica externa efectiva. Los valores grandes de Woods et al. viven en el problema 3D; aqui se usan solo para fijar escalas del modelo 1D
E0_efectivo = 0.0 * meV                                     # Offset energetico ya absorbido/renormalizado en el modelo 1D
V_gate_global = 0.0 * meV                                   # Gate global efectivo usado como control electrostático uniforme
amplitud_VSC_efectiva = -0.04 * meV                         # Desplazamiento efectivo de interfaz SM-SC bajo la cobertura, valor ajustado del modelo 1D
altura_barrera = 0.006 * meV                                 # Barrera muy suave/casi transparente: evita ABS triviales de borde y mantiene un perfil experimental limpio (Pan et al., PRB 2021)
centro_barrera = 75.0 * nm                                   # Centro de la barrera suave cerca de cada extremo; no debe crear puntos cuánticos triviales
anchura_barrera = 55.0 * nm                                  # Anchura espacial de las barreras suaves gaussianas

# Geometria efectiva de la zona cubierta por superconductor
x_inicio_superconductor = 45.0 * nm                          # Cobertura SC casi completa: elimina normal dots triviales en los extremos y deja Majoranas en los bordes físicos
suavidad_superconductor = 18.0 * nm                          # Suavidad de las transiciones normales-SC simétricas

# Hartree 1D efectivo. nu_ij imita el tensor de Green apantallado del problema 3D
U0_hartree = 0.045 * meV                                     # Intensidad Hartree efectiva 1D moderada: realista sin generar ABS triviales de confinamiento
lambda_apantallamiento = 45.0 * nm                          # Longitud de apantallamiento efectiva
mezcla_hartree = 0.18                                       # Mezcla Hartree más rápida; estabilidad numérica para el modelo 1D efectivo
tolerancia_hartree = 5.0e-6                                 # Criterio de convergencia Hartree suficiente para perfiles/figuras de TFG
max_iteraciones_hartree = 55                                # Límite de iteraciones Hartree para evitar cálculos innecesariamente largos
neutralizar_carga_media = True                              # Resta la carga media para que un fondo uniforme no redefina solo mu globalmente

# Casos representativos. Con Z != 1 el cierre estimado usa Delta_ind/Z
delta_inducido_bulk_estimado = gamma0_interfaz * Delta0_padre / (gamma0_interfaz + Delta0_padre)
Z_bulk_estimado = 1.0 / (1.0 + gamma0_interfaz / Delta0_padre)
zeeman_trivial = 0.0 * meV                                  # Caso trivial sin campo Zeeman
zeeman_transicion_estimado = np.sqrt(potencial_quimico**2 + (delta_inducido_bulk_estimado / Z_bulk_estimado)**2)             # Estimacion uniforme de la transicion topologica
zeeman_transicion = zeeman_transicion_estimado              # Caso de transición usado para graficar la zona de cierre/reapertura
zeeman_topologico = 0.62 * meV                              # Campo Zeeman elegido por encima de la transicion estimada para visualizar el regimen topologico





###############################################################################
################## Valores elegidos para la discretización ####################
###############################################################################

# Malla longitudinal. Se separan tamaños porque no todas las figuras necesitan diagonalizar matrices enormes
longitud_hilo = 4000.0 * nm                                  # 4 micras: reduce mucho el solapamiento entre Majoranas y hace más limpio el modo cero
numero_sitios_prueba = 101                                  # Malla pequeña para comprobar rápido el perfil simétrico
numero_sitios_hartree = 251                                 # Malla auxiliar para Hartree denso: suficiente para perfiles suaves sin disparar coste
numero_sitios_barridos = 601                                # Malla de barridos: más larga y fina para que el modo cero topológico no oscile visualmente
numero_sitios_estados = 801                                 # Malla principal de estados/LDOS/Majoranas: compromiso alto entre física y coste
numero_sitios_perfiles = 801                                # Malla de perfiles físicos: figuras suaves y simétricas para el TFG
numero_sitios_final = numero_sitios_estados                 # Resolución alta principal para figuras espaciales finales
numero_sitios = numero_sitios_estados                       # Malla principal del script: usa la resolución alta de estados
paso_espacial = longitud_hilo / (numero_sitios - 1)         # Paso espacial asociado a numero_sitios; con N=10001 y L=2000 nm, a = 0.2 nm

# Escalas discretas para el hamiltoniano 
hopping = hbar2_sobre_2m_e / (masa_efectiva * paso_espacial**2)
hopping_spin_orbita = alpha_rashba / (2.0 * paso_espacial)

# Barridos y visualizacion
numero_energias = 440
ensanchamiento_ldos = 0.025 * meV
numero_puntos_barrido = 71
umbral_cero_visual = 2.5e-3 * meV                         # Solo para gráficas: estados topológicos con |E| menor que esto se dibujan en cero, pero los diagnósticos imprimen el valor real

# Crea la malla longitudinal 1D del nanohilo
def crear_malla_1d(longitud_usada=longitud_hilo, numero_sitios_usado=numero_sitios):
    x = np.linspace(0.0, longitud_usada, numero_sitios_usado)
    paso = longitud_usada / max(numero_sitios_usado - 1, 1)
    return x, paso

# Calcula t y t_so para una malla concreta
def calcular_hoppings_discretos(alpha_rashba_usado, paso_espacial_usado, masa_efectiva_usada=masa_efectiva):
    hopping_usado = hbar2_sobre_2m_e / (masa_efectiva_usada * paso_espacial_usado**2) # t = hbar^2/(2m*a^2)
    hopping_so_usado = alpha_rashba_usado / (2.0 * paso_espacial_usado)               # t_so = alpha_R/(2a)
    return hopping_usado, hopping_so_usado


# Configuración física final: modelo realista suave autoconsistente
configuracion_realista_suave = {
    "nombre": "realista_suave",
    "potencial_quimico": 0.12 * meV,                         # Valor efectivo ajustado para el modelo 1D; controla la ocupación de la subbanda activa
    "x_inicio_superconductor": 45.0 * nm,                    # Perfil simétrico casi completamente cubierto: evita estados Andreev triviales de borde
    "altura_barrera": 0.006 * meV,                           # Barreras suaves casi transparentes: perfil experimental sin generar falsos Majoranas triviales
    "amplitud_VSC_efectiva": -0.04 * meV,                    # Desplazamiento medio de interfaz SM-SC bajo la cobertura superconductora
    "U0_hartree": 0.045 * meV,                               # Hartree efectivo 1D moderado; mantiene autoconsistencia sin sobreatrapar estados
    "zeeman_topologico": 0.62 * meV                           # Campo por encima de la transición: modo cero más localizado y separación topológica más clara
}

# Devuelve una copia independiente de la configuración realista suave
def obtener_configuracion_modelo(nombre="realista_suave"):
    if isinstance(nombre, dict):
        return nombre.copy()
    if nombre != "realista_suave":
        raise ValueError("La versión final solo conserva la configuración realista_suave")
    return configuracion_realista_suave.copy()

# Obtenemos el Zeeman critico de una configuracion usando la formula uniforme: Gamma_c = sqrt(mu^2 + (Delta_ind/Z)^2)
def estimar_zeeman_transicion_configuracion(configuracion):
    gamma0 = configuracion.get("gamma0_interfaz", gamma0_interfaz)
    Delta0 = configuracion.get("Delta0_padre", Delta0_padre)
    mu = configuracion.get("potencial_quimico", potencial_quimico)
    delta_bulk = gamma0 * Delta0 / (gamma0 + Delta0 + 1.0e-30)      # Gap inducido estimado en la region cubierta por superconductor
    Z_bulk = 1.0 / (1.0 + gamma0 / (Delta0 + 1.0e-30))              # Renormalizacion de las escalas normales por proximidad
    return float(np.sqrt(mu**2 + (delta_bulk / Z_bulk)**2))         # Criterio topologico uniforme aproximado

# Convierte campo magnetico en Gamma = g mu_B B / 2
def gamma_zeeman_desde_B(campo_magnetico, factor_g_usado=factor_g):
    return 0.5 * factor_g_usado * magneton_bohr * campo_magnetico

# Convierte Gamma en el campo magnetico equivalente
def campo_desde_gamma_zeeman(gamma_zeeman, factor_g_usado=factor_g):
    return gamma_zeeman / (0.5 * factor_g_usado * magneton_bohr)

# Funcion de Fermi-Dirac estable numericamente: f(E)=1/(exp(E/kBT)+1)
def fermi(energia, temperatura_usada=temperatura):
    energia = np.asarray(energia, dtype=float)                 # Convierte la energia a array de NumPy para poder usar escalares o vectores
    if temperatura_usada <= 0.0:                               # Si T <= 0, se usa el limite de temperatura cero
        return (energia < 0.0).astype(float)                   # A T=0, los estados con E<0 estan ocupados y los de E>0 vacios
    argumento = energia / (constante_boltzmann * temperatura_usada)  # Calcula E/(k_B T), que es el argumento de la exponencia
    argumento = np.clip(argumento, -700.0, 700.0)              # Recorta el argumento para evitar overflow numerico en exp(argumento)
    return 1.0 / (np.exp(argumento) + 1.0)                    

###############################################################################
############################## Matrices de Pauli #############################
###############################################################################

matriz_identidad_2 = np.eye(2, dtype=complex)
matriz_sigma_x_2 = np.array([[0.0, 1.0], [1.0, 0.0]], dtype=complex)
matriz_sigma_y_2 = np.array([[0.0, -1.0j], [1.0j, 0.0]], dtype=complex)
matriz_sigma_z_2 = np.array([[1.0, 0.0], [0.0, -1.0]], dtype=complex)
matriz_tau_x_2 = np.array([[0.0, 1.0], [1.0, 0.0]], dtype=complex)
matriz_tau_y_2 = np.array([[0.0, -1.0j], [1.0j, 0.0]], dtype=complex)
matriz_tau_z_2 = np.array([[1.0, 0.0], [0.0, -1.0]], dtype=complex)

# Base compacta BdG: (u_up, u_down, v_down, -v_up). tau actua en particula-hueco y sigma en espin
sigma_x = np.kron(matriz_identidad_2, matriz_sigma_x_2)
sigma_y = np.kron(matriz_identidad_2, matriz_sigma_y_2)
tau_x = np.kron(matriz_tau_x_2, matriz_identidad_2)
tau_y = np.kron(matriz_tau_y_2, matriz_identidad_2)
tau_z = np.kron(matriz_tau_z_2, matriz_identidad_2)
sigma_y_tau_z = sigma_y @ tau_z
operador_particula_hueco_local = tau_y @ sigma_y


###############################################################################
######################### Perfiles, proximidad y electrostatica ###############
###############################################################################

# Perfil suave simetrico de cobertura SC:
# fSC(x)=1/2[tanh((x-xL)/sSC)-tanh((x-xR)/sSC)]
def perfil_cobertura_superconductor(x, x_inicio=x_inicio_superconductor, suavidad=suavidad_superconductor):
    x = np.asarray(x, dtype=float)                            # Malla longitudinal del nanohilo
    x_izquierda = x_inicio                                    # Inicio suave de la cobertura SC
    x_derecha = np.max(x) - x_inicio                          # Final suave de la cobertura SC, elegido de forma simetrica
    if x_izquierda <= np.min(x) and x_derecha >= np.max(x):    # Si la cobertura abarca todo el hilo
        return np.ones_like(x, dtype=float)                    # Devuelve cobertura uniforme
    fSC = 0.5 * (np.tanh((x - x_izquierda) / suavidad) - np.tanh((x - x_derecha) / suavidad)) # Ventana SC suave y simetrica
    return np.clip(fSC, 0.0, 1.0)                              # Evita pequeños sobrepasos numericos

# Potencial externo efectivo 1D:
# V_ext(x)=E0+Vg+VSC^0 fSC(x)+Vb^0 exp[-(x-xL)^2/(2 sigma_b^2)]+Vb^0 exp[-(x-xR)^2/(2 sigma_b^2)]
def construir_potencial_externo_1d(x, E0=E0_efectivo, Vg=V_gate_global, amplitud_VSC=amplitud_VSC_efectiva, Vb=altura_barrera, xb=centro_barrera, sigma_b=anchura_barrera, x_inicio=x_inicio_superconductor, suavidad=suavidad_superconductor):
    x = np.asarray(x, dtype=float)                            # Malla longitudinal del nanohilo
    fSC = perfil_cobertura_superconductor(x, x_inicio, suavidad)               # Cobertura SC central simetrica
    VSC_x = amplitud_VSC * fSC                                                 # Potencial de interfaz bajo la zona cubierta
    centro_izquierdo = xb                                                       # Centro de la barrera izquierda
    centro_derecho = np.max(x) - xb                                             # Centro de la barrera derecha por simetria
    barrera_izquierda = Vb * np.exp(-0.5 * ((x - centro_izquierdo) / sigma_b)**2) # Barrera gaussiana suave izquierda
    barrera_derecha = Vb * np.exp(-0.5 * ((x - centro_derecho) / sigma_b)**2)     # Barrera gaussiana suave derecha
    barrera = barrera_izquierda + barrera_derecha                              # Potencial simetrico de dos barreras
    V_ext = E0 + Vg + VSC_x + barrera                                          # Potencial externo total
    return V_ext, VSC_x, barrera, fSC                                          
# Parametros de proximidad: gamma(x), Delta_ind(x)=gamma*Delta0/(gamma+Delta0) y Z(x)=1/(1+gamma/Delta0) (modelo efectivo de autoenergia SM-SC usado en literatura de nanohilos Majorana)
def calcular_parametros_proximidad(x, gamma0=gamma0_interfaz, Delta0=Delta0_padre, x_inicio=x_inicio_superconductor, suavidad=suavidad_superconductor):
    fSC = perfil_cobertura_superconductor(x, x_inicio, suavidad)               # Calcula el perfil de cobertura superconductora
    gamma = gamma0 * fSC                                                       # Acoplo local SM-SC: cero en region normal y gamma0 bajo el SC
    Delta_ind = gamma * Delta0 / (gamma + Delta0 + 1.0e-30)                    # Gap inducido local por proximidad
    Z = 1.0 / (1.0 + gamma / (Delta0 + 1.0e-30))                               # Renormalizacion local de las escalas normales
    return gamma, Delta_ind, Z, fSC                                           

# Interaccion Hartree 1D apantallada nu_ij = U0 exp(-|xi-xj|/lambda_scr) (aproximacion 1D al tensor/Green electrostatico apantallado de Woods et al.; no es el calculo 3D literal)
def construir_matriz_interaccion_apantallada(x, U0=U0_hartree, lambda_scr=lambda_apantallamiento):
    distancia = np.abs(x[:, None] - x[None, :])                                # Matriz de distancias |xi-xj| entre todos los sitios
    return U0 * np.exp(-distancia / lambda_scr)                                # Interaccion efectiva que decae exponencialmente con la distancia

# Actualiza el potencial Hartree a partir de la densidad y de un fondo de referencia
def actualizar_hartree(densidad, densidad_referencia, matriz_interaccion, neutralizar_media=neutralizar_carga_media):
    carga_efectiva = np.asarray(densidad - densidad_referencia, dtype=float)    # Calcula la carga relativa respecto a la referencia
    if neutralizar_media:                                                       # Si se activa la neutralizacion global
        carga_efectiva = carga_efectiva - np.mean(carga_efectiva)               # Resta la carga media para evitar un desplazamiento global artificial de mu
    if np.max(np.abs(matriz_interaccion)) < 1.0e-30:                            # Si la interaccion es practicamente nula
        return np.zeros_like(carga_efectiva)                                    # Devuelve Hartree cero
    with np.errstate(over="ignore", invalid="ignore", divide="ignore"):         # Evita warnings numericos no relevantes en la multiplicacion
        return matriz_interaccion @ carga_efectiva                              # Calcula U_H,i = sum_j nu_ij (n_j-n_ref_j)

###############################################################################
############################ Hamiltonianos y diagonalizacion ##################
###############################################################################

# Hamiltoniano normal 1D:
# h_i = (2t - mu + V_eff,i) sigma_0 + Gamma sigma_x
# T_i,i+1 = -t sigma_0 + i t_so sigma_y
def construir_H_normal_1d(numero_sitios_usado, potencial_quimico_usado, energia_zeeman_usada, V_eff, hopping_usado, hopping_spin_orbita_usado):
    V_eff = np.asarray(V_eff, dtype=float)                    # Asegura que el potencial sea un array numerico
    if len(V_eff) != numero_sitios_usado:                     # Comprueba que V_eff tenga un valor por sitio
        raise ValueError("El perfil V_eff no tiene el mismo numero de sitios que el hilo.")
    hamiltoniano = np.zeros((2 * numero_sitios_usado, 2 * numero_sitios_usado), dtype=complex)  # Matriz normal: 2 grados de libertad de espin por sitio
    bloque_hopping = -hopping_usado * matriz_identidad_2 + 1.0j * hopping_spin_orbita_usado * matriz_sigma_y_2
    for i in range(numero_sitios_usado):
        bloque_i = slice(2 * i, 2 * (i + 1))                  # Bloque de espin del sitio i
        bloque_local = (2.0 * hopping_usado - potencial_quimico_usado + V_eff[i]) * matriz_identidad_2 + energia_zeeman_usada * matriz_sigma_x_2
        hamiltoniano[bloque_i, bloque_i] = bloque_local
        if i < numero_sitios_usado - 1:
            bloque_j = slice(2 * (i + 1), 2 * (i + 2))        # Bloque de espin del sitio vecino i+1
            hamiltoniano[bloque_i, bloque_j] = bloque_hopping
            hamiltoniano[bloque_j, bloque_i] = bloque_hopping.conj().T
    return 0.5 * (hamiltoniano + hamiltoniano.conj().T)       # Simetriza para eliminar errores numericos pequeños


# Hamiltoniano BdG realista:
# h_i = Z_i(2t - mu + V_eff,i) tau_z + Z_i Gamma sigma_x + Re(Delta_i) tau_x - Im(Delta_i) tau_y
# T_i,i+1 = sqrt(Z_i Z_i+1)(-t tau_z + i t_so sigma_y tau_z)
def construir_hamiltoniano_bdg_realista(numero_sitios_usado, potencial_quimico_usado, energia_zeeman_usada, V_eff, Delta_ind, Z, hopping_usado, hopping_spin_orbita_usado):
    V_eff = np.asarray(V_eff, dtype=float)                    # Potencial efectivo que entra en el Hamiltoniano
    Delta_ind = np.asarray(Delta_ind, dtype=complex)          # Gap inducido local, permitido complejo
    Z = np.asarray(Z, dtype=float)                            # Factor de renormalizacion local
    hamiltoniano = np.zeros((4 * numero_sitios_usado, 4 * numero_sitios_usado), dtype=complex)  # Matriz BdG: 4 componentes por sitio
    for i in range(numero_sitios_usado):
        bloque_i = slice(4 * i, 4 * (i + 1))                  # Bloque BdG del sitio i
        Zi = Z[i]                                             # Renormalizacion local
        Di = Delta_ind[i]                                     # Gap inducido local
        bloque_local = Zi * (2.0 * hopping_usado - potencial_quimico_usado + V_eff[i]) * tau_z + Zi * energia_zeeman_usada * sigma_x + Di.real * tau_x - Di.imag * tau_y
        hamiltoniano[bloque_i, bloque_i] = bloque_local
        if i < numero_sitios_usado - 1:
            bloque_j = slice(4 * (i + 1), 4 * (i + 2))        # Bloque BdG del sitio vecino i+1
            Zij = np.sqrt(max(Z[i] * Z[i + 1], 0.0))          # Promedio geometrico de la renormalizacion entre vecinos
            bloque_hopping = Zij * (-hopping_usado * tau_z + 1.0j * hopping_spin_orbita_usado * sigma_y_tau_z)
            hamiltoniano[bloque_i, bloque_j] = bloque_hopping
            hamiltoniano[bloque_j, bloque_i] = bloque_hopping.conj().T
    return 0.5 * (hamiltoniano + hamiltoniano.conj().T)       # Fuerza hermiticidad numerica


# Misma matriz BdG que construir_hamiltoniano_bdg_realista, pero guardada como matriz dispersa para N grande
def construir_hamiltoniano_bdg_realista_disperso(numero_sitios_usado, potencial_quimico_usado, energia_zeeman_usada, V_eff, Delta_ind, Z, hopping_usado, hopping_spin_orbita_usado):
    V_eff = np.asarray(V_eff, dtype=float)
    Delta_ind = np.asarray(Delta_ind, dtype=complex)
    Z = np.asarray(Z, dtype=float)
    hamiltoniano = sp.lil_matrix((4 * numero_sitios_usado, 4 * numero_sitios_usado), dtype=complex)  # Formato eficiente para rellenar por bloques
    for i in range(numero_sitios_usado):
        bloque_i = slice(4 * i, 4 * (i + 1))
        Zi = Z[i]
        Di = Delta_ind[i]
        bloque_local = Zi * (2.0 * hopping_usado - potencial_quimico_usado + V_eff[i]) * tau_z + Zi * energia_zeeman_usada * sigma_x + Di.real * tau_x - Di.imag * tau_y
        hamiltoniano[bloque_i, bloque_i] = bloque_local
        if i < numero_sitios_usado - 1:
            bloque_j = slice(4 * (i + 1), 4 * (i + 2))
            Zij = np.sqrt(max(Z[i] * Z[i + 1], 0.0))
            bloque_hopping = Zij * (-hopping_usado * tau_z + 1.0j * hopping_spin_orbita_usado * sigma_y_tau_z)
            hamiltoniano[bloque_i, bloque_j] = bloque_hopping
            hamiltoniano[bloque_j, bloque_i] = bloque_hopping.conj().T
    return hamiltoniano.tocsr()                              # Convierte a CSR, formato eficiente para eigsh


# Diagonalizacion completa de una matriz hermitica: H psi_n = E_n psi_n
def diagonalizar_hamiltoniano(hamiltoniano):
    energias, vectores = la.eigh(hamiltoniano)                # Diagonaliza matriz hermitica y devuelve autovalores/autovectores
    orden = np.argsort(energias)                              # Ordena de menor a mayor energia
    return energias[orden], vectores[:, orden]


# Solo autovalores de una matriz hermitica: H psi_n = E_n psi_n
def diagonalizar_energias(hamiltoniano):
    return np.sort(la.eigvalsh(hamiltoniano))                 # Calcula y ordena solo autovalores, sin autovectores


# Diagonalizacion cerca de cero. Para N grande usa eigsh con shift-invert alrededor de sigma=0; para N pequeño usa diagonalizacion completa
def diagonalizar_cerca_cero_bdg(numero_sitios_usado, potencial_quimico_usado, energia_zeeman_usada, V_eff, Delta_ind, Z, hopping_usado, hopping_spin_orbita_usado, numero_estados=32):
    dimension = 4 * numero_sitios_usado                       # Dimension total BdG
    if scipy_sparse_disponible and dimension > 300:
        hamiltoniano = construir_hamiltoniano_bdg_realista_disperso(numero_sitios_usado, potencial_quimico_usado, energia_zeeman_usada, V_eff, Delta_ind, Z, hopping_usado, hopping_spin_orbita_usado)
        k = min(numero_estados, dimension - 2)                # Numero de estados a calcular
        try:
            energias, vectores = spla.eigsh(hamiltoniano, k=k, sigma=0.0, which="LM")  # Estados mas cercanos a E=0
        except Exception:
            energias, vectores = spla.eigsh(hamiltoniano, k=k, sigma=1.0e-7 * meV, which="LM")  # Pequeño desplazamiento si sigma=0 falla
        orden = np.argsort(energias)
        return energias[orden], vectores[:, orden], hamiltoniano
    if dimension > 2500:
        raise RuntimeError("La diagonalizacion densa seria demasiado grande. Para N alto necesitas SciPy con scipy.sparse.linalg.eigsh.")
    hamiltoniano = construir_hamiltoniano_bdg_realista(numero_sitios_usado, potencial_quimico_usado, energia_zeeman_usada, V_eff, Delta_ind, Z, hopping_usado, hopping_spin_orbita_usado)
    energias, vectores = diagonalizar_hamiltoniano(hamiltoniano)
    indices = np.argsort(np.abs(energias))[:numero_estados]   # Selecciona los estados mas cercanos a cero
    indices = indices[np.argsort(energias[indices])]          # Los reordena por energia
    return energias[indices], vectores[:, indices], hamiltoniano


# Error de hermiticidad: epsilon_H = ||H - H^\dagger|| / ||H||
def error_hermiticidad(hamiltoniano):
    if scipy_sparse_disponible and sp.issparse(hamiltoniano):
        norma = spla.norm(hamiltoniano)
        return float(spla.norm(hamiltoniano - hamiltoniano.getH()) / max(norma, 1.0e-30))
    return float(la.norm(hamiltoniano - hamiltoniano.conj().T) / max(la.norm(hamiltoniano), 1.0e-30))


# Error de simetria particula-hueco: epsilon_C = ||C H^* C^\dagger + H|| / ||H||
def error_simetria_particula_hueco(hamiltoniano, numero_sitios_usado):
    if scipy_sparse_disponible and sp.issparse(hamiltoniano):
        if hamiltoniano.shape[0] > 2500:
            return np.nan                                      # Evita convertir matrices enormes a densas
        hamiltoniano = hamiltoniano.toarray()
    operador_global = np.kron(np.eye(numero_sitios_usado, dtype=complex), operador_particula_hueco_local)
    with np.errstate(over="ignore", invalid="ignore", divide="ignore"):
        prueba = operador_global @ hamiltoniano.conj() @ operador_global.conj().T + hamiltoniano
    return float(la.norm(prueba) / max(la.norm(hamiltoniano), 1.0e-30))


# Error espectral de simetria BdG: epsilon_E = max_n |E_n + E_{-n}|
def error_espectro_particula_hueco(energias):
    energias_ordenadas = np.sort(np.asarray(energias, dtype=float))             # Ordena el espectro
    return float(np.max(np.abs(energias_ordenadas + energias_ordenadas[::-1]))) # Compara cada energia con su opuesta

###############################################################################
############################ Autoconsistencia Hartree #########################
###############################################################################

# Densidad normal usada en Hartree:
# n_i = sum_n f(E_n) ( |psi_{n,up}(i)|^2 + |psi_{n,down}(i)|^2 )
def densidad_normal_1d(energias, vectores, numero_sitios_usado, temperatura_usada=temperatura):
    densidad = np.zeros(numero_sitios_usado, dtype=float)                    # Inicializa la densidad n_i en todos los sitios
    ocupaciones = fermi(energias, temperatura_usada)                         # Calcula la ocupacion f(E_n) de cada autoestado normal
    for indice_estado, ocupacion in enumerate(ocupaciones):                  # Recorre todos los autoestados normales
        if ocupacion < 1.0e-12:                                              # Ignora estados practicamente vacios
            continue
        psi = vectores[:, indice_estado].reshape(numero_sitios_usado, 2)     # Reescribe el autovector como psi(i, spin)
        densidad += ocupacion * (np.abs(psi[:, 0])**2 + np.abs(psi[:, 1])**2) # Suma el peso ocupado de spin up y spin down
    return densidad                                                          # Devuelve la densidad normal n_i


# Densidad BdG usada solo como diagnostico/grafica:
# n_i^BdG = sum_{E_n>0} [ rho_u,n(i) f(E_n) + rho_v,n(i) (1 - f(E_n)) ]
# rho_u = |u_up|^2 + |u_down|^2, rho_v = |v_down|^2 + |-v_up|^2
def densidad_bdg_para_hartree(energias, vectores, numero_sitios_usado, temperatura_usada=temperatura):
    densidad = np.zeros(numero_sitios_usado, dtype=float)                    # Inicializa la densidad BdG en todos los sitios
    for indice_estado, energia in enumerate(energias):                       # Recorre los autoestados BdG calculados
        if energia <= 1.0e-12:                                               # Usa solo E>0 para no duplicar pares particula-hueco
            continue
        psi = vectores[:, indice_estado].reshape(numero_sitios_usado, 4)     # Reescribe el autovector como psi(i, componente BdG)
        rho_u = np.abs(psi[:, 0])**2 + np.abs(psi[:, 1])**2                  # Peso electronico local
        rho_v = np.abs(psi[:, 2])**2 + np.abs(psi[:, 3])**2                  # Peso de hueco local
        ocupacion = float(fermi(np.array([energia]), temperatura_usada)[0])  # Calcula f(E_n) para el estado positivo
        densidad += rho_u * ocupacion + rho_v * (1.0 - ocupacion)            # Suma la contribucion BdG a la densidad
    return densidad                                                          # Devuelve la densidad BdG


# Bucle Hartree normal:
# V_eff^(m)(x) = V_ext(x) + U_H^(m)(x)
# H_N[V_eff^(m)] psi_n = E_n psi_n
# U_H,calc = nu n
# U_H^(m+1) = (1 - eta) U_H^(m) + eta U_H,calc
def resolver_Hartree_normal(x, potencial_quimico_usado, energia_zeeman_usada, V_ext, hopping_usado, hopping_spin_orbita_usado, matriz_interaccion, temperatura_usada=temperatura, eta=mezcla_hartree, tolerancia=tolerancia_hartree, max_iteraciones=max_iteraciones_hartree):
    numero_sitios_usado = len(x)                                             # Numero de puntos de la malla 1D
    densidad_referencia = np.zeros(numero_sitios_usado, dtype=float)         # Referencia fija n_ref=0 en este modelo
    U_H = np.zeros(numero_sitios_usado, dtype=float)                         # Inicializa el potencial Hartree a cero
    historial_error = []                                                     # Guarda el error relativo de cada iteracion
    historial_cambio_maximo = []                                             # Guarda el cambio local maximo para la grafica de convergencia
    convergido = False                                                       # Indicador de convergencia
    densidad = np.zeros(numero_sitios_usado, dtype=float)                    # Inicializa la densidad

    for iteracion in range(max_iteraciones):                                 # Itera hasta converger o alcanzar el maximo
        V_eff = V_ext + U_H                                                  # Potencial total que entra en el Hamiltoniano normal
        hamiltoniano = construir_H_normal_1d(numero_sitios_usado, potencial_quimico_usado, energia_zeeman_usada, V_eff, hopping_usado, hopping_spin_orbita_usado) # Construye H_N
        energias, vectores = diagonalizar_hamiltoniano(hamiltoniano)         # Diagonaliza H_N
        densidad = densidad_normal_1d(energias, vectores, numero_sitios_usado, temperatura_usada) # Calcula n_i
        U_calculado = actualizar_hartree(densidad, densidad_referencia, matriz_interaccion) # Calcula U_H desde la densidad
        U_siguiente = (1.0 - eta) * U_H + eta * U_calculado                  # Mezcla el Hartree anterior con el nuevo
        cambio = U_siguiente - U_H                                           # Diferencia entre dos iteraciones consecutivas
        error = la.norm(cambio) / max(la.norm(U_siguiente), la.norm(U_H), 1.0e-12) # Error relativo de convergencia
        historial_error.append(float(error))                                 # Guarda el error relativo
        historial_cambio_maximo.append(float(np.max(np.abs(cambio))))        # Guarda el cambio maximo local
        U_H = U_siguiente                                                    # Actualiza Hartree

        if error < tolerancia:                                               # Si el cambio ya es suficientemente pequeño
            convergido = True                                                # Marca la convergencia
            break                                                            # Sale del bucle

    return {"energias": energias, "vectores": vectores, "hamiltoniano": hamiltoniano, "U_H": U_H, "V_eff": V_ext + U_H, "densidad": densidad, "historial_error": np.array(historial_error), "historial_cambio_maximo": np.array(historial_cambio_maximo), "convergido": convergido, "iteraciones": iteracion + 1}


###############################################################################
############################ Flujo del dispositivo realista ###################
###############################################################################

# Construye la base electrostatica y de proximidad: malla -> V_ext(x) -> Hartree normal -> V_eff(x) -> gamma(x), Delta_ind(x), Z(x)
def resolver_base_electrostatica(numero_sitios_usado=numero_sitios, V_gate_global_usado=None, configuracion=None, energia_zeeman_Hartree=0.0 * meV):
    configuracion = obtener_configuracion_modelo("realista_suave" if configuracion is None else configuracion) # Carga la configuracion fisica elegida
    if V_gate_global_usado is None: V_gate_global_usado = configuracion.get("V_gate_global", V_gate_global) # Si no se pasa gate, usa la de la configuracion o la global
    x, paso = crear_malla_1d(longitud_hilo, numero_sitios_usado) # Crea la malla 1D del nanohilo
    hopping_usado, hopping_so_usado = calcular_hoppings_discretos(alpha_rashba, paso) # Calcula t y t_so para ese paso espacial
    x_inicio_SC = configuracion.get("x_inicio_superconductor", x_inicio_superconductor) # Inicio de la zona cubierta por SC
    suavidad_SC = configuracion.get("suavidad_superconductor", suavidad_superconductor) # Suavidad de la transicion normal-SC
    V_ext, VSC_x, barrera, fSC = construir_potencial_externo_1d(x, Vg=V_gate_global_usado, amplitud_VSC=configuracion.get("amplitud_VSC_efectiva", amplitud_VSC_efectiva), Vb=configuracion.get("altura_barrera", altura_barrera), x_inicio=x_inicio_SC, suavidad=suavidad_SC) # Construye V_ext = E0 + Vg + VSC(x) + barrera
    U0_usado = configuracion.get("U0_hartree", U0_hartree) # Intensidad Hartree de la configuracion
    if abs(U0_usado) < 1.0e-30: # Si no hay Hartree, se evita resolver el bucle autoconsistente
        matriz_interaccion = None # No se necesita matriz de interaccion
        resultado_hartree = {"energias_normal": np.array([]), "vectores_normal": np.array([[]]), "U_H": np.zeros(numero_sitios_usado), "V_eff": V_ext.copy(), "densidad": np.zeros(numero_sitios_usado), "historial_error": np.array([0.0]), "historial_cambio_maximo": np.array([0.0]), "convergido": True, "iteraciones": 1} # En el caso sin Hartree, V_eff = V_ext
    else:
        matriz_interaccion = construir_matriz_interaccion_apantallada(x, U0_usado, lambda_apantallamiento) # Construye nu_ij
        resultado_hartree = resolver_Hartree_normal(x, configuracion.get("potencial_quimico", potencial_quimico), energia_zeeman_Hartree, V_ext, hopping_usado, hopping_so_usado, matriz_interaccion) # Resuelve Hartree normal
    gamma, Delta_ind, Z, fSC_prox = calcular_parametros_proximidad(x, gamma0=configuracion.get("gamma0_interfaz", gamma0_interfaz), Delta0=configuracion.get("Delta0_padre", Delta0_padre), x_inicio=x_inicio_SC, suavidad=suavidad_SC) # Calcula gamma(x), Delta_ind(x) y Z(x)
    return {"x": x, "paso": paso, "numero_sitios": numero_sitios_usado, "hopping": hopping_usado, "hopping_spin_orbita": hopping_so_usado, "V_ext": V_ext, "VSC_x": VSC_x, "barrera": barrera, "fSC": fSC, "matriz_interaccion": matriz_interaccion, "gamma": gamma, "Delta_ind": Delta_ind, "Z": Z, "resultado_hartree": resultado_hartree, "V_eff": resultado_hartree["V_eff"], "U_H": resultado_hartree["U_H"], "densidad": resultado_hartree["densidad"], "V_gate_global": V_gate_global_usado, "configuracion": configuracion, "potencial_quimico": configuracion.get("potencial_quimico", potencial_quimico), "zeeman_topologico": configuracion.get("zeeman_topologico", zeeman_topologico), "zeeman_transicion_estimado": estimar_zeeman_transicion_configuracion(configuracion)}

# Refinamiento electrostatico: calcula Hartree en una malla manejable e interpola U_H(x) a una malla fina para LDOS/modos
def resolver_base_electrostatica_refinada(numero_sitios_fino=numero_sitios_estados, numero_sitios_hartree=numero_sitios_barridos, V_gate_global_usado=None, configuracion=None, energia_zeeman_Hartree=0.0 * meV):
    configuracion = obtener_configuracion_modelo("realista_suave" if configuracion is None else configuracion) # Carga la configuracion
    base_gruesa = resolver_base_electrostatica(numero_sitios_usado=numero_sitios_hartree, V_gate_global_usado=V_gate_global_usado, configuracion=configuracion, energia_zeeman_Hartree=energia_zeeman_Hartree) # Resuelve Hartree en malla manejable
    configuracion_fina = configuracion.copy() # Copia la configuracion original
    configuracion_fina["U0_hartree"] = 0.0 # Evita recalcular Hartree en la malla fina
    base_fina = resolver_base_electrostatica(numero_sitios_usado=numero_sitios_fino, V_gate_global_usado=base_gruesa["V_gate_global"], configuracion=configuracion_fina, energia_zeeman_Hartree=energia_zeeman_Hartree) # Construye base fina sin Hartree propio
    xg = base_gruesa["x"] # Malla gruesa
    xf = base_fina["x"] # Malla fina
    U_fino = np.interp(xf, xg, base_gruesa["U_H"]) # Interpola U_H a la malla fina
    densidad_fina = np.interp(xf, xg, base_gruesa["densidad"]) # Interpola la densidad a la malla fina
    base_fina["U_H"] = U_fino # Sustituye Hartree fino por el interpolado
    base_fina["V_eff"] = base_fina["V_ext"] + U_fino # Reconstruye V_eff fino
    base_fina["densidad"] = densidad_fina # Guarda densidad interpolada
    base_fina["resultado_hartree"] = base_gruesa["resultado_hartree"] # Mantiene informacion de convergencia real
    base_fina["matriz_interaccion"] = None # No se usa Hartree directo en malla fina
    base_fina["configuracion"] = configuracion # Restaura la configuracion original
    return base_fina # Devuelve la base fina lista para BdG

# Diagonalizacion BdG sobre una base electrostatica ya convergida: H_BdG[V_eff, Delta_ind, Z] Psi_n = E_n Psi_n
def resolver_bdg_desde_base(base, energia_zeeman_usada, numero_estados_cerca_cero=48):
    N = base["numero_sitios"] # Numero de sitios de la base
    mu_usado = base.get("potencial_quimico", potencial_quimico) # Potencial quimico efectivo
    energias, vectores, hamiltoniano = diagonalizar_cerca_cero_bdg(N, mu_usado, energia_zeeman_usada, base["V_eff"], base["Delta_ind"], base["Z"], base["hopping"], base["hopping_spin_orbita"], numero_estados=numero_estados_cerca_cero) # Calcula estados BdG cercanos a cero
    densidad_bdg = densidad_bdg_para_hartree(energias, vectores, N, temperatura) # Densidad BdG subgap para diagnostico/grafica
    return {"energias": energias, "vectores": vectores, "hamiltoniano": hamiltoniano, "densidad_hartree_bdg": densidad_bdg, "densidad_bdg_es_subespacio": len(energias) < 4 * N, "V_eff": base["V_eff"], "U_H": base["U_H"], "convergido": base["resultado_hartree"]["convergido"], "iteraciones": base["resultado_hartree"]["iteraciones"]}

# Calcula los tres regímenes representativos: trivial, transición y topológico
# H_BdG(Gamma) Psi_n = E_n Psi_n
def calcular_resultados_regimenes(base, numero_estados_cerca_cero=48):
    casos = obtener_zeeman_regimenes(base)
    resultados = {}
    for nombre, energia_zeeman in casos.items():
        resultados[nombre] = resolver_bdg_desde_base(base, energia_zeeman, numero_estados_cerca_cero=numero_estados_cerca_cero)
    return resultados


# Calcula cada régimen con su propia base Hartree autoconsistente recalculada a ese Gamma
def calcular_resultados_regimenes_desde_bases(bases_regimenes, numero_estados_cerca_cero=48):
    resultados = {}
    for nombre in ["trivial", "transicion", "topologico"]:
        base = bases_regimenes[nombre]
        energia_zeeman = obtener_zeeman_regimenes(base)[nombre]
        resultados[nombre] = resolver_bdg_desde_base(base, energia_zeeman, numero_estados_cerca_cero=numero_estados_cerca_cero)
    return resultados

# Devuelve la base correspondiente al régimen sin cambiar los perfiles físicos
def seleccionar_base_regimen(base_o_bases, clave):
    if isinstance(base_o_bases, dict) and clave in base_o_bases and isinstance(base_o_bases[clave], dict) and "x" in base_o_bases[clave]:
        return base_o_bases[clave]
    return base_o_bases

###############################################################################
############################# Observables fisicos #############################
###############################################################################

# Separa un autovector BdG en componentes u_up, u_down, v_down y -v_up
def separar_componentes_bdg(vector_estado, numero_sitios_usado):
    psi = vector_estado.reshape(numero_sitios_usado, 4)
    return psi[:, 0].copy(), psi[:, 1].copy(), psi[:, 2].copy(), psi[:, 3].copy()

# Calcula la densidad total normalizada de un estado BdG
def densidad_estado(vector_estado, numero_sitios_usado):
    u_arriba, u_abajo, v_abajo, menos_v_arriba = separar_componentes_bdg(vector_estado, numero_sitios_usado)
    densidad = np.abs(u_arriba)**2 + np.abs(u_abajo)**2 + np.abs(v_abajo)**2 + np.abs(menos_v_arriba)**2
    return densidad / max(float(np.sum(densidad)), 1.0e-30)

# Devuelve por separado densidad electronica, de hueco y total
def densidad_electron_hueco_total(vector_estado, numero_sitios_usado):
    u_arriba, u_abajo, v_abajo, menos_v_arriba = separar_componentes_bdg(vector_estado, numero_sitios_usado)
    densidad_electron = np.abs(u_arriba)**2 + np.abs(u_abajo)**2
    densidad_hueco = np.abs(v_abajo)**2 + np.abs(menos_v_arriba)**2
    norma = max(float(np.sum(densidad_electron + densidad_hueco)), 1.0e-30)
    return densidad_electron / norma, densidad_hueco / norma, (densidad_electron + densidad_hueco) / norma

# Encuentra el indice del estado mas cercano a cero
def indice_estado_mas_cercano_a_cero(energias):
    return int(np.argmin(np.abs(energias)))

# Encuentra el indice del estado positivo mas cercano a cero
def indice_estado_positivo_mas_cercano_a_cero(energias):
    indices_positivos = np.where(energias >= 0.0)[0]
    if len(indices_positivos) == 0:
        return int(np.argmin(np.abs(energias)))
    return int(indices_positivos[np.argmin(np.abs(energias[indices_positivos]))])

# Devuelve min |E|
def energia_minima_absoluta(energias):
    return float(np.min(np.abs(energias)))

# Calcula el peso de un estado en ambos bordes
def pesos_bordes(vector_estado, x, longitud_borde=150.0 * nm):
    densidad = densidad_estado(vector_estado, len(x))
    mascara_izquierda = x <= longitud_borde
    mascara_derecha = x >= x[-1] - longitud_borde
    return float(np.sum(densidad[mascara_izquierda])), float(np.sum(densidad[mascara_derecha]))

# Aplica la simetria particula-hueco
def transformar_particula_hueco(vector_estado, numero_sitios_usado):
    psi = vector_estado.reshape(numero_sitios_usado, 4)
    psi_ph = np.zeros_like(psi)
    for i in range(numero_sitios_usado):
        psi_ph[i, :] = operador_particula_hueco_local @ np.conj(psi[i, :])
    return psi_ph.reshape(4 * numero_sitios_usado)

# Elimina la fase global entre dos estados
def alinear_fase(estado_referencia, estado_objetivo):
    solapamiento = np.vdot(estado_referencia, estado_objetivo)
    if abs(solapamiento) < 1.0e-30:
        return estado_objetivo.copy()
    return estado_objetivo * np.exp(-1.0j * np.angle(solapamiento))

# Construye dos modos localizados izquierda/derecha desde el subespacio casi cero
def construir_majoranas_desde_pareja(energias, vectores, numero_sitios_usado):
    indices_cercanos = np.argsort(np.abs(energias))[:2]                    # Selecciona la pareja casi cero
    indices_cercanos = indices_cercanos[np.argsort(energias[indices_cercanos])]
    estados = vectores[:, indices_cercanos].copy()                         # Subespacio de dos estados casi degenerados
    mascara_izquierda = np.zeros(4 * numero_sitios_usado, dtype=float)
    mascara_izquierda[:4 * (numero_sitios_usado // 2)] = 1.0               # Operador de peso en la mitad izquierda
    matriz_izquierda = np.empty((2, 2), dtype=complex)
    for a in range(2):
        for b in range(2):
            matriz_izquierda[a, b] = np.vdot(estados[:, a], mascara_izquierda * estados[:, b])
    valores, coeficientes = la.eigh(matriz_izquierda)                       # Diagonaliza el peso izquierdo dentro del subespacio
    estado_derecha = estados @ coeficientes[:, 0]
    estado_izquierda = estados @ coeficientes[:, 1]
    densidad_izquierda = densidad_estado(estado_izquierda, numero_sitios_usado)
    densidad_derecha = densidad_estado(estado_derecha, numero_sitios_usado)
    mitad = numero_sitios_usado // 2
    if np.sum(densidad_izquierda[:mitad]) < np.sum(densidad_derecha[:mitad]):
        densidad_izquierda, densidad_derecha = densidad_derecha, densidad_izquierda
    energia_1 = float(energias[indices_cercanos[0]])
    energia_2 = float(energias[indices_cercanos[1]])
    return densidad_izquierda, densidad_derecha, energia_1, energia_2

# Lorentziana para suavizar la LDOS
def lorentziana(energia, centro, anchura):
    return anchura / (np.pi * ((energia - centro)**2 + anchura**2))



# Calcula una LDOS BdG simetrizada: rho(x,E)=sum_{E_n>0} rho_n(x)[L(E-E_n)+L(E+E_n)]
def calcular_ldos(numero_sitios_usado, energias, vectores, energias_ldos, anchura):
    ldos = np.zeros((len(energias_ldos), numero_sitios_usado), dtype=float)
    for indice_estado, energia_n in enumerate(energias):
        if energia_n <= 1.0e-12:
            continue
        densidad_total = densidad_estado(vectores[:, indice_estado], numero_sitios_usado)
        peso_positivo = lorentziana(energias_ldos, energia_n, anchura)
        peso_negativo = lorentziana(energias_ldos, -energia_n, anchura)
        ldos += (peso_positivo[:, None] + peso_negativo[:, None]) * densidad_total[None, :]
    return ldos


###############################################################################
############################## Barridos fisicos ###############################
###############################################################################

# Barre Gamma y calcula el espectro finito

def barrer_zeeman_desde_base(base, valores_zeeman, recalcular_hartree_por_zeeman=False):
    espectros = []
    energias_minimas = []
    pesos_izquierdos = []
    pesos_derechos = []
    for energia_zeeman in valores_zeeman:
        base_usada = resolver_base_electrostatica(numero_sitios_usado=base["numero_sitios"], V_gate_global_usado=base["V_gate_global"], configuracion=base["configuracion"], energia_zeeman_Hartree=energia_zeeman) if recalcular_hartree_por_zeeman else base        
        resultado = resolver_bdg_desde_base(base_usada, energia_zeeman)
        energias = resultado["energias"]
        indice = indice_estado_mas_cercano_a_cero(energias)
        peso_izq, peso_der = pesos_bordes(resultado["vectores"][:, indice], base_usada["x"])
        espectros.append(energias)
        energias_minimas.append(energia_minima_absoluta(energias))
        pesos_izquierdos.append(peso_izq)
        pesos_derechos.append(peso_der)
    return {"valores_zeeman": np.asarray(valores_zeeman), "espectros": np.asarray(espectros), "energias_minimas": np.asarray(energias_minimas), "pesos_izquierdos": np.asarray(pesos_izquierdos), "pesos_derechos": np.asarray(pesos_derechos)}

# Barre mu y calcula el espectro finito para un Gamma fijo
def barrer_mu_desde_base(base, valores_mu, energia_zeeman_usada, numero_estados_cerca_cero=48):
    espectros = []
    energias_minimas = []
    for mu_usado in valores_mu:
        energias, vectores, hamiltoniano = diagonalizar_cerca_cero_bdg(base["numero_sitios"], mu_usado, energia_zeeman_usada, base["V_eff"], base["Delta_ind"], base["Z"], base["hopping"], base["hopping_spin_orbita"], numero_estados=numero_estados_cerca_cero)
        espectros.append(energias)
        energias_minimas.append(energia_minima_absoluta(energias))
    return {"valores_mu": np.asarray(valores_mu), "espectros": np.asarray(espectros), "energias_minimas": np.asarray(energias_minimas), "energia_zeeman": energia_zeeman_usada}



# Usa N alto solo si hay diagonalizacion dispersa; en caso contrario mantiene un modo seguro
def elegir_numero_sitios_para_estados():
    return numero_sitios_estados if scipy_sparse_disponible else numero_sitios_prueba


# Calcula bases autoconsistentes separadas para ver como cambia V_eff con Gamma
def calcular_bases_autoconsistentes_por_regimen(configuracion="realista_suave", numero_sitios_fino=numero_sitios_perfiles, numero_sitios_hartree=numero_sitios_barridos):
    config = obtener_configuracion_modelo(configuracion)
    base_referencia = resolver_base_electrostatica(numero_sitios_usado=numero_sitios_hartree, configuracion=config)
    zeeman_regimenes = obtener_zeeman_regimenes(base_referencia)
    bases = {}
    for nombre, energia_zeeman in zeeman_regimenes.items():
        bases[nombre] = resolver_base_electrostatica_refinada(numero_sitios_fino=numero_sitios_fino, numero_sitios_hartree=numero_sitios_hartree, configuracion=config, energia_zeeman_Hartree=energia_zeeman)
    return bases

# Hamiltoniano bulk uniforme H(k) correspondiente al modelo 1D realista
def construir_hamiltoniano_bulk_k_realista(k, potencial_quimico_usado, energia_zeeman_usada, V_uniforme, Delta_uniforme, Z_uniforme, hopping_usado, hopping_spin_orbita_usado, paso_usado):
    fase = np.exp(1.0j * k * paso_usado)
    bloque_local = Z_uniforme * (2.0 * hopping_usado - potencial_quimico_usado + V_uniforme) * tau_z + Z_uniforme * energia_zeeman_usada * sigma_x + Delta_uniforme.real * tau_x - Delta_uniforme.imag * tau_y
    bloque_hopping = Z_uniforme * (-hopping_usado * tau_z + 1.0j * hopping_spin_orbita_usado * sigma_y_tau_z)
    hamiltoniano_k = bloque_local + bloque_hopping * fase + bloque_hopping.conj().T * np.conj(fase) # Sistema infinito: el hopping se vuelve dependiente de k
    return 0.5 * (hamiltoniano_k + hamiltoniano_k.conj().T)

# Calcula las bandas bulk de la region cubierta efectiva
def calcular_bandas_bulk_realistas(base, energia_zeeman_usada, numero_puntos_k=420, k_maximo=0.16 / nm):
    mascara = base["fSC"] > 0.90
    if not np.any(mascara):
        mascara = np.ones_like(base["x"], dtype=bool)
    V_uniforme = float(np.mean(base["V_eff"][mascara]))
    Delta_uniforme = complex(np.mean(base["Delta_ind"][mascara]))
    Z_uniforme = float(np.mean(base["Z"][mascara]))
    k_limite = min(np.pi / base["paso"], k_maximo)
    valores_k = np.linspace(-k_limite, k_limite, numero_puntos_k)
    bandas = np.zeros((numero_puntos_k, 4), dtype=float)
    for i, k in enumerate(valores_k):
        Hk = construir_hamiltoniano_bulk_k_realista(k, base.get("potencial_quimico", potencial_quimico), energia_zeeman_usada, V_uniforme, Delta_uniforme, Z_uniforme, base["hopping"], base["hopping_spin_orbita"], base["paso"])
        bandas[i, :] = la.eigvalsh(Hk)
    parametros_bulk = {"V_uniforme": V_uniforme, "Delta_uniforme": Delta_uniforme, "Z_uniforme": Z_uniforme}
    return valores_k, bandas, parametros_bulk



###############################################################################
################################ Diagnósticos #################################
###############################################################################

# Estima parámetros uniformes de la región cubierta del hilo
def parametros_bulk_desde_base(base):
    mascara = base["fSC"] > 0.90
    if not np.any(mascara):
        mascara = np.ones_like(base["x"], dtype=bool)
    V_uniforme = float(np.mean(base["V_eff"][mascara]))
    Delta_uniforme = float(np.mean(np.real(base["Delta_ind"][mascara])))
    Z_uniforme = float(np.mean(base["Z"][mascara]))
    mu_uniforme = float(base.get("potencial_quimico", potencial_quimico) - V_uniforme)
    return {"V_uniforme": V_uniforme, "Delta_uniforme": Delta_uniforme, "Z_uniforme": Z_uniforme, "mu_uniforme": mu_uniforme}

# Gamma_c = sqrt(mu_eff^2+(Delta_ind/Z)^2)
def estimar_gamma_critico_desde_base(base):
    p = parametros_bulk_desde_base(base)
    return float(np.sqrt(p["mu_uniforme"]**2 + (p["Delta_uniforme"] / max(p["Z_uniforme"], 1.0e-30))**2))

# Devuelve los tres valores de Gamma usados en las figuras finales
def obtener_zeeman_regimenes(base):
    gamma_c = estimar_gamma_critico_desde_base(base)
    gamma_top = base.get("zeeman_topologico", zeeman_topologico)
    if gamma_top <= 1.10 * gamma_c:
        gamma_top = 1.35 * gamma_c
    return {"trivial": zeeman_trivial, "transicion": gamma_c, "topologico": gamma_top}

# mu_c = +-sqrt(Gamma^2-(Delta_ind/Z)^2)
def fronteras_mu_criticas_desde_base(base, energia_zeeman_usada):
    p = parametros_bulk_desde_base(base)
    delta_sobre_Z = p["Delta_uniforme"] / max(p["Z_uniforme"], 1.0e-30)
    argumento = energia_zeeman_usada**2 - delta_sobre_Z**2
    if argumento <= 0.0:
        return None
    valor = float(np.sqrt(argumento))
    V_uniforme = p["V_uniforme"]
    return np.array([V_uniforme - valor, V_uniforme + valor])

# Error relativo de simetría de la LDOS: max|rho(E,x)-rho(-E,x)|/max|rho|
def error_simetria_ldos(ldos):
    return float(np.max(np.abs(ldos - ldos[::-1, :])) / max(np.max(np.abs(ldos)), 1.0e-30))

# Diagnóstico de localización en bordes frente al centro
def diagnostico_localizacion_vector(vector_estado, x, etiqueta, longitud_borde=150.0 * nm):
    densidad = densidad_estado(vector_estado, len(x))
    mascara_izquierda = x <= longitud_borde
    mascara_derecha = x >= x[-1] - longitud_borde
    mascara_centro = (x >= 0.40 * x[-1]) & (x <= 0.60 * x[-1])
    peso_izquierda = float(np.sum(densidad[mascara_izquierda]))
    peso_derecha = float(np.sum(densidad[mascara_derecha]))
    peso_bordes = peso_izquierda + peso_derecha
    peso_centro = float(np.sum(densidad[mascara_centro]))
    print(f"Localización {etiqueta}: borde izquierdo={peso_izquierda:.3f}, borde derecho={peso_derecha:.3f}, bordes={peso_bordes:.3f}, centro={peso_centro:.3f}")
    if etiqueta == "topológico" and (peso_bordes < 0.45 or peso_bordes < 2.0 * max(peso_centro, 1.0e-30)):
        print("ADVERTENCIA: el estado cercano a cero topológico no está claramente dominado por los bordes; revisar parámetros, tamaño del hilo, barrera o construcción de la LDOS")
    return {"izquierda": peso_izquierda, "derecha": peso_derecha, "bordes": peso_bordes, "centro": peso_centro}

# Busca la pareja +-E más cercana a cero y devuelve las dos energías
def pareja_mas_cercana_cero(energias):
    energias = np.asarray(energias, dtype=float)
    indices_pos = np.where(energias >= 0.0)[0]
    indices_neg = np.where(energias < 0.0)[0]
    if len(indices_pos) == 0 or len(indices_neg) == 0:
        indice = int(np.argmin(np.abs(energias)))
        return energias[indice], np.nan
    ip = indices_pos[np.argmin(np.abs(energias[indices_pos]))]
    ineg = indices_neg[np.argmin(np.abs(energias[indices_neg] + energias[ip]))]
    return float(energias[ineg]), float(energias[ip])

# Imprime comprobaciones después de obtener los resultados principales
def ejecutar_diagnosticos_finales(base_estados, resultados_regimenes, datos_zeeman, ldos_por_regimen, residuos_bandas):
    print("\nComprobaciones finales:")
    for clave, etiqueta in [("trivial", "trivial"), ("transicion", "transición"), ("topologico", "topológico")]:
        energias = resultados_regimenes[clave]["energias"]
        err = error_espectro_particula_hueco(energias)
        print(f"  Simetría espectral {etiqueta}: {err:.3e} meV")
        if err > 1.0e-6:
            print("  ADVERTENCIA: el subespacio calculado no conserva bien pares ±E; aumentar numero_estados_cerca_cero o revisar eigsh")
    energias_top = np.sort(np.abs(resultados_regimenes["topologico"]["energias"]))
    e0 = float(energias_top[0])
    gap_excitado = float(energias_top[2]) if len(energias_top) > 2 else float(energias_top[-1])
    print(f"  Cero topológico: |E0|={e0:.3e} meV, gap excitado≈{gap_excitado:.3e} meV")
    if e0 > 0.25 * max(gap_excitado, 1.0e-30):
        print("  ADVERTENCIA: el estado más cercano a cero no está suficientemente separado del gap excitado")
    gamma_c = estimar_gamma_critico_desde_base(base_estados)
    valores = datos_zeeman["valores_zeeman"]
    minimos = datos_zeeman["energias_minimas"]
    indice_minimo = int(np.argmin(minimos))
    print(f"  Gamma_c teórico={gamma_c / meV:.4f} meV, cierre numérico mínimo en Gamma={valores[indice_minimo] / meV:.4f} meV con min|E|={minimos[indice_minimo] / meV:.3e} meV")
    if abs(valores[indice_minimo] - gamma_c) > 0.18 * max(gamma_c, 1.0e-30):
        print("  ADVERTENCIA: el cierre numérico no cae cerca de Gamma_c; revisar mu efectivo, Z, Delta_ind, V_eff o resolución del barrido")
    for clave, ldos in ldos_por_regimen.items():
        err_ldos = error_simetria_ldos(ldos)
        print(f"  Simetría LDOS {clave}: {err_ldos:.3e}")
        if err_ldos > 1.0e-8:
            print("  ADVERTENCIA: la LDOS no es simétrica respecto a E=0 dentro de la tolerancia")
    indice_top = indice_estado_mas_cercano_a_cero(resultados_regimenes["topologico"]["energias"])
    diagnostico_localizacion_vector(resultados_regimenes["topologico"]["vectores"][:, indice_top], base_estados["x"], "topológico")
    if residuos_bandas is not None:
        print(f"  Compatibilidad bandas analíticas/diag. H(k): residuo máximo={max(residuos_bandas.values()):.3e} meV")
        if max(residuos_bandas.values()) > 0.08:
            print("  ADVERTENCIA: las bandas analíticas y la diagonalización de H(k) no coinciden; revisar Hamiltoniano bulk o fórmula analítica")
    else:
        print("  Bandas bulk H(k): comprobación omitida porque ahora se grafican autovalores finitos con bordes abiertos")

###############################################################################
################################## Gráficas ###################################
###############################################################################

# Estilo común de las figuras
def preparar_estilo_graficas():
    plt.style.use("default")
    plt.rcParams.update({
        "figure.figsize": (9.4, 5.3),
        "axes.grid": True,
        "grid.alpha": 0.23,
        "font.size": 11,
        "axes.titlesize": 11,
        "axes.labelsize": 12,
        "legend.fontsize": 8.5,
        "lines.linewidth": 2.0,
        "figure.dpi": 120,
        "figure.titlesize": 13
    })

# Coloca cada autovalor finito en el k donde la banda bulk tiene energía más cercana
def puntos_finitos_sobre_bandas(valores_k, bandas, energias_finitas):
    k_puntos = []
    e_puntos = []
    residuos = []
    for energia in energias_finitas:
        diferencia = np.abs(bandas - energia)
        indice_k, indice_banda = np.unravel_index(np.argmin(diferencia), diferencia.shape)
        k_puntos.append(valores_k[indice_k])
        e_puntos.append(energia)
        residuos.append(float(diferencia[indice_k, indice_banda]))
    return np.asarray(k_puntos), np.asarray(e_puntos), np.asarray(residuos)

# Bandas bulk analíticas y diagonalización directa de H(k) para trivial/transición/topológico
def bandas_bulk_realistas_analiticas_q(base, energia_zeeman_usada, puntos_q):
    mascara = base["fSC"] > 0.90
    if not np.any(mascara):
        mascara = np.ones_like(base["x"], dtype=bool)
    V_uniforme = float(np.mean(base["V_eff"][mascara]))
    Delta_uniforme = float(np.mean(np.real(base["Delta_ind"][mascara])))
    Z_uniforme = float(np.mean(base["Z"][mascara]))
    mu_usado = float(base.get("potencial_quimico", potencial_quimico))
    bandas = np.zeros((len(puntos_q), 4), dtype=float)
    for i, q in enumerate(puntos_q):
        xi = Z_uniforme * (2.0 * base["hopping"] - mu_usado + V_uniforme - 2.0 * base["hopping"] * np.cos(q))
        rashba = -2.0 * Z_uniforme * base["hopping_spin_orbita"] * np.sin(q)
        gamma = Z_uniforme * energia_zeeman_usada
        Delta_abs = abs(Delta_uniforme)
        A = xi**2 + rashba**2 + gamma**2 + Delta_abs**2
        B = 2.0 * np.sqrt(gamma**2 * (xi**2 + Delta_abs**2) + xi**2 * rashba**2)
        energia_menor = np.sqrt(max(A - B, 0.0))
        energia_mayor = np.sqrt(max(A + B, 0.0))
        bandas[i, :] = [-energia_mayor, -energia_menor, energia_menor, energia_mayor]
    parametros = {"V_uniforme": V_uniforme, "Delta_uniforme": Delta_uniforme, "Z_uniforme": Z_uniforme}
    return bandas, parametros

# Diagonalización directa de H(k) usando el mismo Hamiltoniano bulk efectivo
def bandas_bulk_realistas_diagonalizando_q(base, energia_zeeman_usada, puntos_q):
    mascara = base["fSC"] > 0.90
    if not np.any(mascara):
        mascara = np.ones_like(base["x"], dtype=bool)
    V_uniforme = float(np.mean(base["V_eff"][mascara]))
    Delta_uniforme = complex(np.mean(base["Delta_ind"][mascara]))
    Z_uniforme = float(np.mean(base["Z"][mascara]))
    mu_usado = float(base.get("potencial_quimico", potencial_quimico))
    bandas = np.zeros((len(puntos_q), 4), dtype=float)
    for i, q in enumerate(puntos_q):
        bloque_local = Z_uniforme * (2.0 * base["hopping"] - mu_usado + V_uniforme) * tau_z + Z_uniforme * energia_zeeman_usada * sigma_x + Delta_uniforme.real * tau_x - Delta_uniforme.imag * tau_y
        bloque_hopping = Z_uniforme * (-base["hopping"] * tau_z + 1.0j * base["hopping_spin_orbita"] * sigma_y_tau_z)
        Hq = bloque_local + bloque_hopping * np.exp(1.0j * q) + bloque_hopping.conj().T * np.exp(-1.0j * q)
        bandas[i, :] = la.eigvalsh(0.5 * (Hq + Hq.conj().T))
    return bandas

# Bandas bulk analíticas + cruces de diagonalización de H(k) para trivial/transición/topológico
def bandas_bulk_realistas_analiticas_k(base, energia_zeeman_usada, valores_k):
    mascara = base["fSC"] > 0.90
    if not np.any(mascara):
        mascara = np.ones_like(base["x"], dtype=bool)
    V_uniforme = float(np.mean(base["V_eff"][mascara]))
    Delta_uniforme = float(np.mean(np.real(base["Delta_ind"][mascara])))
    Z_uniforme = float(np.mean(base["Z"][mascara]))
    mu_usado = float(base.get("potencial_quimico", potencial_quimico))
    bandas = np.zeros((len(valores_k), 4), dtype=float)
    for i, k in enumerate(valores_k):
        q = k * base["paso"]
        xi = Z_uniforme * (2.0 * base["hopping"] - mu_usado + V_uniforme - 2.0 * base["hopping"] * np.cos(q))
        rashba = -2.0 * Z_uniforme * base["hopping_spin_orbita"] * np.sin(q)
        gamma = Z_uniforme * energia_zeeman_usada
        Delta_abs = abs(Delta_uniforme)
        A = xi**2 + rashba**2 + gamma**2 + Delta_abs**2
        B = 2.0 * np.sqrt(gamma**2 * (xi**2 + Delta_abs**2) + xi**2 * rashba**2)
        energia_menor = np.sqrt(max(A - B, 0.0))
        energia_mayor = np.sqrt(max(A + B, 0.0))
        bandas[i, :] = [-energia_mayor, -energia_menor, energia_menor, energia_mayor]
    parametros = {"V_uniforme": V_uniforme, "Delta_uniforme": Delta_uniforme, "Z_uniforme": Z_uniforme}
    return bandas, parametros
"""
# Bandas bulk analíticas + cruces de diagonalización de H(k) para trivial/transición/topológico
def graficar_bandas_bulk_teoria_y_autovalores(base_o_bases, resultados_regimenes=None, numero_autovalores=18):
    claves = ["trivial", "transicion", "topologico"]
    etiquetas = ["trivial", "transición", "topológico"]
    figura, ejes = plt.subplots(1, 3, figsize=(13.8, 4.25), sharey=True)
    residuos_maximos = {}
    for eje, clave, etiqueta in zip(ejes, claves, etiquetas):
        base = seleccionar_base_regimen(base_o_bases, clave)
        energia_zeeman = obtener_zeeman_regimenes(base)[clave]
        valores_k, bandas_diag_continuas, parametros_bulk = calcular_bandas_bulk_realistas(base, energia_zeeman, numero_puntos_k=420, k_maximo=0.16 / nm)
        bandas_teoricas, _ = bandas_bulk_realistas_analiticas_k(base, energia_zeeman, valores_k)
        paso_marca = max(len(valores_k) // 34, 1)
        k_marcas = valores_k[::paso_marca]
        bandas_diag_marcas = bandas_diag_continuas[::paso_marca, :]
        bandas_teoricas_marcas, _ = bandas_bulk_realistas_analiticas_k(base, energia_zeeman, k_marcas)
        residuos_maximos[clave] = float(np.max(np.abs(bandas_diag_marcas - bandas_teoricas_marcas)))
        for indice_banda in range(4):
            eje.plot(valores_k / (1.0 / nm), bandas_teoricas[:, indice_banda] / meV, color="0.25", linewidth=1.65, alpha=0.86)
            eje.scatter(k_marcas / (1.0 / nm), bandas_diag_marcas[:, indice_banda] / meV, s=18, marker="x", alpha=0.78)
        eje.axhline(0.0, color="black", linewidth=1.0, alpha=0.55)
        eje.set_title(f"{etiqueta}\n$\\Gamma={energia_zeeman / meV:.2f}$ meV", fontsize=10.5)
        eje.set_xlabel(r"$k$ (nm$^{-1}$)")
        eje.set_ylim(-0.85, 0.85)
        eje.text(0.03, 0.94, f"$\\Delta={parametros_bulk['Delta_uniforme'].real / meV:.3f}$ meV, $Z={parametros_bulk['Z_uniforme']:.2f}$", transform=eje.transAxes, fontsize=8.0, va="top")
        if residuos_maximos[clave] > 1.0e-7:
            eje.text(0.03, 0.08, "ADVERTENCIA\nbulk no coincide", transform=eje.transAxes, fontsize=8.0, color="crimson", va="bottom")
    ejes[0].set_ylabel("Energía BdG (meV)")
    leyenda = [
        Line2D([0], [0], color="0.25", linewidth=1.65, label="expresión analítica bulk"),
        Line2D([0], [0], marker="x", linestyle="", markersize=6, label="diag. de $H(k)$")
    ]
    ejes[1].legend(handles=leyenda, loc="upper center")
    figura.suptitle("Bandas bulk del nanohilo BdG", fontsize=13)
    figura.tight_layout()
    print("Compatibilidad bandas bulk analíticas vs diagonalización H(k):")
    for clave in claves:
        print(f"  {clave}: residuo={residuos_maximos[clave]:.3e} meV")
        if residuos_maximos[clave] > 1.0e-7:
            print("  ADVERTENCIA: las bandas analíticas y la diagonalización de H(k) no coinciden")
    plt.show()
    return figura, residuos_maximos
"""

# Autovalores finitos del Hamiltoniano BdG realista en los tres regímenes
def graficar_autovalores_finitos_subgap_tres_regimenes(
    resultados_regimenes,
    ventana_energia=0.60 * meV,
    numero_autovalores=80,
    tamano_punto=18.0
):
    claves = ["trivial", "transicion", "topologico"]
    etiquetas = ["trivial", "transición", "topológico"]

    figura, ejes = plt.subplots(1, 3, figsize=(13.8, 4.25), sharey=True)

    for eje, clave, etiqueta in zip(ejes, claves, etiquetas):
        resultado = resultados_regimenes[clave]
        energias = np.asarray(resultado["energias"], dtype=float)

        indices_ordenados = np.argsort(np.abs(energias))
        indices_mostrar = indices_ordenados[:numero_autovalores]
        energias_mostrar = energias[indices_mostrar]

        orden_por_energia = np.argsort(energias_mostrar)
        energias_mostrar = energias_mostrar[orden_por_energia]

        niveles = np.arange(len(energias_mostrar))

        eje.scatter(
            niveles,
            energias_mostrar / meV,
            s=tamano_punto,
            color="black",
            alpha=0.78,
            linewidths=0.0
        )

        energia_minima = float(np.min(np.abs(energias)))
        eje.axhline(0.0, color="black", linewidth=1.0, alpha=0.60)

        eje.set_ylim(-ventana_energia / meV, ventana_energia / meV)
        eje.set_title(
            f"{etiqueta}\n"
            rf"$\min |E|={energia_minima / meV:.3e}$ meV",
            fontsize=10.0
        )
        eje.set_xlabel("Índice de nivel ordenado")

        if clave == "topologico":
            eje.text(
                0.04,
                0.92,
                "modo casi cero\nfinito",
                transform=eje.transAxes,
                fontsize=8.5,
                color="crimson",
                va="top"
            )

    ejes[0].set_ylabel("Energía BdG finita (meV)")

    figura.suptitle(
        "Autovalores finitos del Hamiltoniano BdG realista",
        fontsize=13
    )
    figura.tight_layout()

    print("Autovalores finitos subgap:")
    for clave in claves:
        energias = np.asarray(resultados_regimenes[clave]["energias"], dtype=float)
        print(
            f"  {clave}: min|E| = {np.min(np.abs(energias)) / meV:.6e} meV"
        )


    plt.show()
    return figura


# Selecciona las ramas BdG positivas y negativas del sistema finito
# Se usan para construir una representación visual tipo bandas del sistema real-space con bordes abiertos

def obtener_ramas_bdg_finitas(energias, tolerancia=1.0e-12 * meV):
    energias = np.sort(np.asarray(energias, dtype=float))
    rama_hueco = energias[energias < -tolerancia]
    rama_cuasiparticula = energias[energias > tolerancia]
    if len(rama_cuasiparticula) == 0:
        rama_cuasiparticula = energias[energias >= 0.0]
    if len(rama_hueco) == 0:
        rama_hueco = energias[energias <= 0.0]
    rama_cuasiparticula = np.sort(rama_cuasiparticula)
    rama_hueco = np.sort(rama_hueco)
    return rama_hueco, rama_cuasiparticula

# Convierte autovalores finitos en una representación tipo bandas usando un pseudo-k
# Se separan la rama cuasipartícula E>0 y la rama cuasihueco E<0 para mantener el formato visual anterior

def preparar_puntos_tipo_bandas_finitas_eh(energias, k_maximo, numero_niveles=None):
    rama_hueco, rama_cuasiparticula = obtener_ramas_bdg_finitas(energias)
    if numero_niveles is not None:
        cantidad_hueco = min(numero_niveles, len(rama_hueco))
        cantidad_particula = min(numero_niveles, len(rama_cuasiparticula))
        rama_hueco = rama_hueco[-cantidad_hueco:]
        rama_cuasiparticula = rama_cuasiparticula[:cantidad_particula]

    def construir_rama(energias_rama):
        if len(energias_rama) == 0:
            return np.array([]), np.array([])
        k_positivo = np.linspace(0.0, k_maximo, len(energias_rama))
        k_pseudo = np.concatenate((-k_positivo[:0:-1], k_positivo))
        energias_pseudo = np.concatenate((energias_rama[:0:-1], energias_rama))
        return k_pseudo, energias_pseudo

    k_h, e_h = construir_rama(rama_hueco)
    k_p, e_p = construir_rama(rama_cuasiparticula)
    return {
        "k_hueco": k_h,
        "E_hueco": e_h,
        "k_cuasiparticula": k_p,
        "E_cuasiparticula": e_p,
        "rama_hueco": rama_hueco,
        "rama_cuasiparticula": rama_cuasiparticula,
    }

# Autovalores finitos del Hamiltoniano BdG realista completo
# Esta función diagonaliza el sistema real-space con bordes abiertos, no H(k)

def calcular_espectro_finito_completo_desde_base(base, clave_regimen, dimension_maxima_densa=2600):
    N = base["numero_sitios"]
    dimension = 4 * N
    energia_zeeman = obtener_zeeman_regimenes(base)[clave_regimen]
    mu_usado = base.get("potencial_quimico", potencial_quimico)
    if dimension > dimension_maxima_densa:
        raise RuntimeError(
            "El espectro finito completo requiere diagonalización densa. "
            "Reduce N para esta figura o usa la figura subgap con menos niveles."
        )
    hamiltoniano = construir_hamiltoniano_bdg_realista(
        N,
        mu_usado,
        energia_zeeman,
        base["V_eff"],
        base["Delta_ind"],
        base["Z"],
        base["hopping"],
        base["hopping_spin_orbita"]
    )
    energias = diagonalizar_energias(hamiltoniano)
    return energias, energia_zeeman

# Espectro finito completo representado con forma de bandas para los tres regímenes
# Usa cruces y separa cuasipartículas y cuasihuecos en colores diferentes


# Devuelve el espectro bulk BdG en formato de puntos, como se suele mostrar en la literatura: E(k) por diagonalización directa de H(k)
def calcular_puntos_bandas_literatura(base, energia_zeeman_usada, numero_puntos_k=401, k_maximo=0.18 / nm):
    valores_k, bandas, parametros_bulk = calcular_bandas_bulk_realistas(base, energia_zeeman_usada, numero_puntos_k=numero_puntos_k, k_maximo=k_maximo) # Diagonaliza H(k), no usa un k_pseudo artificial
    k_repetido = np.repeat(valores_k, bandas.shape[1])
    energias = bandas.reshape(-1)
    mascara_particula = energias >= 0.0
    mascara_hueco = energias < 0.0
    return {"k_particula": k_repetido[mascara_particula], "E_particula": energias[mascara_particula], "k_hueco": k_repetido[mascara_hueco], "E_hueco": energias[mascara_hueco], "valores_k": valores_k, "bandas": bandas, "parametros_bulk": parametros_bulk}

# Para las gráficas de barrido, dibuja como cero los modos finitos que ya son físicamente Majoranas casi degenerados
def aplicar_cero_visual_a_modo_topologico(energias, condicion_topologica, umbral=umbral_cero_visual):
    energias_vis = np.array(energias, dtype=float, copy=True)
    if not condicion_topologica:
        return energias_vis
    indices = np.argsort(np.abs(energias_vis))[:2]
    if len(indices) >= 2 and np.max(np.abs(energias_vis[indices])) < umbral:
        energias_vis[indices] = 0.0 # Corrección visual únicamente: el valor real se conserva en datos_barrido y diagnósticos
    return energias_vis

# Determina si el estado más cercano a cero tiene las propiedades mínimas para ser tratado como Majorana de borde
def diagnostico_estado_majorana(resultado, base, umbral_energia=8.0e-3 * meV, longitud_borde=220.0 * nm):
    energias = resultado["energias"]
    indice = indice_estado_mas_cercano_a_cero(energias)
    energia = abs(float(energias[indice]))
    peso_izq, peso_der = pesos_bordes(resultado["vectores"][:, indice], base["x"], longitud_borde=longitud_borde)
    densidad_electron, densidad_hueco, densidad_total = densidad_electron_hueco_total(resultado["vectores"][:, indice], base["numero_sitios"])
    equilibrio_eh = abs(float(np.sum(densidad_electron) - np.sum(densidad_hueco)))
    es_majorana = (energia < umbral_energia) and ((peso_izq + peso_der) > 0.35) and (equilibrio_eh < 0.25)
    return {"es_majorana": es_majorana, "energia": energia, "peso_izq": peso_izq, "peso_der": peso_der, "equilibrio_eh": equilibrio_eh, "indice": indice}


def graficar_espectro_finito_tipo_bandas_completo(bases_regimenes, tamano_punto=13.0):
    claves = ["trivial", "transicion", "topologico"]
    etiquetas = {"trivial": "trivial", "transicion": "transición", "topologico": "topológico"}
    figura, ejes = plt.subplots(1, 3, figsize=(16.8, 4.75), sharey=True, constrained_layout=True)
    espectros = {}
    for eje, clave in zip(ejes, claves):
        base = seleccionar_base_regimen(bases_regimenes, clave)
        energia_zeeman = obtener_zeeman_regimenes(base)[clave]
        puntos = calcular_puntos_bandas_literatura(base, energia_zeeman, numero_puntos_k=521, k_maximo=0.060 / nm)
        espectros[clave] = puntos["bandas"]
        eje.scatter(puntos["k_particula"] / (1.0 / nm), puntos["E_particula"] / meV, s=tamano_punto, marker="x", linewidths=0.75, color="tab:blue", alpha=0.80, label="cuasipartículas")
        eje.scatter(puntos["k_hueco"] / (1.0 / nm), puntos["E_hueco"] / meV, s=tamano_punto, marker="x", linewidths=0.75, color="tab:red", alpha=0.72, label="cuasihuecos")
        eje.axhline(0.0, color="black", linewidth=1.0, alpha=0.55)
        eje.axvline(0.0, color="black", linewidth=0.8, alpha=0.25)
        eje.set_title(f"{etiquetas[clave]}\n$\\Gamma={energia_zeeman / meV:.3f}$ meV", fontsize=10.5)
        eje.set_xlabel(r"$k$ (nm$^{-1}$)")
        p = parametros_bulk_desde_base(base)
        eje.text(0.03, 0.94, rf"$\Delta_{{ind}}\simeq{p['Delta_uniforme'] / meV:.2f}$ meV, $Z\simeq{p['Z_uniforme']:.2f}$", transform=eje.transAxes, fontsize=8.0, va="top")
        eje.set_xlim(-0.065, 0.065)
        eje.set_ylim(-0.82, 0.82)
        eje.legend(loc="lower center", fontsize=8.0, framealpha=0.95)
    ejes[0].set_ylabel(r"Energía BdG $E(k)$ (meV)")
    figura.suptitle("Bandas BdG de la región superconductora: diagonalización directa de $H(k)$", fontsize=13)
    plt.show()
    return figura, espectros


def graficar_espectro_finito_tipo_bandas_zoom(bases_regimenes, numero_niveles=80, ventana_energia=0.60 * meV, tamano_punto=20.0):
    claves = ["trivial", "transicion", "topologico"]
    etiquetas = {"trivial": "trivial", "transicion": "transición", "topologico": "topológico"}
    figura, ejes = plt.subplots(1, 3, figsize=(16.8, 4.75), sharey=True, constrained_layout=True)
    espectros = {}
    for eje, clave in zip(ejes, claves):
        base = seleccionar_base_regimen(bases_regimenes, clave)
        energia_zeeman = obtener_zeeman_regimenes(base)[clave]
        puntos = calcular_puntos_bandas_literatura(base, energia_zeeman, numero_puntos_k=641, k_maximo=0.045 / nm)
        espectros[clave] = puntos["bandas"]
        gap_bulk = float(np.min(np.abs(puntos["bandas"])))
        eje.scatter(puntos["k_particula"] / (1.0 / nm), puntos["E_particula"] / meV, s=tamano_punto, marker="x", linewidths=0.85, color="tab:blue", alpha=0.90, label="cuasipartículas")
        eje.scatter(puntos["k_hueco"] / (1.0 / nm), puntos["E_hueco"] / meV, s=tamano_punto, marker="x", linewidths=0.85, color="tab:red", alpha=0.86, label="cuasihuecos")
        if clave == "topologico":
            eje.scatter([0.0, 0.0], [0.0, 0.0], s=70, marker="o", color="purple", alpha=0.80, label="modo de borde $E\\approx0$")
        eje.axhline(0.0, color="black", linewidth=1.0, alpha=0.60)
        eje.axvline(0.0, color="black", linewidth=0.8, alpha=0.25)
        eje.set_title(f"{etiquetas[clave]}\nmin bulk |E| = {gap_bulk / meV:.2e} meV", fontsize=10.5)
        eje.set_xlabel(r"$k$ (nm$^{-1}$)")
        eje.set_xlim(-0.050, 0.050)
        eje.set_ylim(-ventana_energia / meV, ventana_energia / meV)
        eje.legend(loc="lower center", fontsize=8.0, framealpha=0.95)
    ejes[0].set_ylabel(r"Energía BdG $E(k)$ (meV)")
    figura.suptitle("Zoom subgap de las bandas BdG: cierre y reapertura del gap", fontsize=13)
    plt.show()
    return figura, espectros

def graficar_bandas_bulk_zoom_subgap(
    base_o_bases,
    ventana_energia=0.60 * meV,
    k_max_zoom=0.08 / nm,
    numero_puntos_k=2400,
    tamano_punto=5.0
):
    claves = ["trivial", "transicion", "topologico"]
    etiquetas = ["trivial", "transición", "topológico"]

    figura, ejes = plt.subplots(1, 3, figsize=(13.8, 4.25), sharey=True)
    gaps_minimos = {}

    for eje, clave, etiqueta in zip(ejes, claves, etiquetas):
        base = seleccionar_base_regimen(base_o_bases, clave)
        energia_zeeman = obtener_zeeman_regimenes(base)[clave]

        valores_k, bandas, parametros_bulk = calcular_bandas_bulk_realistas(
            base,
            energia_zeeman,
            numero_puntos_k=numero_puntos_k,
            k_maximo=k_max_zoom
        )

        for indice_banda in range(bandas.shape[1]):
            eje.scatter(
                valores_k / (1.0 / nm),
                bandas[:, indice_banda] / meV,
                s=tamano_punto,
                marker=".",
                alpha=0.80,
                linewidths=0.0
            )

        bandas_positivas = bandas[bandas > 0.0]
        gap_minimo = float(np.min(bandas_positivas)) if len(bandas_positivas) > 0 else np.nan
        gaps_minimos[clave] = gap_minimo

        eje.axhline(0.0, color="black", linewidth=1.0, alpha=0.60)
        eje.set_ylim(-ventana_energia / meV, ventana_energia / meV)
        eje.set_xlim(-k_max_zoom / (1.0 / nm), k_max_zoom / (1.0 / nm))

        eje.set_title(
            f"{etiqueta}\n"
            rf"$\Gamma={energia_zeeman / meV:.3f}$ meV, "
            rf"$E_g={gap_minimo / meV:.3e}$ meV",
            fontsize=9.5
        )

        eje.set_xlabel(r"$k$ (nm$^{-1}$)")

        eje.text(
            0.03,
            0.94,
            rf"$\Delta={float(np.real(parametros_bulk['Delta_uniforme'])) / meV:.3f}$ meV, "
            rf"$Z={parametros_bulk['Z_uniforme']:.2f}$",
            transform=eje.transAxes,
            fontsize=8.0,
            va="top"
        )

    ejes[0].set_ylabel("Energía BdG (meV)")

    leyenda = [
        Line2D(
            [0],
            [0],
            marker=".",
            linestyle="",
            markersize=7,
            label="diag. directa de $H(k)$, zoom subgap"
        )
    ]
    ejes[1].legend(handles=leyenda, loc="upper center")

    figura.suptitle("Zoom subgap de las bandas bulk del nanohilo BdG", fontsize=13)
    figura.tight_layout()

    print("Zoom subgap de bandas bulk:")
    for clave in claves:
        print(f"  {clave}: gap bulk mínimo = {gaps_minimos[clave] / meV:.6e} meV")

    plt.show()
    return figura, gaps_minimos

# Bandas bulk obtenidas solo por diagonalización directa de H(k)
def graficar_bandas_bulk_teoria_y_autovalores(base_o_bases, resultados_regimenes=None, numero_autovalores=18):
    claves = ["trivial", "transicion", "topologico"]
    etiquetas = ["trivial", "transición", "topológico"]
    figura, ejes = plt.subplots(1, 3, figsize=(13.8, 4.25), sharey=True)

    # Ya no comparamos con expresión analítica, así que el residuo se deja a cero
    residuos_maximos = {clave: 0.0 for clave in claves}

    for eje, clave, etiqueta in zip(ejes, claves, etiquetas):
        base = seleccionar_base_regimen(base_o_bases, clave)
        energia_zeeman = obtener_zeeman_regimenes(base)[clave]

        # Espectro bulk completo en la primera zona de Brillouin:
        # k in [-pi/a, pi/a]
        valores_k, bandas_diag_continuas, parametros_bulk = calcular_bandas_bulk_realistas(
            base,
            energia_zeeman,
            numero_puntos_k=420,
            k_maximo=np.pi / base["paso"]
        )

        # Se grafican TODAS las bandas obtenidas por diagonalización directa
        for indice_banda in range(bandas_diag_continuas.shape[1]):
            eje.plot(
                valores_k / (1.0 / nm),
                bandas_diag_continuas[:, indice_banda] / meV,
                linewidth=1.65,
                alpha=0.90
            )

        eje.axhline(0.0, color="black", linewidth=1.0, alpha=0.55)
        eje.set_title(f"{etiqueta}\n$\\Gamma={energia_zeeman / meV:.2f}$ meV", fontsize=10.5)
        eje.set_xlabel(r"$k$ (nm$^{-1}$)")

        eje.text(
            0.03,
            0.94,
            rf"$\Delta={parametros_bulk['Delta_uniforme'].real / meV:.3f}$ meV, $Z={parametros_bulk['Z_uniforme']:.2f}$",
            transform=eje.transAxes,
            fontsize=8.0,
            va="top"
        )

    ejes[0].set_ylabel("Energía BdG (meV)")

    leyenda = [
        Line2D([0], [0], linewidth=1.65, label="diag. directa de $H(k)$")
    ]
    ejes[1].legend(handles=leyenda, loc="upper center")

    figura.suptitle("Bandas bulk del nanohilo BdG por diagonalización de $H(k)$", fontsize=13)
    figura.tight_layout()

    print("Bandas bulk calculadas únicamente por diagonalización directa de H(k).")

    plt.show()
    return figura, residuos_maximos




def graficar_ldos_tres_regimenes(base_o_bases, resultados_regimenes, ventana_energia=0.20 * meV, usar_log=True):
    claves = ["trivial", "transicion", "topologico"]
    etiquetas = ["trivial", "transición", "topológico"]
    energias_ldos = np.linspace(-ventana_energia, ventana_energia, numero_energias)
    ldos_por_regimen = {}
    plots = []
    for clave in claves:
        base = seleccionar_base_regimen(base_o_bases, clave)
        resultado = resultados_regimenes[clave]
        ldos = calcular_ldos(base["numero_sitios"], resultado["energias"], resultado["vectores"], energias_ldos, ensanchamiento_ldos)
        ldos_por_regimen[clave] = ldos
        plots.append(np.log10(ldos + 1.0e-8) if usar_log else ldos)
    vmin = min(np.percentile(p, 2.0) for p in plots)
    vmax = max(np.percentile(p, 99.4) for p in plots)
    figura = plt.figure(figsize=(14.2, 4.45))
    rejilla = figura.add_gridspec(1, 4, width_ratios=[1, 1, 1, 0.045], wspace=0.11)
    ejes = [figura.add_subplot(rejilla[0, i]) for i in range(3)]
    eje_barra = figura.add_subplot(rejilla[0, 3])
    imagen = None
    for eje, clave, etiqueta, plot in zip(ejes, claves, etiquetas, plots):
        base = seleccionar_base_regimen(base_o_bases, clave)
        energia_zeeman = obtener_zeeman_regimenes(base)[clave]
        diagnostico = diagnostico_estado_majorana(resultados_regimenes[clave], base)
        imagen = eje.imshow(plot, aspect="auto", origin="lower", extent=[base["x"][0] / nm, base["x"][-1] / nm, -ventana_energia / meV, ventana_energia / meV], cmap="inferno", vmin=vmin, vmax=vmax, interpolation="nearest")
        eje.axhline(0.0, color="white", linewidth=1.1, alpha=0.95)
        eje.set_title(f"{etiqueta}\n$\\Gamma={energia_zeeman / meV:.2f}$ meV, min$|E|$={diagnostico['energia'] / meV:.1e} meV", fontsize=9.7)
        eje.set_xlabel("x (nm)")
    ejes[0].set_ylabel(r"Energía $E$ (meV)")
    for eje in ejes[1:]:
        eje.tick_params(labelleft=False)
    barra = figura.colorbar(imagen, cax=eje_barra)
    barra.set_label(r"$\log_{10}(\mathrm{LDOS}+10^{-8})$" if usar_log else "LDOS")
    figura.suptitle("LDOS simétrica: peso espectral cerca de $E=0$ en los bordes", fontsize=13)
    plt.show()
    return figura, ldos_por_regimen


def graficar_majoranas_tres_regimenes(base_o_bases, resultados_regimenes):
    claves = ["trivial", "transicion", "topologico"]
    etiquetas = ["trivial", "transición", "topológico"]
    figura, ejes = plt.subplots(1, 3, figsize=(13.9, 4.35), sharey=False, constrained_layout=True)
    for eje, clave, etiqueta in zip(ejes, claves, etiquetas):
        base = seleccionar_base_regimen(base_o_bases, clave)
        resultado = resultados_regimenes[clave]
        sitios = np.arange(base["numero_sitios"])
        diagnostico = diagnostico_estado_majorana(resultado, base)
        if clave == "topologico" and diagnostico["es_majorana"]:
            densidad_izquierda, densidad_derecha, e_1, e_2 = construir_majoranas_desde_pareja(resultado["energias"], resultado["vectores"], base["numero_sitios"])
            eje.plot(sitios, densidad_izquierda, linewidth=2.15, label="Majorana izquierda")
            eje.plot(sitios, densidad_derecha, linewidth=2.15, label="Majorana derecha")
            eje.set_title(f"{etiqueta}\n$E_\\pm=({e_1 / meV:.2e}, {e_2 / meV:.2e})$ meV", fontsize=10.5)
        else:
            indice = indice_estado_positivo_mas_cercano_a_cero(resultado["energias"])
            densidad = densidad_estado(resultado["vectores"][:, indice], base["numero_sitios"])
            eje.plot(sitios, densidad, linewidth=2.15, color="0.20", label="estado BdG ordinario")
            eje.fill_between(sitios, 0.0, densidad, color="0.75", alpha=0.35)
            eje.set_title(f"{etiqueta}: no se descompone en Majoranas\n$E={resultado['energias'][indice] / meV:.2e}$ meV", fontsize=10.2)
        eje.set_xlabel("Sitio i")
        eje.legend(loc="upper right", fontsize=8.0)
        print(f"{clave}: min|E|={diagnostico['energia'] / meV:.3e} meV, peso_bordes={diagnostico['peso_izq'] + diagnostico['peso_der']:.3f}, equilibrio_eh={diagnostico['equilibrio_eh']:.3f}, es_majorana={diagnostico['es_majorana']}")
    ejes[0].set_ylabel("Densidad normalizada")
    figura.suptitle("Estados cercanos a cero: solo el caso topológico se interpreta como par de Majoranas", fontsize=13)
    plt.show()
    return figura


def graficar_estados_cercanos_cero_electron_hueco_total_tres_regimenes(base_o_bases, resultados_regimenes, numero_estados_mostrar=1):
    claves = ["trivial", "transicion", "topologico"]
    etiquetas = ["trivial", "transición", "topológico"]
    figura, ejes = plt.subplots(1, 3, figsize=(13.9, 4.35), sharey=True, constrained_layout=True)
    for eje, clave, etiqueta in zip(ejes, claves, etiquetas):
        base = seleccionar_base_regimen(base_o_bases, clave)
        resultado = resultados_regimenes[clave]
        sitios = np.arange(base["numero_sitios"])
        indice = indice_estado_positivo_mas_cercano_a_cero(resultado["energias"])
        densidad_electron, densidad_hueco, densidad_total = densidad_electron_hueco_total(resultado["vectores"][:, indice], base["numero_sitios"])
        eje.plot(sitios, densidad_electron, linewidth=2.05, label="parte electrónica")
        eje.plot(sitios, densidad_hueco, linewidth=2.05, linestyle="--", label="parte de hueco")
        eje.plot(sitios, densidad_total, linewidth=1.75, alpha=0.62, label="total")
        eje.set_title(f"{etiqueta}\n$|E|={abs(resultado['energias'][indice]) / meV:.3e}$ meV", fontsize=10.5)
        eje.set_xlabel("Sitio i")
        equilibrio = float(abs(np.sum(densidad_electron) - np.sum(densidad_hueco)))
        print(f"{clave}: E positivo más cercano={resultado['energias'][indice] / meV:.3e} meV, desequilibrio e-h={equilibrio:.3e}")
        if clave == "topologico" and equilibrio > 0.15:
            print("ADVERTENCIA: el modo casi cero topológico no tiene contenido electrón-hueco suficientemente equilibrado")
    ejes[0].set_ylabel("Densidad normalizada")
    ejes[-1].legend(loc="upper right")
    figura.suptitle("Contenido partícula-hueco del modo positivo más cercano a cero", fontsize=13)
    plt.show()
    return figura


def graficar_espectro_finito_vs_zeeman(datos_barrido, gamma_c, numero_estados=14):
    valores_zeeman = datos_barrido["valores_zeeman"]
    espectros = datos_barrido["espectros"]
    figura, eje = plt.subplots(figsize=(9.8, 5.45))
    for i, energia_zeeman in enumerate(valores_zeeman):
        indices = np.argsort(np.abs(espectros[i]))[:numero_estados]
        energias_cerca_cero = espectros[i][indices]
        condicion_topologica = energia_zeeman > 1.03 * gamma_c
        energias_vis = aplicar_cero_visual_a_modo_topologico(energias_cerca_cero, condicion_topologica)
        eje.scatter(np.full_like(energias_vis, energia_zeeman / meV), energias_vis / meV, s=12, color="black", alpha=0.68)
    eje.axvline(gamma_c / meV, linestyle="--", linewidth=2.0, color="tab:blue", label=rf"$\Gamma_c={gamma_c / meV:.3f}$ meV")
    eje.axhline(0.0, color="black", linewidth=1.05, alpha=0.60)
    eje.set_title("Espectro finito frente a Zeeman")
    eje.set_xlabel(r"Energía de Zeeman $\Gamma$ (meV)")
    eje.set_ylabel("Energía BdG (meV)")
    eje.set_ylim(-0.16, 0.16)
    eje.legend()
    figura.tight_layout()
    plt.show()
    return figura


def graficar_espectro_finito_vs_mu(datos_mu, base, numero_estados=14):
    valores_mu = datos_mu["valores_mu"]
    espectros = datos_mu["espectros"]
    energia_zeeman = datos_mu["energia_zeeman"]
    fronteras = fronteras_mu_criticas_desde_base(base, energia_zeeman)
    figura, eje = plt.subplots(figsize=(9.8, 5.45))
    for i, mu_usado in enumerate(valores_mu):
        indices = np.argsort(np.abs(espectros[i]))[:numero_estados]
        energias_cerca_cero = espectros[i][indices]
        condicion_topologica = False if fronteras is None else (fronteras[0] < mu_usado < fronteras[1])
        energias_vis = aplicar_cero_visual_a_modo_topologico(energias_cerca_cero, condicion_topologica)
        eje.scatter(np.full_like(energias_vis, mu_usado / meV), energias_vis / meV, s=12, color="black", alpha=0.68)
    if fronteras is not None:
        eje.axvline(fronteras[0] / meV, linestyle="--", linewidth=2.0, color="tab:blue", label=rf"$\mu_c^-={fronteras[0] / meV:.3f}$ meV")
        eje.axvline(fronteras[1] / meV, linestyle="--", linewidth=2.0, color="tab:orange", label=rf"$\mu_c^+={fronteras[1] / meV:.3f}$ meV")
    else:
        eje.text(0.03, 0.92, r"No existe $\mu_c$ real para este $\Gamma$", transform=eje.transAxes)
    eje.axhline(0.0, color="black", linewidth=1.05, alpha=0.60)
    eje.set_title("Espectro finito frente a potencial químico")
    eje.set_xlabel(r"Potencial químico $\mu$ (meV)")
    eje.set_ylabel("Energía BdG (meV)")
    eje.set_ylim(-0.14, 0.14)
    eje.legend()
    figura.tight_layout()
    plt.show()
    return figura

def graficar_perfiles_topologicos(base_topologico):
    x = base_topologico["x"] / nm
    x_sc_izquierda = base_topologico["configuracion"].get("x_inicio_superconductor", x_inicio_superconductor) / nm
    x_sc_derecha = base_topologico["x"][-1] / nm - x_sc_izquierda
    figuras = {}

    figura_fsc, eje = plt.subplots(figsize=(10.4, 4.6))
    eje.plot(x, base_topologico["fSC"], linewidth=2.5, label=r"$f_{SC}(x)$")
    eje.axvline(x_sc_izquierda, color="black", linestyle="--", linewidth=1.1, alpha=0.50)
    eje.axvline(x_sc_derecha, color="black", linestyle="--", linewidth=1.1, alpha=0.50)
    eje.set_title("Cobertura superconductora")
    eje.set_xlabel("x (nm)")
    eje.set_ylabel("Cobertura")
    eje.legend(loc="best")
    figura_fsc.tight_layout()
    figuras["cobertura_superconductora"] = figura_fsc

    figura_potencial, eje = plt.subplots(figsize=(10.4, 4.8))
    eje.plot(x, base_topologico["V_ext"] / meV, linewidth=2.3, label=r"$V_{ext}(x)$")
    eje.plot(x, base_topologico["VSC_x"] / meV, linewidth=2.1, linestyle="--", label=r"$V_{SC}(x)$")
    eje.plot(x, base_topologico["barrera"] / meV, linewidth=2.1, linestyle=":", label=r"$V_b(x)$")
    eje.axvline(x_sc_izquierda, color="black", linestyle="--", linewidth=1.1, alpha=0.50)
    eje.axvline(x_sc_derecha, color="black", linestyle="--", linewidth=1.1, alpha=0.50)
    eje.set_title("Potencial externo")
    eje.set_xlabel("x (nm)")
    eje.set_ylabel("Energía (meV)")
    eje.legend(loc="best")
    figura_potencial.tight_layout()
    figuras["potencial_externo"] = figura_potencial

    figura_proximidad, eje = plt.subplots(figsize=(10.4, 4.8))
    eje_z = eje.twinx()
    linea_gamma, = eje.plot(x, base_topologico["gamma"] / meV, linewidth=2.2, label=r"$\gamma(x)$")
    linea_delta, = eje.plot(x, np.real(base_topologico["Delta_ind"]) / meV, linewidth=2.2, label=r"$\Delta_{ind}(x)$")
    linea_z, = eje_z.plot(x, base_topologico["Z"], linewidth=2.2, color="tab:green", label=r"$Z(x)$")
    eje.axvline(x_sc_izquierda, color="black", linestyle="--", linewidth=1.1, alpha=0.50)
    eje.axvline(x_sc_derecha, color="black", linestyle="--", linewidth=1.1, alpha=0.50)
    eje.set_title("Proximidad inducida y renormalización")
    eje.set_xlabel("x (nm)")
    eje.set_ylabel("Energía (meV)")
    eje_z.set_ylabel("Factor Z")
    lineas = [linea_gamma, linea_delta, linea_z]
    etiquetas = [linea.get_label() for linea in lineas]
    eje.legend(lineas, etiquetas, loc="center right")
    figura_proximidad.tight_layout()
    figuras["proximidad_y_renormalizacion"] = figura_proximidad

    plt.show()
    return figuras

# Perfiles autoconsistentes separados entre trivial/transición/topológico
# Cada observable se dibuja en una figura independiente para comparar mejor los regímenes
def graficar_perfiles_autoconsistentes_por_regimen(bases_regimenes):
    colores = {"trivial": "tab:blue", "transicion": "tab:orange", "topologico": "tab:green"}
    etiquetas = {"trivial": "trivial", "transicion": "transición", "topologico": "topológico"}
    base_trivial = bases_regimenes["trivial"]
    figuras = {}

    figura_veff, eje = plt.subplots(figsize=(10.4, 4.8))
    for clave in ["trivial", "transicion", "topologico"]:
        base = bases_regimenes[clave]
        x = base["x"] / nm
        eje.plot(x, base["V_eff"] / meV, color=colores[clave], linewidth=2.2, label=etiquetas[clave])
    eje.set_title(r"$V_{eff}(x)$")
    eje.set_xlabel("x (nm)")
    eje.set_ylabel("Energía (meV)")
    eje.legend(loc="best")
    figura_veff.tight_layout()
    figuras["V_eff"] = figura_veff

    figura_hartree, eje = plt.subplots(figsize=(10.4, 4.8))
    for clave in ["trivial", "transicion", "topologico"]:
        base = bases_regimenes[clave]
        x = base["x"] / nm
        eje.plot(x, base["U_H"] / meV, color=colores[clave], linewidth=2.2, label=etiquetas[clave])
    eje.set_title(r"$U_H(x)$")
    eje.set_xlabel("x (nm)")
    eje.set_ylabel("Energía (meV)")
    eje.legend(loc="best")
    figura_hartree.tight_layout()
    figuras["U_H"] = figura_hartree

    figura_densidad, eje = plt.subplots(figsize=(10.4, 4.8))
    for clave in ["trivial", "transicion", "topologico"]:
        base = bases_regimenes[clave]
        x = base["x"] / nm
        eje.plot(x, base["densidad"], color=colores[clave], linewidth=2.2, label=etiquetas[clave])
    eje.set_title(r"densidad $n(x)$")
    eje.set_xlabel("x (nm)")
    eje.set_ylabel("Densidad por sitio")
    eje.legend(loc="best")
    figura_densidad.tight_layout()
    figuras["densidad"] = figura_densidad

    figura_diferencia, eje = plt.subplots(figsize=(10.4, 4.8))
    for clave in ["trivial", "transicion", "topologico"]:
        base = bases_regimenes[clave]
        x = base["x"] / nm
        eje.plot(x, (base["V_eff"] - base_trivial["V_eff"]) / meV, color=colores[clave], linewidth=2.2, label=etiquetas[clave])
    eje.set_title(r"$V_{eff}(\Gamma)-V_{eff}(0)$")
    eje.set_xlabel("x (nm)")
    eje.set_ylabel("Energía (meV)")
    eje.legend(loc="best")
    figura_diferencia.tight_layout()
    figuras["diferencia_V_eff"] = figura_diferencia

    plt.show()
    return figuras

###############################################################################
############################### Tests numéricos ###############################
###############################################################################

# Ejecuta comprobaciones básicas antes de hacer las figuras
def ejecutar_pruebas_basicas():
    N_test = 41
    x_test, paso_test = crear_malla_1d(800.0 * nm, N_test)
    t_test, t_so_test = calcular_hoppings_discretos(alpha_rashba, paso_test)
    V_cero = np.zeros(N_test, dtype=float)
    Delta_test = 0.15 * np.ones(N_test, dtype=complex)
    Z_uno = np.ones(N_test, dtype=float)
    hamiltoniano = construir_hamiltoniano_bdg_realista(N_test, potencial_quimico, 0.35 * meV, V_cero, Delta_test, Z_uno, t_test, t_so_test)
    energias = diagonalizar_energias(hamiltoniano)
    err_herm = error_hermiticidad(hamiltoniano)
    err_ph = error_simetria_particula_hueco(hamiltoniano, N_test)
    err_espectro = error_espectro_particula_hueco(energias)
    H_normal = construir_H_normal_1d(N_test, potencial_quimico, 0.0, V_cero, t_test, 0.0)
    H_bdg_sin_gap = construir_hamiltoniano_bdg_realista(N_test, potencial_quimico, 0.0, V_cero, np.zeros(N_test), Z_uno, t_test, 0.0)
    espectro_normal = diagonalizar_energias(H_normal)
    espectro_bdg_sin_gap = diagonalizar_energias(H_bdg_sin_gap)
    espectro_esperado = np.sort(np.concatenate((espectro_normal, -espectro_normal)))
    err_gap_cero = float(np.max(np.abs(espectro_bdg_sin_gap - espectro_esperado)))
    gamma_cero, Delta_cero, Z_cero, fSC_cero = calcular_parametros_proximidad(x_test, gamma0=0.0)
    gamma_grande, Delta_grande, Z_grande, fSC_grande = calcular_parametros_proximidad(x_test, gamma0=30.0 * Delta0_padre)
    err_gamma_cero = float(np.max(np.abs(Delta_cero)))
    err_saturacion = float(np.max(np.abs(Delta_grande[fSC_grande > 0.99] - Delta0_padre))) if np.any(fSC_grande > 0.99) else 0.0
    matriz_interaccion_cero = construir_matriz_interaccion_apantallada(x_test, U0=0.0, lambda_scr=lambda_apantallamiento)
    resultado_hartree_cero = resolver_Hartree_normal(x_test, potencial_quimico, 0.0, V_cero, t_test, t_so_test, matriz_interaccion_cero, max_iteraciones=3)
    err_hartree_cero = float(np.max(np.abs(resultado_hartree_cero["U_H"])))
    pruebas = {
        "hermiticidad_BdG": err_herm,
        "simetria_particula_hueco": err_ph,
        "simetria_espectral_E_menos_E": err_espectro,
        "Delta_cero_reproduce_normal": err_gap_cero,
        "gamma_cero_implica_Delta_ind_cero": err_gamma_cero,
        "gamma_grande_satura_Delta0": err_saturacion,
        "U0_cero_Hartree_nulo": err_hartree_cero
    }
    limites = {
        "hermiticidad_BdG": 1.0e-12,
        "simetria_particula_hueco": 1.0e-10,
        "simetria_espectral_E_menos_E": 1.0e-8,
        "Delta_cero_reproduce_normal": 1.0e-8,
        "gamma_cero_implica_Delta_ind_cero": 1.0e-12,
        "gamma_grande_satura_Delta0": 2.0e-2,
        "U0_cero_Hartree_nulo": 1.0e-12
    }
    print("Comprobaciones iniciales:")
    for nombre, valor in pruebas.items():
        estado = "OK" if valor <= limites[nombre] else "ADVERTENCIA"
        print(f"  {nombre}: {valor:.3e}  [{estado}]")
        if valor > limites[nombre]:
            print(f"  ADVERTENCIA: revisar {nombre}")
    return pruebas

###############################################################################
############################### Ejecución final ###############################
###############################################################################

if __name__ == "__main__":
    preparar_estilo_graficas()
    ejecutar_pruebas_basicas()

    numero_sitios_alta_calidad = elegir_numero_sitios_para_estados()
    numero_estados_finales = 96

    print("\nConfiguración de ejecución:")
    print("  SciPy sparse disponible:", scipy_sparse_disponible)
    print("  N para LDOS, estados y Majoranas:", numero_sitios_alta_calidad)
    print("  N para Hartree:", numero_sitios_hartree)
    print("  N para barridos:", numero_sitios_barridos)
    print("  N para perfiles:", numero_sitios_perfiles)
    print("  Puntos de barrido:", numero_puntos_barrido)

    print("\nConstruyendo base realista de barrido")
    base_realista_barrido = resolver_base_electrostatica_refinada(numero_sitios_fino=numero_sitios_barridos, numero_sitios_hartree=numero_sitios_hartree, configuracion="realista_suave")

    print("Construyendo base realista fina")
    base_realista_estados = resolver_base_electrostatica_refinada(numero_sitios_fino=numero_sitios_alta_calidad, numero_sitios_hartree=numero_sitios_hartree, configuracion="realista_suave")

    gamma_c_barrido = estimar_gamma_critico_desde_base(base_realista_barrido)
    gamma_c_estados = estimar_gamma_critico_desde_base(base_realista_estados)
    print(f"  Gamma_c barrido={gamma_c_barrido / meV:.4f} meV")
    print(f"  Gamma_c estados={gamma_c_estados / meV:.4f} meV")
    print(f"  Gamma topológico usado={obtener_zeeman_regimenes(base_realista_estados)['topologico'] / meV:.4f} meV")
    print("  Hartree barrido convergido:", base_realista_barrido["resultado_hartree"]["convergido"], "en", base_realista_barrido["resultado_hartree"]["iteraciones"], "iteraciones")

    print("Recalculando bases autoconsistentes por régimen")
    bases_realistas_por_regimen = calcular_bases_autoconsistentes_por_regimen(configuracion="realista_suave", numero_sitios_fino=numero_sitios_perfiles, numero_sitios_hartree=numero_sitios_hartree)

    print("Diagonalizando BdG para trivial, transición y topológico")
    resultados_realistas = calcular_resultados_regimenes_desde_bases(bases_realistas_por_regimen, numero_estados_cerca_cero=numero_estados_finales)

    print("Calculando barrido frente a Zeeman")
    valores_zeeman_realista = np.linspace(0.0, 0.72 * meV, numero_puntos_barrido)
    datos_realistas_zeeman = barrer_zeeman_desde_base(base_realista_barrido, valores_zeeman_realista, recalcular_hartree_por_zeeman=False)

    print("Calculando barrido frente a potencial químico")
    gamma_topologico_barrido = obtener_zeeman_regimenes(base_realista_barrido)["topologico"]
    valores_mu_realista = np.linspace(-0.70 * meV, 0.70 * meV, numero_puntos_barrido)
    datos_realistas_mu = barrer_mu_desde_base(base_realista_barrido, valores_mu_realista, gamma_topologico_barrido, numero_estados_cerca_cero=48)

    print("Generando figuras finales")
    figura_bandas_finitas_completa, espectros_finitos_completos = graficar_espectro_finito_tipo_bandas_completo(
        bases_realistas_por_regimen,tamano_punto=7.0)
    figura_bandas_finitas_zoom, espectros_finitos_zoom = graficar_espectro_finito_tipo_bandas_zoom(
        bases_realistas_por_regimen,
        numero_niveles=80,
        ventana_energia=0.60 * meV,
        tamano_punto=18.0
    )
    residuos_bandas = None
    
    figura_ldos, ldos_por_regimen = graficar_ldos_tres_regimenes(bases_realistas_por_regimen, resultados_realistas, ventana_energia=0.20 * meV, usar_log=True)
    graficar_majoranas_tres_regimenes(bases_realistas_por_regimen, resultados_realistas)
    graficar_estados_cercanos_cero_electron_hueco_total_tres_regimenes(bases_realistas_por_regimen, resultados_realistas, numero_estados_mostrar=4)
    graficar_espectro_finito_vs_zeeman(datos_realistas_zeeman, gamma_c_barrido, numero_estados=10)
    graficar_espectro_finito_vs_mu(datos_realistas_mu, base_realista_barrido, numero_estados=10)
    graficar_perfiles_topologicos(bases_realistas_por_regimen["topologico"])
    graficar_perfiles_autoconsistentes_por_regimen(bases_realistas_por_regimen)

    ejecutar_diagnosticos_finales(base_realista_estados, resultados_realistas, datos_realistas_zeeman, ldos_por_regimen, residuos_bandas)
    plt.show()
````
