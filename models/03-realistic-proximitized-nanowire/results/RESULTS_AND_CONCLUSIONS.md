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
x_i=linspace(0,L,N),\qquad a=\frac{L}{N-1}.
$$

<p align="right"><strong>(3.1)</strong></p>

The endpoint-inclusive mesh makes the physical wire length exactly $L$.

### Discrete hopping scales

`calcular_hoppings_discretos` evaluates 

$$
t=\frac{\hbar^2}{2m^*a^2},\qquad t_{SO}=\frac{\alpha_R}{2a}.
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

