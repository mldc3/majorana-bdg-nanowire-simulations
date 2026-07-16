# Self-Consistent Bogoliubov–de Gennes Nanowire

## Complete theoretical development from superconducting mean field to the lattice equations used in the simulation

Continuation of the Kitaev-chain documentation for the repository `majorana-bdg-nanowire-simulations`. The material starts from the second-quantized framework already established there and develops the Bogoliubov–de Gennes formalism and the self-consistent nanowire model.

**Technical study document in English**

**Source model:** `src/self_consistent_bdg_nanowire.py`

---

## Contents

1. [Scope and relation to the Kitaev-chain section](#1-scope-and-relation-to-the-kitaev-chain-section)
   1. [Starting point](#11-starting-point)
   2. [What is self-consistent in this model?](#12-what-is-self-consistent-in-this-model)
   3. [Physical system](#13-physical-system)
2. [Why the Bogoliubov–de Gennes formalism is needed](#2-why-the-bogoliubovde-gennes-formalism-is-needed)
   1. [Superconductors do not conserve particle number at mean-field level](#21-superconductors-do-not-conserve-particle-number-at-mean-field-level)
   2. [BCS interpretation](#22-bcs-interpretation)
   3. [BdG is an eigenvalue problem plus a closure condition](#23-bdg-is-an-eigenvalue-problem-plus-a-closure-condition)
3. [Microscopic interacting Hamiltonian](#3-microscopic-interacting-hamiltonian)
   1. [One-body and interaction terms](#31-one-body-and-interaction-terms)
   2. [Local singlet specialization](#32-local-singlet-specialization)
4. [Mean-field decoupling step by step](#4-mean-field-decoupling-step-by-step)
   1. [Fluctuation decomposition](#41-fluctuation-decomposition)
   2. [Definition of the order parameter](#42-definition-of-the-order-parameter)
   3. [Why the result is quadratic](#43-why-the-result-is-quadratic)
5. [Equations of motion](#5-equations-of-motion)
   1. [Fermionic commutator identity](#51-fermionic-commutator-identity)
   2. [Compact matrix form before diagonalization](#52-compact-matrix-form-before-diagonalization)
6. [Bogoliubov transformation](#6-bogoliubov-transformation)
   1. [Canonical quasiparticles](#61-canonical-quasiparticles)
   2. [Diagonalization condition](#62-diagonalization-condition)
7. [Bogoliubov–de Gennes equations](#7-bogoliubovde-gennes-equations)
   1. [Component equations](#71-component-equations)
   2. [Nambu matrix](#72-nambu-matrix)
   3. [Nambu doubling is a redundancy](#73-nambu-doubling-is-a-redundancy)
8. [Particle–hole symmetry in the code basis](#8-particlehole-symmetry-in-the-code-basis)
   1. [Antiunitary operator](#81-antiunitary-operator)
   2. [Matrix check](#82-matrix-check)
9. [Derivation of the self-consistency equation](#9-derivation-of-the-self-consistency-equation)
   1. [Inverse Bogoliubov transformation](#91-inverse-bogoliubov-transformation)
   2. [Spin-singlet anomalous amplitude](#92-spin-singlet-anomalous-amplitude)
   3. [Positive-energy and all-energy forms](#93-positive-energy-and-all-energy-forms)
   4. [Energy cutoff](#94-energy-cutoff)
   5. [Gauge freedom and phase fixing](#95-gauge-freedom-and-phase-fixing)
10. [Continuum model of the Rashba–Zeeman nanowire](#10-continuum-model-of-the-rashbazeeman-nanowire)
    1. [Normal Hamiltonian](#101-normal-hamiltonian)
    2. [BdG Hamiltonian](#102-bdg-hamiltonian)
    3. [Topological criterion in the continuum limit](#103-topological-criterion-in-the-continuum-limit)
11. [Spatial discretization without omitted steps](#11-spatial-discretization-without-omitted-steps)
    1. [Lattice points and Taylor expansions](#111-lattice-points-and-taylor-expansions)
    2. [Discrete kinetic term](#112-discrete-kinetic-term)
    3. [Discrete Rashba term](#113-discrete-rashba-term)
    4. [Local and hopping blocks](#114-local-and-hopping-blocks)
12. [Homogeneous bulk Hamiltonian and bands](#12-homogeneous-bulk-hamiltonian-and-bands)
    1. [Fourier transform](#121-fourier-transform)
    2. [Analytical dispersion](#122-analytical-dispersion)
    3. [Bulk gap](#123-bulk-gap)
13. [Observables derived from BdG eigenvectors](#13-observables-derived-from-bdg-eigenvectors)
    1. [Total BdG density](#131-total-bdg-density)
    2. [Electron and hole sectors](#132-electron-and-hole-sectors)
    3. [Electronic local density of states](#133-electronic-local-density-of-states)
    4. [Majorana combinations](#134-majorana-combinations)
    5. [Finite-length splitting](#135-finite-length-splitting)
14. [Complete equation-to-code map](#14-complete-equation-to-code-map)
15. [Conclusions](#15-conclusions)
16. [Bibliography](#bibliography)

---

# 1 Scope and relation to the Kitaev-chain section

## 1.1 Starting point

The preceding Kitaev-chain documentation introduced Fock space, fermionic creation and annihilation operators, canonical anticommutation relations, fermion parity, Majorana operators, Nambu doubling, and the basic meaning of particle–hole symmetry. Those foundations will not be repeated here. They are assumed throughout this document.

The present objective is different. We begin with an interacting spinful electronic system, derive a superconducting mean-field Hamiltonian, obtain the Bogoliubov–de Gennes (BdG) equations without omitting the intermediate algebra, specialize them to a one-dimensional Rashba nanowire in a Zeeman field, and discretize the continuum operator into the exact block-tridiagonal matrix assembled by the Python script.

> **Logical dependency**
>
> Second quantization supplies the operator language; the BdG construction supplies the quasiparticle eigenproblem; self-consistency closes the problem by requiring the pair potential used in the Hamiltonian to equal the anomalous expectation value calculated from the resulting eigenstates.

## 1.2 What is self-consistent in this model?

The pairing field is not inserted once and held fixed. A trial spatial profile $\Delta^{(0)}(x)$ is used to construct a BdG Hamiltonian. Its eigenvectors determine a new anomalous pair amplitude and hence a new profile $\Delta_{\mathrm{calc}}(x)$. The Hamiltonian is rebuilt with a mixed profile, and the procedure is repeated until successive profiles differ by less than a prescribed tolerance. Thus the superconducting order parameter and the quasiparticle spectrum are solved together.

## 1.3 Physical system

The model describes a one-dimensional spinful electronic channel with effective mass $m^*$, chemical potential $\mu$, a possible electrostatic profile $V(x)$, Rashba spin–orbit coupling $\alpha_R$, a Zeeman energy $E_Z$, and local spin-singlet pairing $\Delta(x)$. In continuum form the single-particle part is

$$
h_0=\frac{p_x^2}{2m^*}-\mu+V(x)+\frac{\alpha_R}{\hbar}p_x\sigma_y+E_Z\sigma_x.
$$

**(1.1)**

The pairing field couples particles and holes. The combination of spin–orbit coupling, Zeeman splitting, and s-wave pairing can produce an effective low-energy spinless superconducting sector and, after a bulk gap closing and reopening, end-localized near-zero modes.

# 2 Why the Bogoliubov–de Gennes formalism is needed

## 2.1 Superconductors do not conserve particle number at mean-field level

A normal quadratic Hamiltonian contains terms of the form $\psi^\dagger\psi$. It moves particles between one-particle states but does not mix sectors with different particle number. A superconducting mean-field Hamiltonian also contains $\psi^\dagger\psi^\dagger$ and $\psi\psi$. These terms create and annihilate pairs. The quasiparticles that diagonalize such a Hamiltonian must therefore mix annihilation and creation operators.

The mean-field Hamiltonian does preserve fermion parity because every term changes particle number by an even integer. This is precisely the algebraic setting in which BdG particle–hole redundancy and Majorana self-conjugacy arise.

## 2.2 BCS interpretation

An effective attraction in the spin-singlet channel correlates time-reversed electronic states. The condensate is characterized by the anomalous average $\langle\psi_\downarrow\psi_\uparrow\rangle$. A nonzero value defines an order parameter $\Delta$ and opens an excitation gap. The elementary excitations are Bogoliubov quasiparticles: coherent electron–hole superpositions whose amplitudes are denoted $u$ and $v$.

## 2.3 BdG is an eigenvalue problem plus a closure condition

For a prescribed $\Delta$, BdG theory is a linear Hermitian eigenproblem. When $\Delta$ is generated by the same electronic degrees of freedom, the eigenproblem is nonlinear as a whole because $\Delta$ depends on its own eigenvectors. The standard solution is fixed-point iteration.

# 3 Microscopic interacting Hamiltonian

## 3.1 One-body and interaction terms

We start from

$$
H=H_0+H_{\mathrm{int}},
$$

**(3.1)**

$$
H_0=\sum_{\alpha\beta}\int dx\,dx'\,\psi_\alpha^\dagger(x)h_{\alpha\beta}(x,x')\psi_\beta(x'),
$$

**(3.2)**

$$
H_{\mathrm{int}}=-\int dx\,dx'\,V_{\mathrm{eff}}(x,x')\psi_\uparrow^\dagger(x)\psi_\downarrow^\dagger(x')\psi_\downarrow(x')\psi_\uparrow(x).
$$

**(3.3)**

The effective interaction is written with $V_{\mathrm{eff}}>0$ and an explicit minus sign, so the selected channel is attractive. The spin indices of the one-body operator may include spin mixing; this is required for Rashba and Zeeman terms.

## 3.2 Local singlet specialization

For a local s-wave model,

$$
V_{\mathrm{eff}}(x,x')=V_{\mathrm{eff}}\delta(x-x'),
$$

**(3.4)**

and the relevant pair operator is $\psi_\downarrow(x)\psi_\uparrow(x)$. The orbital part is even under exchange, so the spin part is antisymmetric, i.e. a singlet. The corresponding order parameter is

$$
\Delta(x)=V_{\mathrm{eff}}\langle\psi_\downarrow(x)\psi_\uparrow(x)\rangle.
$$

**(3.5)**

# 4 Mean-field decoupling step by step

## 4.1 Fluctuation decomposition

Define

$$
A=\psi_\uparrow^\dagger(x)\psi_\downarrow^\dagger(x'),\qquad B=\psi_\downarrow(x')\psi_\uparrow(x).
$$

**(4.1)**

Write $A=\langle A\rangle+\delta A$ and $B=\langle B\rangle+\delta B$. Then

$$
AB=\langle A\rangle\langle B\rangle+\langle A\rangle\delta B+\delta A\langle B\rangle+\delta A\delta B
$$

**(4.2)**

$$
AB\simeq A\langle B\rangle+\langle A\rangle B-\langle A\rangle\langle B\rangle,
$$

**(4.3)**

where the product of fluctuations has been neglected. This is the only approximation in the algebraic reduction from a quartic interaction to a quadratic Hamiltonian.

## 4.2 Definition of the order parameter

Set

$$
\Delta(x,x')=V_{eff}(x,x')\langle\psi_\downarrow(x')\psi_\uparrow(x)\rangle,
$$

**(4.4)**

with

$$
\Delta^*(x,x')=V_{\mathrm{eff}}(x,x')\langle\psi_\uparrow^\dagger(x)\psi_\downarrow^\dagger(x')\rangle.
$$

**(4.5)**

Substitution gives

$$
H_{\mathrm{MF}}=H_0+\int dx\,dx'\left[\Delta(x,x')\psi_\uparrow^\dagger(x)\psi_\downarrow^\dagger(x')+\Delta^*(x,x')\psi_\downarrow(x')\psi_\uparrow(x)\right]
$$

**(4.6)**

$$
H_{MF}=H_0+\int dx\,dx'\frac{|\Delta(x,x')|^2}{V_{eff}(x,x')}.
$$

**(4.7)**

The final term is a c-number. It is essential when comparing thermodynamic energies, but it drops out of the operator equations of motion and of the BdG eigenvectors.

## 4.3 Why the result is quadratic

Every operator-dependent term now contains exactly two field operators. Consequently the Hamiltonian can be diagonalized by a linear canonical transformation. The price is that $\Delta$ must later satisfy the expectation-value definition used to introduce it.

# 5 Equations of motion

## 5.1 Fermionic commutator identity

For fermionic operators one repeatedly uses

$$
[A,BC]=\{A,B\}C-B\{A,C\}.
$$

**(5.1)**

For the normal term,

$$
[\psi_\uparrow(x),H_0]=\int dx'\,[\psi_\uparrow(x),\psi_\alpha^\dagger(x')h_{\alpha\beta}\psi_\beta(x')]
$$

**(5.2)**

$$
[\psi_\uparrow(x),H_0]=h_{\uparrow\beta}\psi_\beta(x).
$$

**(5.3)**

For the pair-creation term,

$$
\left[\psi_\uparrow(x),\int dx_1\,dx_2\,\Delta(x_1,x_2)\psi_\uparrow^\dagger(x_1)\psi_\downarrow^\dagger(x_2)\right]
$$

**(5.4)**

$$
\left[\psi_\uparrow(x),\int dx_1\,dx_2\,\Delta(x_1,x_2)\psi_\uparrow^\dagger(x_1)\psi_\downarrow^\dagger(x_2)\right]=\int dx_2\,\Delta(x,x_2)\psi_\downarrow^\dagger(x_2).
$$

**(5.5)**

The pair-annihilation term commutes with $\psi_\uparrow$ in the required ordinary commutator after the fermionic anticommutators are applied. Hence

$$
i\hbar\partial_t\psi_\uparrow(x)=h_{\uparrow\beta}\psi_\beta(x)+\int dx'\,\Delta(x,x')\psi_\downarrow^\dagger(x').
$$

**(5.6)**

Similarly,

$$
i\hbar\partial_t\psi_\downarrow(x)=h_{\downarrow\beta}\psi_\beta(x)-\int dx'\,\Delta(x',x)\psi_\uparrow^\dagger(x'),
$$

**(5.7)**

and Hermitian conjugation supplies the equations for the creation fields.

## 5.2 Compact matrix form before diagonalization

Introduce the field Nambu spinor

$$
\widehat{\Psi}(x)=\begin{pmatrix}\psi_\uparrow(x)\\\psi_\downarrow(x)\\\psi_\downarrow^\dagger(x)\\-\psi_\uparrow^\dagger(x)\end{pmatrix}.
$$

**(5.8)**

The chosen ordering is exactly the convention used by the code. It is related to other common Nambu bases by a unitary permutation and sign change; physical predictions are unchanged, but every matrix representation must remain consistent with the chosen basis.

# 6 Bogoliubov transformation

## 6.1 Canonical quasiparticles

A quasiparticle operator is sought as

$$
\gamma_n=\int dx\left[u_{n\uparrow}^\ast (x)\psi_\uparrow(x)+u_{n\downarrow}^\ast (x)\psi_\downarrow(x)+v_{n\downarrow}^\ast(x)\psi_\downarrow^\dagger(x)-v_{n\uparrow}^\ast(x)\psi_\uparrow^\dagger(x)\right].
$$

**(6.1)**

Its adjoint is obtained by conjugating coefficients and operators. Canonical anticommutation requires

$$
\{\gamma_n,\gamma_m^\dagger\}=\delta_{nm},
$$
$$
\{\gamma_n,\gamma_m\}=0,
$$

**(6.2)**

which implies the normalization

$$
\int dx\left(|u_{n\uparrow}|^2+|u_{n\downarrow}|^2+|v_{n\uparrow}|^2+|v_{n\downarrow}|^2\right)=1.
$$

**(6.3)**

## 6.2 Diagonalization condition

If

$$
H_{\mathrm{MF}}=E_{\mathrm{GS}}+\sum_n E_n\gamma_n^\dagger\gamma_n,
$$

**(6.4)**

then

$$
[\gamma_n,H_{\mathrm{MF}}]=E_n\gamma_n.
$$

**(6.5)**

Substituting the transformation and matching the independent field operators yields the BdG equations. This is equivalent to inserting the quasiparticle expansion into the Heisenberg equations of motion.

# 7 Bogoliubov–de Gennes equations

## 7.1 Component equations

Before adding spin-mixing terms, the spin-diagonal form used in the microscopic derivation is

$$
h_\uparrow u_{n\uparrow}+\Delta v_{n\downarrow}=E_nu_{n\uparrow},
$$

**(7.1)**

$$
h_\downarrow u_{n\downarrow}-\Delta v_{n\uparrow}=E_nu_{n\downarrow},
$$

**(7.2)**

$$
\Delta^*u_{n\downarrow}-h_\uparrow^*v_{n\uparrow}=E_nv_{n\uparrow},
$$

**(7.3)**

$$
\Delta^*u_{n\uparrow}-h_\downarrow^*v_{n\downarrow}=E_nv_{n\downarrow}.
$$

**(7.4)**

The opposite signs in the two singlet pair couplings express the antisymmetric spin structure. Rashba and transverse Zeeman terms add off-diagonal spin entries to the one-body operator; the compact Nambu matrix below incorporates those entries without changing the singlet pairing convention.

## 7.2 Nambu matrix

Define

$$
\Phi_n(x)=\begin{pmatrix}u_{n\uparrow} \\ u_{n\downarrow} \\ v_{n\downarrow} \\ -v_{n\uparrow}\end{pmatrix}.
$$

**(7.5)**

Then

$$
H_{\mathrm{BdG}}\Phi_n=E_n\Phi_n.
$$

**(7.6)**

For local singlet pairing and a one-body spin matrix $h_0$, a basis-independent block form is

$$
H_{\mathrm{BdG}}=\begin{pmatrix}h_0&\Delta i\sigma_y \\ -\Delta^* i\sigma_y&-h_0^*\end{pmatrix}
$$

**(7.7)**

in the conventional $(\psi_\uparrow,\psi_\downarrow,\psi_\uparrow^\dagger,\psi_\downarrow^\dagger)$ ordering. After transforming to the code ordering, the pairing term becomes $Re\Delta\,\tau_x-Im\Delta\,\tau_y$.

## 7.3 Nambu doubling is a redundancy

The BdG matrix has twice as many components as the electronic single-particle problem. This does not double the number of independent physical excitations. Particle–hole symmetry maps every positive-energy eigenvector to a negative-energy partner. One may construct the quasiparticle spectrum from one member of each pair.

# 8 Particle–hole symmetry in the code basis

## 8.1 Antiunitary operator

In the basis $(u_\uparrow,u_\downarrow,v_\downarrow,-v_\uparrow)^T$, the local particle–hole operation is

$$
\mathcal{C}=\tau_y\sigma_yK,
$$

**(8.1)**

where $K$ denotes complex conjugation. It satisfies

$$
\mathcal{C}H_{\mathrm{BdG}}\mathcal{C}^{-1}=-H_{\mathrm{BdG}}.
$$

**(8.2)**

Thus, if $H_{\mathrm{BdG}}\Phi_n=E_n\Phi_n$, then

$$
H_{\mathrm{BdG}}(\mathcal{C}\Phi_n)=-E_n(\mathcal{C}\Phi_n).
$$

**(8.3)**

At exactly zero energy, a state can be chosen to satisfy $\mathcal{C}\Phi_0=e^{i\chi}\Phi_0$; after a phase choice, it is self-conjugate and represents a Majorana solution at the BdG level.

## 8.2 Matrix check

For $N$ sites the global unitary part is

$$
U_C=I_N\otimes(\tau_y\sigma_y).
$$

**(8.4)**

The script evaluates the relative residual

$$
\epsilon_C=\frac{\left\|U_CH^*U_C^\dagger+H\right\|}{\|H\|}.
$$

**(8.5)**

This is a structural test of basis consistency, not a topological invariant.

# 9 Derivation of the self-consistency equation

## 9.1 Inverse Bogoliubov transformation

The field operators can be expanded in a complete set of BdG modes. In the selected convention, the anomalous expectation value at one point is built from products of electron and hole amplitudes. Thermal quasiparticle averages are

$$
\langle\gamma_n^\dagger\gamma_m\rangle=\delta_{nm}f(E_n),\qquad \langle\gamma_n\gamma_m^\dagger\rangle=\delta_{nm}[1-f(E_n)],
$$

**(9.1)**

where

$$
f(E)=\frac{1}{e^{E/(k_BT)}+1}.
$$

**(9.2)**

Using $1-2f(E)=\tanh[E/(2k_BT)]$ converts the occupation factors into the form used in the code.

## 9.2 Spin-singlet anomalous amplitude

For the basis vector

$$
\Phi_n(i)=\begin{pmatrix}u_{n\uparrow}(i) \\ u_{n\downarrow}(i) \\ v_{n\downarrow}(i) \\ -v_{n\uparrow}(i)\end{pmatrix},
$$

**(9.3)**

the local singlet combination is

$$
F_n(i)=u_{n\uparrow}(i)v_{n\downarrow}^\ast (i)-u_{n\downarrow}(i)v_{n\uparrow}^\ast (i).
$$

**(9.4)**

Because the fourth stored component is $-v_{n\uparrow}$, the same quantity is written computationally as

$$
F_n(i)=u_{n\uparrow}(i)v_{n\downarrow}^\ast (i)+u_{n\downarrow}(i)[-v_{n\uparrow}(i)]^\ast .
$$

**(9.5)**

## 9.3 Positive-energy and all-energy forms

When the sum is restricted to positive energies,

$$
\Delta_i=V_{\mathrm{eff}}\sum_{|E_n|<E_c \\ E_n>0}F_n(i)\tanh\left(\frac{E_n}{2k_BT}\right).
$$

**(9.6)**

If both members of every $\pm E_n$ pair are included, the equivalent expression contains a factor one half:

$$
\Delta_i=\frac{V_{eff}}{2}\sum_{|E_n|<E_c}F_n(i)\tanh\left(\frac{E_n}{2k_BT}\right).
$$

**(9.7)**

This is the form implemented by `calcular_delta_autoconsistente`.

## 9.4 Energy cutoff

The effective attraction is assumed to act only in a finite energy window $|E_n|<E_c$. The cutoff is part of the effective model. Changing it while keeping $V_{eff}$ fixed changes the converged gap; in a microscopic treatment the interaction and cutoff must be calibrated together.

## 9.5 Gauge freedom and phase fixing

A global transformation $\psi\mapsto e^{i\theta/2}\psi$ changes $\Delta\mapsto e^{i\theta}\Delta$ but leaves gauge-invariant observables unchanged. During iteration, numerical eigenvectors carry arbitrary phases. The script removes the average phase of $\Delta_i$ and chooses its mean real part positive. This does not alter the spectrum; it selects a stable gauge representative.

# 10 Continuum model of the Rashba–Zeeman nanowire

## 10.1 Normal Hamiltonian

The continuum electronic Hamiltonian is

$$
h_0(x)=\left[\frac{p_x^2}{2m^\ast}-\mu+V(x)\right]\sigma_0+\frac{\alpha_R}{\hbar}p_x\sigma_y+E_Z\sigma_x.
$$

**(10.1)**

The terms have distinct roles: kinetic energy defines the dispersion; $\mu$ selects the occupied portion of the band; $V(x)$ permits inhomogeneity; Rashba coupling locks spin to momentum; and the Zeeman term breaks time-reversal symmetry and separates the spin-split branches.

## 10.2 BdG Hamiltonian

In the code basis,

$$
H_{BdG}(x)=\left[\frac{p_x^2}{2m^\ast}-\mu+V(x)\right]\tau_z+\frac{\alpha_R}{\hbar}p_x\sigma_y\tau_z+E_Z\sigma_x+Re\Delta(x)\tau_x-Im\Delta(x)\tau_y.
$$

**(10.2)**

The factor $\tau_z$ on normal kinetic and spin–orbit terms expresses the opposite sign with which particles and holes inherit the normal-state Hamiltonian. The Zeeman term has the displayed form in this Nambu convention.

## 10.3 Topological criterion in the continuum limit

At $k=0$ the Rashba term vanishes. The gap closes when

$$
E_Z^2=\mu^2+|\Delta|^2.
$$

**(10.3)**

For a homogeneous infinite wire, the low-density continuum criterion is therefore

$$
E_Z>\sqrt{\mu^2+|\Delta|^2}.
$$

**(10.4)**

In a self-consistent calculation $\Delta$ itself depends on $E_Z$, so the phase boundary is not obtained by inserting one immutable gap value. Finite length, a discrete lattice, and spatial inhomogeneity further distinguish the numerical near-zero-energy map from a strict bulk invariant.

# 11 Spatial discretization without omitted steps

## 11.1 Lattice points and Taylor expansions

Let $x_i=ia$ and $\Phi_i=\Phi(x_i)$. Taylor expansion gives

$$
\Phi_{i+1}=\Phi_i+a\partial_x\Phi_i+\frac{a^2}{2}\partial_x^2\Phi_i+\mathcal{O}(a^3),
$$

**(11.1)**

$$
\Phi_{i-1}=\Phi_i-a\partial_x\Phi_i+\frac{a^2}{2}\partial_x^2\Phi_i+\mathcal{O}(a^3).
$$

**(11.2)**

Subtracting and adding yield

$$
\partial_x\Phi_i\simeq\frac{\Phi_{i+1}-\Phi_{i-1}}{2a},\qquad \partial_x^2\Phi_i\simeq\frac{\Phi_{i+1}-2\Phi_i+\Phi_{i-1}}{a^2}.
$$

**(11.3)**

## 11.2 Discrete kinetic term

Since $p_x^2/(2m^*)=-\hbar^2\partial_x^2/(2m^*)$, define

$$
t=\frac{\hbar^2}{2m^*a^2}.
$$

**(11.4)**

Then

$$
-\frac{\hbar^2}{2m^*}\partial_x^2\Phi_i\simeq2t\Phi_i-t\Phi_{i+1}-t\Phi_{i-1}.
$$

**(11.5)**

The $2t$ contribution enters the on-site block and the two $-t$ contributions enter the nearest-neighbour blocks.

## 11.3 Discrete Rashba term

With $p_x=-i\hbar\partial_x$,

$$
\frac{\alpha_R}{\hbar}p_x\sigma_y\tau_z=-i\alpha_R\partial_x\sigma_y\tau_z.
$$

**(11.6)**

The sign written in the final matrix depends on the chosen direction of the hopping block. Define

$$
t_{\mathrm{so}}=\frac{\alpha_R}{2a}.
$$

**(11.7)**

The action on site $i$ is represented by opposite imaginary amplitudes to the right and left, ensuring Hermiticity.

## 11.4 Local and hopping blocks

Collecting all terms,

$$
h_i=(2t-\mu+V_i)\tau_z+E_Z\sigma_x+Re\Delta_i\tau_x-Im\Delta_i\tau_y,
$$

**(11.8)**

$$
H_{i,i+1}=-t\tau_z+it_{\mathrm{so}}\sigma_y\tau_z,\qquad H_{i+1,i}=H_{i,i+1}^\dagger.
$$

**(11.9)**

The complete open-chain matrix is block tridiagonal:

$$
H=\begin{pmatrix}h_0&T&0& \cdots \\ T^\dagger&h_1 & T & \cdots \\ 0&T^\dagger&h_2&\ddots \\ \vdots & \vdots & \ddots & \ddots \end{pmatrix},\qquad T=-t\tau_z+it_{\mathrm{so}}\sigma_y\tau_z.
$$

**(11.10)**

Each block is $4\times4$; the full matrix is $4N\times4N$.

# 12 Homogeneous bulk Hamiltonian and bands

## 12.1 Fourier transform

For uniform $\Delta$ and $V=0$, write $\Phi_i=N^{-1/2}\sum_k e^{iki}\Phi_k$. The two hopping directions produce

$$
Te^{ik}+T^\dagger e^{-ik}=-2t\cos k\,\tau_z-2t_{\mathrm{so}}\sin k\,\sigma_y\tau_z.
$$

**(12.1)**

Hence

$$
H(k)=\xi_k\tau_z+\alpha_k\sigma_y\tau_z+E_Z\sigma_x+\Delta_R\tau_x-\Delta_I\tau_y,
$$

**(12.2)**

with

$$
\xi_k=2t-\mu-2t\cos k,\qquad \alpha_k=-2t_{\mathrm{so}}\sin k.
$$

**(12.3)**

## 12.2 Analytical dispersion

Squaring and diagonalizing the remaining spin structure gives four bands $\pm E_-(k)$ and $\pm E_+(k)$, where

$$
E_\pm^2(k)=\xi_k^2+\alpha_k^2+E_Z^2+|\Delta|^2
$$

**(12.4)**

$$
E_\pm^2(k)=\pm2\sqrt{E_Z^2|\Delta|^2+\xi_k^2(E_Z^2+\alpha_k^2)}.
$$

**(12.5)**

The script evaluates this expression and independently diagonalizes the $4\times4$ matrix at each $k$. Agreement of both curves checks the analytical formula and the matrix convention.

## 12.3 Bulk gap

The positive bulk excitation gap is

$$
E_{\mathrm{gap}}^{\mathrm{bulk}}=\min_k E_-(k).
$$

**(12.6)**

A bulk topological transition requires this quantity to close and reopen. By contrast, $\min_n|E_n|$ in an open finite wire can be exponentially small because of hybridized end states even while the bulk remains gapped.

# 13 Observables derived from BdG eigenvectors

## 13.1 Total BdG density

For a normalized eigenvector,

$$
\rho_n(i)=|u_{n\uparrow}(i)|^2+|u_{n\downarrow}(i)|^2+|v_{n\downarrow}(i)|^2+|-v_{n\uparrow}(i)|^2.
$$

**(13.1)**

The spatial sum is normalized to one. This quantity shows whether a state is extended, bulk-localized, or concentrated at the boundaries.

## 13.2 Electron and hole sectors

Define

$$
\rho_n^{(e)}(i)=|u_{n\uparrow}|^2+|u_{n\downarrow}|^2,\qquad \rho_n^{(h)}(i)=|v_{n\downarrow}|^2+|v_{n\uparrow}|^2.
$$

**(13.2)**

A Majorana-like zero mode has balanced integrated electron and hole weights, although balance alone does not establish topology.

## 13.3 Electronic local density of states

Using only one representative from each positive-energy pair, the broadened electronic LDOS is

$$
\rho(i,E)=\sum_{E_n>0}\left[\rho_n^{(e)}(i)L_\eta(E-E_n)+\rho_n^{(h)}(i)L_\eta(E+E_n)\right],
$$

**(13.3)**

where

$$
L_\eta(x)=\frac{\eta}{\pi(x^2+\eta^2)}.
$$

**(13.4)**

This is the quantity visualized as a function of energy and position.

## 13.4 Majorana combinations

Given a resolved particle–hole pair $\Phi_{+E}$ and $\Phi_{-E}$ with aligned phases, define

$$
\Gamma_1=\frac{\Phi_{+E}+\Phi_{-E}}{\sqrt{2}},\qquad \Gamma_2=-i\frac{\Phi_{+E}-\Phi_{-E}}{\sqrt{2}}.
$$

**(13.5)**

In a long topological wire these combinations can localize predominantly at opposite ends. Their overlap produces the finite-size energy splitting.

## 13.5 Finite-length splitting

For overlapping end modes,

$$
E_{\mathrm{split}}(L)\sim Ae^{-L/\xi}\cos(k_FL+\phi).
$$

**(13.6)**

The exponential envelope reflects localization, while the cosine produces the nonmonotonic nodes visible in finite-length sweeps.

# 14 Complete equation-to-code map

| Theory object | Implementation object |
|---|---|
| $t=\hbar^2/(2m^*a^2)$ | `hopping` and `calcular_hoppings_discretos` |
| $t_{\mathrm{so}}=\alpha_R/(2a)$ | `hopping_spin_orbita` |
| $h_i$ | `bloque_local` in `construir_hamiltoniano_bdg` |
| $H_{i,i+1}$ | `bloque_hopping` |
| $H_{\mathrm{BdG}}\Phi_n=E_n\Phi_n$ | `diagonalizar_hamiltoniano` |
| Self-consistent gap equation | `calcular_delta_autoconsistente` |
| Fixed-point update | `diagonalizar_autoconsistente` |
| $\rho_n(i)$ | `densidad_estado` |
| $\rho_n^{(e)},\rho_n^{(h)}$ | `densidad_electron_y_hueco` |
| $\rho(i,E)$ | `calcular_ldos` |
| $H(k)$ and $E_\pm(k)$ | `bandas_bulk_por_diagonalizacion` and `bandas_bulk_analiticas` |
| $\min_k E_-(k)$ | `gap_bulk_minimo` |
| $\Gamma_1,\Gamma_2$ | `construir_majoranas_desde_pareja` |

# 15 Conclusions

The self-consistent BdG model begins with an interacting spinful Hamiltonian and replaces the pairing interaction by an order parameter determined from an anomalous expectation value. The Bogoliubov transformation converts the quadratic mean-field Hamiltonian into a Hermitian eigenproblem in Nambu space. The same eigenvectors then regenerate the pairing field, closing the nonlinear problem.

For the nanowire, central finite differences yield a $4N\times4N$ block-tridiagonal matrix containing kinetic hopping, Rashba hopping, Zeeman splitting, electrostatic potential, and a local complex pair potential. The resulting spectrum supports bulk, finite-size, local-density, particle–hole, and Majorana-decomposition diagnostics. The essential conceptual distinction is that self-consistency makes the superconducting gap an output that responds to the Zeeman field and interaction strength rather than a passive external constant.

# Bibliography

[1] P. G. de Gennes, *Superconductivity of Metals and Alloys*.

[2] M. Tinkham, *Introduction to Superconductivity*.

[3] J.-X. Zhu, *Bogoliubov-de Gennes Method and Its Applications*, Lecture Notes in Physics 924 (2016).

[4] R. M. Lutchyn, J. D. Sau, and S. Das Sarma, semiconductor nanowire proposal for topological superconductivity, *Physical Review Letters* 105 (2010).

[5] Y. Oreg, G. Refael, and F. von Oppen, helical-liquid nanowire proposal, *Physical Review Letters* 105 (2010).
