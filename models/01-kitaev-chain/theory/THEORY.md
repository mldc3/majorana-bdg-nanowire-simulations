# The Kitaev Chain and Majorana Zero Modes

*A step-by-step theoretical derivation aligned with the numerical repository*

Technical documentation for `majorana-bdg-nanowire-simulations`

July 2026

# Abstract

This document derives, without skipping the algebraic steps needed by the numerical implementation, the spinless one-dimensional Kitaev chain in real space and momentum space. The discussion begins with Fock space and fermionic second quantization, constructs Majorana operators and fermion parity, derives the general Bogoliubov–de Gennes (BdG) matrix from a quadratic superconducting Hamiltonian, and then specializes it to the exact convention used by `src/kitaev_chain_1d.py`. The infinite-chain dispersion, gap-closing conditions, topological invariant, finite-chain matrix, particle–hole symmetry, edge-mode localization, finite-size splitting, local density of states, and the numerical observables used in the plots are all obtained explicitly. The purpose is to make every matrix entry and every plotted quantity traceable to a physical operator and a stated convention.

# Contents

1. Scope, physical system, and conventions
   1. What system is being modelled?
   2. Units and notation
   3. The real-space Hamiltonian convention used by the code
2. Second quantization from first principles
   1. Why first quantization becomes inconvenient
   2. Fock space
   3. Canonical anticommutation relations
   4. Number operators
   5. From one-particle operators to second-quantized operators
3. Majorana operators and fermion parity
   1. Decomposing one complex fermion into two Hermitian operators
   2. Majorana algebra, derived explicitly
   3. Occupation expressed as a Majorana bilinear
   4. Fermion parity and physical even operators
      1. Local and global parity
      2. Action of parity on fermionic operators
      3. Why observables of an isolated fermionic system are even
      4. Parity conservation in a quadratic superconductor
   5. Non-local fermions and boundary Majoranas
4. Superconducting quadratic Hamiltonians and BdG form
   1. Why superconductivity mixes particles and holes
   2. Nambu doubling and the factor of one half
   3. Hermiticity of the BdG matrix
   4. Bogoliubov transformation
   5. Particle-hole symmetry
   6. The zero-energy Majorana condition
5. The Kitaev chain in real space
   1. Normal and pairing matrices
   2. Open boundaries
   3. Majorana representation of the chain
   4. Two limiting pictures
      1. Trivial atomic limit
      2. Topological “sweet spot”
6. Infinite chain: Fourier transformation and bulk bands
   1. Fourier transform
   2. Chemical potential and hopping
   3. Odd-parity pairing in momentum space
   4. Bloch BdG Hamiltonian
   5. Dispersion relation
   6. Gap-closing conditions
7. Topological invariant and bulk–boundary correspondence
   1. Class-D $\mathbb{Z}_2$ invariant
   2. Geometric winding picture for real parameters
   3. What bulk–boundary correspondence predicts
8. Analytical edge localization and finite-size splitting
   1. Zero-mode recurrence
   2. Finite-size energy splitting
9. Finite-chain observables used in the repository
   1. BdG probability density
   2. Electron and hole weights
   3. Constructing two Majorana combinations from a resolved $\pm E$ pair
   4. Bulk bands
   5. Finite spectrum versus chemical potential
   6. Minimum absolute energy
   7. Local density of states
   8. Splitting versus length
   9. Pairing strength and localization
   10. Phase diagram
10. Exact mapping from equations to the current Python model
11. Interpretive cautions
    1. A zero eigenvalue is necessary but not sufficient in realistic devices
    2. Finite-chain $E_{\min}$ is not the bulk gap
    3. Dense diagonalization and scaling
    4. Floating-point degeneracy
12. Conclusions

# Chapter 1

## Scope, physical system, and conventions

### 1.1 What system is being modelled?

The Kitaev chain is an idealized one-dimensional lattice of spinless fermionic modes. Each lattice site contains a single complex fermion operator $c_j$. Nearest-neighbour hopping allows a fermion to move between adjacent sites, the chemical potential controls the filling, and a nearest-neighbour odd-parity superconducting pair potential mixes particle and hole degrees of freedom. The model is not, by itself, a microscopic semiconductor nanowire. It is the minimal effective model that isolates the topology and boundary Majorana physics later reproduced by more realistic spinful nanowire Hamiltonians.

Two geometries are used and must not be confused:

1. An infinite or periodic chain, where translation symmetry permits a momentum-space Hamiltonian $H_{\mathrm{BdG}}(k)$. This geometry determines the bulk bands, the bulk gap, and the topological invariant.
2. A finite open chain of $N$ sites, where the bond between sites $N$ and $1$ is absent. This geometry can possess physical ends and therefore supports exponentially localized boundary modes.

> **Central physical question**
>
> For which values of the hopping $t$, chemical potential $\mu$, and nearest-neighbour pairing amplitude $\Delta$ is the bulk gapped and topologically non-trivial, and how does that bulk topology produce near-zero-energy states localized at the two ends of a finite chain?

### 1.2 Units and notation

The numerical script sets $t=1$ and therefore measures all energies in units of the hopping. Site labels are $j=1,\ldots,N$. For the bulk calculation, the lattice spacing is taken to be $a=1$, so the Brillouin zone is $k \in [-\pi,\pi]$.

The finite-chain Nambu basis used by the code is

$$
\Psi = \begin{pmatrix}c_1&\cdots&c_N&c_1^{\dagger}&\cdots&c_N^{\dagger}\end{pmatrix}^{T}.
$$

A BdG eigenvector therefore has the form

$$
\Phi_n = \begin{pmatrix}u_{n1}&\cdots&u_{nN}&v_{n1}&\cdots&v_{nN}\end{pmatrix}^{T},\qquad \sum_{j=1}^{N}\left(\lvert u_{nj}\rvert^2+\lvert v_{nj}\rvert^2\right)=1. 
$$

The first $N$ entries are electron amplitudes and the final $N$ entries are hole amplitudes.

### 1.3 The real-space Hamiltonian convention used by the code

The matrix implemented in the script corresponds, up to a physically irrelevant additive constant, to

$$
H=-\mu\sum_{j=1}^{N}c_j^{\dagger} c_j-t\sum_{j=1}^{N-1}\left(c_j^{\dagger} c_{j+1}+c_{j+1}^{\dagger} c_j\right)+\sum_{j=1}^{N-1}\left(\Delta c_j^{\dagger} c_{j+1}^{\dagger}+\Delta^{*}c_{j+1}c_j\right). 
$$

The more symmetric chemical-potential convention $-\mu\sum_j(c_j^{\dagger} c_j-1/2)$ differs only by the constant $N\mu/2$. Constants do not affect eigenvectors, excitation energies, topology, or any plot generated by diagonalizing the BdG matrix.

The ordering of operators in the pairing term fixes a phase convention for $\Delta$. Other texts often write $\Delta c_jc_{j+1}+\mathrm{H.c.}$. The two forms are related by a global superconducting phase and lead to the same spectrum and phase boundary. Throughout this document, signs are tracked using Eq. (1.3) so that the final matrix coincides with the code.

# Chapter 2

## Second quantization from first principles

### 2.1 Why first quantization becomes inconvenient

For a fixed number $N_p$ of identical fermions, a first-quantized wavefunction $\psi(x_1,\ldots,x_{N_p})$ must be antisymmetric under exchange of any two particle labels. The number of coordinates changes when particles are added or removed, and superconductivity explicitly couples sectors whose particle numbers differ by two. A formalism based on occupation numbers is therefore more natural.

Second quantization does not quantize a system “twice.” It changes the description. Instead of labelling indistinguishable particles, one specifies which one-particle modes are occupied. Creation and annihilation operators automatically enforce the exchange statistics.

### 2.2 Fock space

Let $\mathcal{H}_1$ be the one-particle Hilbert space spanned by orthonormal modes $\{\lvert i\rangle\}$. The fermionic Fock space is

$$
\mathcal{F}=\bigoplus_{N_p=0}^{\infty}\mathcal{H}_{N_p}^{(-)}, 
$$

where $\mathcal{H}_{N_p}^{(-)}$ denotes the antisymmetric $N_p$-particle sector. For $M$ available fermionic modes, each occupation number satisfies $n_i\in\{0,1\}$, and a basis state is

$$
\lvert n_1,n_2,\ldots,n_M\rangle=(c_1^{\dagger})^{n_1}(c_2^{\dagger})^{n_2}\cdots(c_M^{\dagger})^{n_M}\lvert0\rangle.
$$

The chosen order of the creation operators matters because exchanging two fermionic operators introduces a minus sign. There are two choices for every mode, so $\dim\mathcal{F}=2^M$.

### 2.3 Canonical anticommutation relations

Fermionic creation and annihilation operators satisfy

$$
\{c_i,c_j^{\dagger}\}=\delta_{ij},\qquad \{c_i,c_j\}=0,\qquad \{c_i^{\dagger},c_j^{\dagger}\}=0,
$$

where $\{A,B\}=AB+BA$. Setting $i=j$ in the last two relations gives

$$
2c_i^2=0,\qquad 2(c_i^{\dagger})^2=0,
$$

so

$$
c_i^2=(c_i^{\dagger})^2=0.
$$

Equation (2.5) is the operator form of Pauli exclusion: applying $c_i^{\dagger}$ twice to the same mode gives zero.

The local actions are

$$
c_i\lvert 0\rangle_i=0,\qquad c_i^{\dagger}\lvert 0\rangle_i=\lvert 1\rangle_i,\qquad c_i\lvert 1\rangle_i=\lvert 0\rangle_i,\qquad c_i^{\dagger}\lvert 1\rangle_i=0.
$$

For a many-mode occupation basis, additional signs appear because $c_i$ must be anticommuted through all earlier creation operators. Explicitly,

$$
c_i\lvert n_1\ldots n_M\rangle=(-1)^{\sum_{m=1}^{i-1}n_m}n_i\lvert n_1\ldots 0_i\ldots n_M\rangle.
$$

$$
c_i^{\dagger}\lvert n_1\ldots n_M\rangle=(-1)^{\sum_{m=1}^{i-1}n_m}(1-n_i)\lvert n_1\ldots 1_i\ldots n_M\rangle.
$$

### 2.4 Number operators

The occupation operator of mode $i$ is

$$
n_i=c_i^{\dagger} c_i. 
$$

Using $c_ic_i^{\dagger}=1-c_i^{\dagger} c_i$,

$$
n_i^2=c_i^{\dagger} c_ic_i^{\dagger} c_i=c_i^{\dagger}(1-c_i^{\dagger} c_i)c_i
$$

$$
=c_i^{\dagger} c_i-c_i^{\dagger} c_i^{\dagger} c_ic_i=c_i^{\dagger} c_i=n_i.
$$

Thus $n_i$ is a projector and its eigenvalues are $0$ and $1$. The total fermion number is

$$
\hat N_F=\sum_i n_i. 
$$

A normal hopping Hamiltonian commutes with $\hat N_F$. A superconducting mean-field Hamiltonian generally does not, because terms such as $c_i^{\dagger} c_j^{\dagger}$ change the particle number by two.

### 2.5 From one-particle operators to second-quantized operators

If a one-particle operator has matrix elements $h_{ij}=\langle i\lvert\hat h\rvert j\rangle$, its many-body second-quantized form is

$$
\hat H_1=\sum_{ij}h_{ij}c_i^{\dagger} c_j.
$$

To see why, apply $c_j$ to remove a fermion from mode $j$, apply the one-body amplitude $h_{ij}$, and then apply $c_i^{\dagger}$ to place it in mode $i$. Hermiticity requires $h=h^{\dagger}$.

A general two-body interaction is

$$
\hat H_2=\frac{1}{2}\sum_{ijkl}V_{ij;kl}c_i^{\dagger} c_j^{\dagger} c_lc_k.
$$

The factor $1/2$ prevents double counting of equivalent pairs. In BCS mean-field theory, a four-operator interaction is decoupled in a pairing channel and becomes an effective quadratic Hamiltonian containing $c^{\dagger} c^{\dagger}$ and $cc$ terms.

# Chapter 3

## Majorana operators and fermion parity

### 3.1 Decomposing one complex fermion into two Hermitian operators

For each complex mode $c_i$, define

$$
\gamma_{2i-1}=c_i+c_i^{\dagger},\qquad \gamma_{2i}=-i(c_i-c_i^{\dagger}).
$$

Taking Hermitian conjugates,

$$
\gamma_{2i-1}^{\dagger}=c_i^{\dagger}+c_i=\gamma_{2i-1},
$$

$$
\gamma_{2i}^{\dagger}=i(c_i^{\dagger}-c_i)=-i(c_i-c_i^{\dagger})=\gamma_{2i}.
$$

Thus each $\gamma_a$ is Hermitian. This is the condensed-matter analogue of self-conjugacy.

The inverse transformation follows by adding and subtracting Eq. (3.1):

$$
c_i=\frac{\gamma_{2i-1}+i\gamma_{2i}}{2},\qquad c_i^{\dagger}=\frac{\gamma_{2i-1}-i\gamma_{2i}}{2}.
$$

Two Majorana operators are therefore the two real components of one ordinary complex fermion.

### 3.2 Majorana algebra, derived explicitly

For the first Majorana associated with mode $i$,

$$
\gamma_{2i-1}^2=(c_i+c_i^{\dagger})^2 
$$

$$
=c_i^2+c_ic_i^{\dagger}+c_i^{\dagger} c_i+(c_i^{\dagger})^2 
$$

$$
=\{c_i,c_i^{\dagger}\}=1.
$$

For the second,

$$
\gamma_{2i}^2=-\left(c_i-c_i^{\dagger}\right)^2
$$

$$
=-c_i^2+c_ic_i^{\dagger}+c_i^{\dagger} c_i-(c_i^{\dagger})^2=1.
$$

Their mixed anticommutator is

$$
\{\gamma_{2i-1},\gamma_{2i}\}=-i \left[ (c_i+c_i^{\dagger})(c_i-c_i^{\dagger})+(c_i-c_i^{\dagger})(c_i+c_i^{\dagger}) \right] =0. 
$$

For different modes, the canonical anticommutation relations similarly give zero. Hence

$$
\{\gamma_a,\gamma_b\}=2\delta_{ab}.
$$

In particular, $\gamma_a^2=1$. This does not violate Pauli exclusion because a single $\gamma_a$ is not an annihilation operator for an independently occupiable complex fermion.

### 3.3 Occupation expressed as a Majorana bilinear

Substituting Eq. (3.4) into $n_i=c_i^{\dagger} c_i$,

$$
n_i=\frac{1}{4}(\gamma_{2i-1}-i\gamma_{2i})(\gamma_{2i-1}+i\gamma_{2i}) =\frac{1}{4}\left[\gamma_{2i-1}^2+\gamma_{2i}^2+i\gamma_{2i-1}\gamma_{2i}-i\gamma_{2i}\gamma_{2i-1}\right]
$$

$$
=\frac{1}{4}\left[2+2i\gamma_{2i-1}\gamma_{2i}\right] =\frac{1}{2}\left(1+i\gamma_{2i-1}\gamma_{2i}\right) \rightarrow i\gamma_{2i-1}\gamma_{2i}=2n_i-1.
$$

The occupation of the complex mode is encoded nonlinearly in a bilinear of its two Majoranas.

### 3.4 Fermion parity and physical even operators

#### 3.4.1 Local and global parity

The local parity of mode $i$ is

$$
P_i=(-1)^{n_i}=1-2n_i=-i\gamma_{2i-1}\gamma_{2i}.
$$

The equality $(-1)^{n_i}=1-2n_i$ follows because the only eigenvalues of $n_i$ are $0$ and $1$. Thus an empty mode has $P_i=+1$ and an occupied mode has $P_i=-1$.

The total parity operator is

$$
\hat P=(-1)^{\hat N_F}=e^{i\pi\hat N_F}=\prod_{i=1}^{M}P_i=\prod_{i=1}^{M}\left(-i\gamma_{2i-1}\gamma_{2i}\right). 
$$

With the natural ordered product of Majoranas, this may also be written

$$
\hat P=(-i)^M\gamma_1\gamma_2\cdots\gamma_{2M}.
$$

#### 3.4.2 Action of parity on fermionic operators

The total number operator satisfies

$$
[\hat N_F,c_i]=-c_i,\qquad [\hat N_F,c_i^{\dagger}]=c_i^{\dagger}.
$$

Using $e^{A}Be^{-A}=\sum_{m=0}^{\infty}\frac{1}{m!}\mathrm{ad}_{A}^{m}(B)$, with $\mathrm{ad}_{A}(B)=[A,B]$, one obtains

$$
\hat{P}\,c_i\,\hat{P}^{-1}=e^{i\pi\hat{N}_F}c_i e^{-i\pi\hat{N}_F}=e^{-i\pi}c_i=-c_i,
$$

$$
\hat{P}\,c_i^{\dagger}\,\hat{P}^{-1}=e^{+i\pi}c_i^{\dagger}=-c_i^{\dagger}.
$$

Every Majorana is likewise odd:

$$
\hat{P}\,\gamma_a\,\hat{P}^{-1}=-\gamma_a.
$$

An operator containing an odd number of fermionic factors changes sign under parity; a product containing an even number is parity invariant.

#### 3.4.3 Why observables of an isolated fermionic system are even

If global fermion parity is conserved, the Fock space separates into superselection sectors,

$$
\mathcal{F}=\mathcal{F}_{\mathrm{even}}\oplus\mathcal{F}_{\mathrm{odd}}.
$$

For $M$ complex fermionic modes, exactly half of the $2^M$ occupation states have even particle number and half have odd particle number, so

$$
\dim\mathcal{F}_{\mathrm{even}}=\dim\mathcal{F}_{\mathrm{odd}}=2^{M-1}. 
$$

Internal observables of a closed parity-conserving system must commute with $\hat P$. A single $c_i$, $c_i^{\dagger}$, or $\gamma_a$ connects opposite parity sectors and is therefore not by itself an observable of the isolated system. Even bilinears such as

$$
i\gamma_a\gamma_b 
$$

commute with parity and can appear in Hamiltonians and measurable response functions. This is why the occupation information of two Majoranas is carried by their bilinear rather than by either Majorana alone.

This statement does not prohibit single-electron tunnelling between a device and an external lead. The tunnelling operator may be odd when restricted to the device, but the full device-plus-reservoir operator is even and preserves the total parity of the enlarged closed system.

#### 3.4.4 Parity conservation in a quadratic superconductor

A superconducting mean-field Hamiltonian does not generally conserve $\hat N_F$, but every term contains an even number of fermionic operators. Explicitly,

$$
\hat P(c_i^{\dagger} c_j)\hat P^{-1}=(-c_i^{\dagger})(-c_j)=c_i^{\dagger} c_j,
$$

$$
\hat P(c_i^{\dagger} c_j^{\dagger})\hat P^{-1}=(-c_i^{\dagger})(-c_j^{\dagger})=c_i^{\dagger} c_j^{\dagger},
$$

$$
\hat P(c_jc_i)\hat P^{-1}=(-c_j)(-c_i)=c_jc_i. 
$$

Consequently,

$$
[H,\hat P]=0,\qquad [H,\hat N_F]\neq0\quad\text{in general}.
$$

Pairing changes the number by two and therefore leaves its parity unchanged. This distinction is essential in Majorana systems: the BdG description breaks particle-number conservation at mean-field level while retaining the parity structure that labels the two occupations of a non-local fermionic mode.

### 3.5 Non-local fermions and boundary Majoranas

Suppose two spatially separated Majorana operators $\gamma_L$ and $\gamma_R$ remain at low energy. They define a non-local complex fermion

$$
f=\frac{\gamma_L+i\gamma_R}{2},\qquad f^{\dagger}=\frac{\gamma_L-i\gamma_R}{2}, 
$$

with occupation

$$
f^{\dagger} f=\frac{1}{2}\left(1+i\gamma_L\gamma_R\right).
$$

No operator confined to only the left end can determine the full occupation without coupling to the right Majorana. This spatial separation is the origin of the non-local character often associated with topological Majorana modes. It does not make a finite device perfectly immune to all perturbations: any process that couples the two ends, changes global parity, or closes the bulk gap can remove the protection.

# Chapter 4

## Superconducting quadratic Hamiltonians and BdG form

### 4.1 Why superconductivity mixes particles and holes

A normal quadratic Hamiltonian contains $c_i^{\dagger}c_j$ and conserves particle number. A superconducting mean-field Hamiltonian also contains pair creation and annihilation:

$$
H=\sum_{i,j}c_i^{\dagger}h_{ij}c_j+\frac{1}{2}\sum_{i,j}\left(c_i^{\dagger}\Delta_{ij}c_j^{\dagger}+c_j\Delta_{ij}^{*}c_i\right).
$$

Hermiticity requires $h=h^{\dagger}$. Fermionic antisymmetry requires

$$
\Delta^{T}=-\Delta.
$$

To prove Eq. (4.2), note that only the antisymmetric part of a pair matrix contributes:

$$
\sum_{i,j}c_i^{\dagger}\Delta_{ij}c_j^{\dagger}
=
\frac{1}{2}\sum_{i,j}\left(c_i^{\dagger}\Delta_{ij}c_j^{\dagger}+c_j^{\dagger}\Delta_{ji}c_i^{\dagger}\right)
=
\frac{1}{2}\sum_{i,j}c_i^{\dagger}\left(\Delta_{ij}-\Delta_{ji}\right)c_j^{\dagger}.
$$

The symmetric part cancels because $c_j^{\dagger}c_i^{\dagger}=-c_i^{\dagger}c_j^{\dagger}$.

### 4.2 Nambu doubling and the factor of one half

Define the $2N$-component Nambu operator

$$
\Psi=
\begin{pmatrix}
c \\
c^{\dagger}
\end{pmatrix},
\qquad
c=
\begin{pmatrix}
c_1 & \cdots & c_N
\end{pmatrix}^{T}.
$$

Consider the matrix

$$
H_{\mathrm{BdG}}=
\begin{pmatrix}
h & \Delta \\
-\Delta^{*} & -h^{*}
\end{pmatrix}.
$$

Expanding $\frac{1}{2}\Psi^{\dagger}H_{\mathrm{BdG}}\Psi$ gives four terms:

$$
\frac{1}{2}\Psi^{\dagger}H_{\mathrm{BdG}}\Psi
=
\frac{1}{2}c^{\dagger}hc
+
\frac{1}{2}c^{\dagger}\Delta c^{\dagger}
-
\frac{1}{2}c\Delta^{*}c
-
\frac{1}{2}ch^{*}c^{\dagger}.
$$

The last normal term must be reordered. In indices,

$$
-\frac{1}{2}\sum_{i,j}c_i h_{ij}^{*}c_j^{\dagger}
=
-\frac{1}{2}\sum_{i,j}h_{ij}^{*}\left(\delta_{ij}-c_j^{\dagger}c_i\right).
$$

$$
=
-\frac{1}{2}\mathrm{Tr}\,h^{*}
+
\frac{1}{2}\sum_{i,j}c_j^{\dagger}h_{ij}^{*}c_i
=
-\frac{1}{2}\mathrm{Tr}\,h
+
\frac{1}{2}c^{\dagger}hc.
$$

where $h=h^{\dagger}$ was used in the last step. The pair-annihilation term similarly reproduces the Hermitian conjugate of the pair-creation term when $\Delta^{T}=-\Delta$. Therefore

$$
H=\frac{1}{2}\Psi^{\dagger}H_{\mathrm{BdG}}\Psi+\frac{1}{2}\mathrm{Tr}\,h.
$$

The factor $1/2$ is necessary because Nambu space contains both a degree of freedom and its Hermitian conjugate. The trace term is a constant and does not alter the BdG eigenproblem.

### 4.3 Hermiticity of the BdG matrix

Taking the adjoint of Eq. (4.6),

$$
H_{\mathrm{BdG}}^{\dagger}
=
\begin{pmatrix}
h^{\dagger} & \left(-\Delta^{*}\right)^{\dagger} \\
\Delta^{\dagger} & \left(-h^{*}\right)^{\dagger}
\end{pmatrix}
=
\begin{pmatrix}
h & -\Delta^{T} \\
\Delta^{\dagger} & -h^{*}
\end{pmatrix}.
$$

Since $\Delta^{T}=-\Delta$ and $\Delta^{\dagger}=\left(\Delta^{*}\right)^{T}=-\Delta^{*}$, one obtains $H_{\mathrm{BdG}}^{\dagger}=H_{\mathrm{BdG}}$.

### 4.4 Bogoliubov transformation

A quasiparticle annihilation operator is a linear combination

$$
\eta_n=\sum_j\left(u_{nj}^{*}c_j+v_{nj}^{*}c_j^{\dagger}\right).
$$

Its adjoint is

$$
\eta_n^{\dagger}=\sum_j\left(u_{nj}c_j^{\dagger}+v_{nj}c_j\right).
$$

Canonical fermionic anticommutation requires

$$
\left\{\eta_n,\eta_m^{\dagger}\right\}
=
\sum_j\left(u_{nj}^{*}u_{mj}+v_{nj}^{*}v_{mj}\right)
=
\delta_{nm},
$$

$$
\left\{\eta_n,\eta_m\right\}
=
\sum_j\left(u_{nj}^{*}v_{mj}^{*}+v_{nj}^{*}u_{mj}^{*}\right)
=
0.
$$

For one normalized mode this gives

$$
\sum_j\left(\lvert u_{nj}\rvert^2+\lvert v_{nj}\rvert^2\right)=1.
$$

Demanding $[\eta_n,H]=E_n\eta_n$ yields the BdG eigenvalue equation

$$
H_{\mathrm{BdG}}\Phi_n=E_n\Phi_n,
\qquad
\Phi_n=
\begin{pmatrix}
u_n \\
v_n
\end{pmatrix}.
$$

After collecting only one member of every $\pm E$ pair, the Hamiltonian can be written

$$
H=E_{\mathrm{GS}}+\sum_{E_n>0}E_n\eta_n^{\dagger}\eta_n.
$$

### 4.5 Particle-hole symmetry

In the basis of Eq. (1.1), define the antiunitary operator

$$
\mathcal{C}=\tau_xK,
\qquad
\tau_x=
\begin{pmatrix}
0 & I_N \\
I_N & 0
\end{pmatrix},
$$

where $K$ denotes complex conjugation. Acting on a vector,

$$
\mathcal{C}
\begin{pmatrix}
u \\
v
\end{pmatrix}
=
\begin{pmatrix}
v^{*} \\
u^{*}
\end{pmatrix}.
$$

The BdG matrix obeys

$$
\tau_xH_{\mathrm{BdG}}^{*}\tau_x=-H_{\mathrm{BdG}}.
$$

If $H_{\mathrm{BdG}}\Phi_E=E\Phi_E$

$$
H_{\mathrm{BdG}}(\mathcal{C}\Phi_E)=-E(\mathcal{C}\Phi_E).
$$

The spectrum is therefore symmetric about zero. Positive- and negative-energy eigenvectors are not independent physical quasiparticles; they are particle–hole partners generated by Nambu doubling.

### 4.6 The zero-energy Majorana condition

At exactly zero energy, particle–hole symmetry maps the zero-energy subspace into itself. A vector can be chosen to satisfy

$$
\mathcal{C}\Phi_0=e^{i\chi}\Phi_0. 
$$

A global phase $\Phi_0\mapsto e^{-i\chi/2}\Phi_0$ sets $\chi=0$, so

$$
\mathcal{C}\Phi_0=\Phi_0\quad\Longleftrightarrow\quad v=u^{*}.
$$

The corresponding operator satisfies

$$
\eta_0=\eta_0^{\dagger},
$$

which is precisely the Majorana self-conjugacy condition.

# Chapter 5

## The Kitaev chain in real space

### 5.1 Normal and pairing matrices

Comparing Eq. (1.3) with Eq. (4.1), the normal matrix is

$$
h_{ij}=-\mu\delta_{ij}-t(\delta_{i,j+1}+\delta_{i+1,j}), 
$$

and the pairing matrix is

$$
\Delta_{ij}=\Delta\delta_{i+1,j}-\Delta\delta_{i,j+1}.
$$

Equivalently,

$$
\Delta_{j,j+1}=+\Delta,\qquad \Delta_{j+1,j}=-\Delta.
$$

This explicitly verifies $\Delta^{T}=-\Delta$, as required for spinless odd-parity pairing.

For $N=4$ and real parameters,

$$
h=
\begin{pmatrix}
-\mu & -t & 0 & 0 \\
-t & -\mu & -t & 0 \\
0 & -t & -\mu & -t \\
0 & 0 & -t & -\mu
\end{pmatrix},
\qquad
\Delta=
\begin{pmatrix}
0 & \Delta & 0 & 0 \\
-\Delta & 0 & \Delta & 0 \\
0 & -\Delta & 0 & \Delta \\
0 & 0 & -\Delta & 0
\end{pmatrix}.
$$

The $8\times8$ BdG matrix is then

$$
H_{\mathrm{BdG}}=
\begin{pmatrix}
h & \Delta \\
-\Delta & -h
\end{pmatrix},
$$

which is exactly the block construction used by the script for real $t$, $\mu$, and $\Delta$.

### 5.2 Open boundaries

The loops in the code run only to $N-1$. Therefore the only hopping and pairing bonds are $(1,2),(2,3),\ldots,(N-1,N)$. No matrix entry connects $N$ back to $1$. This is an open chain. The omission is not a numerical convenience; it is the physical boundary required for end modes.

### 5.3 Majorana representation of the chain

Define at every site

$$
a_j=c_j+c_j^{\dagger},\qquad b_j=-i(c_j-c_j^{\dagger}). 
$$

Then $c_j=(a_j+ib_j)/2$. First, the onsite term is

$$
c_j^{\dagger} c_j=\frac{1}{4}(a_j-ib_j)(a_j+ib_j)
$$

$$
=\frac{1}{4}(a_j^2+b_j^2+ia_jb_j-ib_ja_j)
$$

$$
=\frac{1}{2}+\frac{i}{2}a_jb_j. 
$$

Thus

$$
-\mu\left(c_j^{\dagger} c_j-\frac{1}{2}\right)=-\frac{i\mu}{2}a_jb_j. 
$$

For neighbouring sites $j$ and $j+1$, direct substitution gives

$$
-t(c_j^{\dagger} c_{j+1}+c_{j+1}^{\dagger} c_j)=\frac{it}{2}(b_ja_{j+1}-a_jb_{j+1}), 
$$

and, for the code convention with real $\Delta$,

$$
\Delta(c_j^{\dagger} c_{j+1}^{\dagger}+c_{j+1}c_j)=-\frac{i\Delta}{2}(a_jb_{j+1}+b_ja_{j+1}). 
$$

Adding Eqs. (5.11) and (5.12),

$$
H=\frac{i}{2}\sum_{j=1}^{N}(-\mu a_jb_j)+\frac{i}{2}\sum_{j=1}^{N-1}\left[(t-\Delta)b_ja_{j+1}-(t+\Delta)a_jb_{j+1}\right]. 
$$

### 5.4 Two limiting pictures

#### 5.4.1 Trivial atomic limit

Set $t=\Delta=0$ and $\mu\neq0$. Then

$$
H=-\frac{i\mu}{2}\sum_j a_jb_j.
$$

The two Majoranas on each site are paired locally. There are no unpaired boundary operators. This is the trivial phase.

#### 5.4.2 Topological “sweet spot”

Set $\mu=0$ and $\Delta=t>0$. Equation (5.13) becomes

$$
H=-it\sum_{j=1}^{N-1}a_jb_{j+1}.
$$

Every $a_j$ for $j=1,\ldots,N-1$ is paired with $b_{j+1}$, but $b_1$ and $a_N$ do not appear in $H$. They commute with the Hamiltonian and are exact zero-energy Majorana operators localized on the first and last sites. The opposite edge labels often shown in textbooks result from the alternative superconducting phase convention; the physics is identical.

> **Bulk-boundary intuition**
>
> In the trivial limit, Majoranas are paired on the same site. In the topological limit, the energetically paired Majoranas lie on neighbouring sites, leaving one unpaired Majorana at each physical boundary. A continuous interpolation between the two patterns requires the bulk gap to close.

# Chapter 6

## Infinite chain: Fourier transformation and bulk bands

### 6.1 Fourier transform

For a periodic chain of $N$ sites,

$$
c_j=\frac{1}{\sqrt N}\sum_ke^{ikj}c_k,\qquad c_j^{\dagger}=\frac{1}{\sqrt N}\sum_ke^{-ikj}c_k^{\dagger}.
$$

The inverse follows from $\sum_je^{i(k-q)j}=N\delta_{kq}$.

### 6.2 Chemical potential and hopping

The onsite term is

$$
-\mu\sum_jc_j^{\dagger} c_j=-\frac{\mu}{N}\sum_{j,k,q}e^{-ikj}e^{iqj}c_k^{\dagger} c_q
$$

$$
=-\mu\sum_kc_k^{\dagger} c_k.
$$

For the forward hopping,

$$
\sum_jc_j^{\dagger} c_{j+1}=\frac{1}{N}\sum_{j,k,q}e^{-ikj}e^{iq(j+1)}c_k^{\dagger} c_q
$$

$$
=\sum_ke^{ik}c_k^{\dagger} c_k. 
$$

The Hermitian conjugate contributes $e^{-ik}$, so

$$
-t\sum_j(c_j^{\dagger} c_{j+1}+\mathrm{H.c.})=-2t\sum_k\cos k\,c_k^{\dagger} c_k.
$$

Hence the normal-state dispersion measured from the chemical potential is

$$
\xi_k=-\mu-2t\cos k. 
$$

### 6.3 Odd-parity pairing in momentum space

The creation part of the code’s pairing term is

$$
\sum_j\Delta c_j^{\dagger} c_{j+1}^{\dagger}=\frac{\Delta}{N}\sum_{j,k,q}e^{-ikj}e^{-iq(j+1)}c_k^{\dagger} c_q^{\dagger}
$$

$$
=\Delta\sum_ke^{ik}c_k^{\dagger} c_{-k}^{\dagger}.
$$

Because $c_{-k}^{\dagger} c_k^{\dagger}=-c_k^{\dagger} c_{-k}^{\dagger}$, only the odd part of $e^{ik}$ contributes when the sum is symmetrized over $k$ and $-k$:

$$
\sum_ke^{ik}c_k^{\dagger} c_{-k}^{\dagger}=\frac{1}{2}\sum_k\left(e^{ik}-e^{-ik}\right)c_k^{\dagger} c_{-k}^{\dagger}
$$

$$
=i\sum_k\sin k\,c_k^{\dagger} c_{-k}^{\dagger}.
$$

In the conventional $\frac{1}{2}\sum_k$ Nambu representation, the off-diagonal pair amplitude is therefore proportional to $2i\Delta\sin k$. It vanishes at $k=0$ and $k=\pi$, as required for odd-parity $p$-wave pairing.

### 6.4 Bloch BdG Hamiltonian

Use the momentum-space Nambu spinor

$$
\Psi_k=
\begin{pmatrix}
c_k \\
c_{-k}^{\dagger}
\end{pmatrix}.
$$

A matrix consistent with the real-space convention is

$$
H_{\mathrm{BdG}}(k)
=
\begin{pmatrix}
\xi_k & 2i\Delta\sin k \\
-2i\Delta\sin k & -\xi_k
\end{pmatrix}
=
\xi_k\tau_z-2\Delta\sin k\,\tau_y.
$$

The sign of the $\tau_y$ term changes under a global phase change of $\Delta$ and has no effect on the eigenvalues.

### 6.5 Dispersion relation

The characteristic equation is

$$
0
=
\det\!\left[H_{\mathrm{BdG}}(k)-EI_2\right]
=
\det
\begin{pmatrix}
\xi_k-E & 2i\Delta\sin k \\
-2i\Delta\sin k & -\xi_k-E
\end{pmatrix}
=
(\xi_k-E)(-\xi_k-E)-(2i\Delta\sin k)(-2i\Delta\sin k).
$$

$$
=
E^2-\xi_k^2-4\lvert\Delta\rvert^2\sin^2k
\rightarrow
E_{\pm}(k)
=
\pm\sqrt{(-\mu-2t\cos k)^2+(2\lvert\Delta\rvert\sin k)^2}.
$$

This is the exact formula implemented by `bandas_bulk_kitaev`.

### 6.6 Gap-closing conditions

A sum of non-negative terms can vanish only if both vanish:

$$
\xi_k=0,\qquad 2\Delta\sin k=0.
$$

For $\Delta\neq0$, $\sin k=0$ implies $k=0$ or $k=\pi$. At $k=0$,

$$
\xi_0=-\mu-2t=0\quad\Longrightarrow\quad\mu=-2t. 
$$

At $k=\pi$,

$$
\xi_\pi=-\mu+2t=0\quad\Longrightarrow\quad\mu=+2t. 
$$

Thus the superconducting bulk gap closes at


$$
\lvert\mu\rvert=2\lvert t\rvert. 
$$

If $\Delta=0$, the system is not a gapped topological superconductor. For $\lvert\mu\rvert<2\lvert t\rvert$ it is a gapless normal metal in the thermodynamic limit. This is why the code correctly excludes $\Delta=0$ from the topological region, although a binary phase map should label that line as gapless rather than ordinary trivial.

# Chapter 7

## Topological invariant and bulk–boundary correspondence

### 7.1 Class-D $\mathbb{Z}_2$ invariant

At the particle–hole-invariant momenta $k=0,\pi$, the pairing vanishes and $H_{\mathrm{BdG}}(k)=\xi_k\tau_z$. The one-dimensional class-D Majorana number can be expressed as the sign of the product of the two normal-state masses:

$$
\mathcal{M}=\mathrm{sgn}\!\left(\xi_0\xi_\pi\right).
$$

Using $\xi_0=-\mu-2t$ and $\xi_\pi=-\mu+2t$,

$$
\mathcal{M}
=
\mathrm{sgn}\!\left[(-\mu-2t)(-\mu+2t)\right]
=
\mathrm{sgn}\!\left(\mu^2-4t^2\right).
$$

Therefore

$$
\mathcal{M}
=
\begin{cases}
-1, & \lvert\mu\rvert<2\lvert t\rvert\quad\text{topological}, \\
+1, & \lvert\mu\rvert>2\lvert t\rvert\quad\text{trivial}.
\end{cases}
$$

The invariant can change only when one factor crosses zero, which is exactly the bulk gap closing found above.

### 7.2 Geometric winding picture for real parameters

For real $t$ and $\Delta$, Eq. (6.13) contains only $\tau_y$ and $\tau_z$. Define

$$
q(k)=\xi_k+i(2\Delta\sin k).
$$

As $k$ traverses the Brillouin zone, $q(k)$ traces an ellipse centered at $(-\mu,0)$ with horizontal semiaxis $2\lvert t\rvert$ and vertical semiaxis $2\lvert\Delta\rvert$. The winding number is

$$
\nu=\frac{1}{2\pi}\int_{-\pi}^{\pi}dk\,\partial_k\arg q(k).
$$

The origin lies inside the ellipse precisely when $\lvert\mu\rvert<2\lvert t\rvert$. Then $\lvert\nu\rvert=1$; otherwise $\nu=0$. The sign depends on the orientation, and hence on the sign and convention of $\Delta$. Its parity reproduces the class-D $\mathbb{Z}_2$ result.

### 7.3 What bulk–boundary correspondence predicts

For an interface between regions with different topological invariants, the bulk gap cannot remain open everywhere while the Hamiltonian changes continuously from one phase to the other. A localized boundary state must appear. A finite topological chain adjacent to vacuum has two interfaces, so it supports one Majorana mode near each end. In a finite chain, the two exponentially decaying wavefunctions overlap and combine into an ordinary fermionic state with small energies $\pm E_M$ rather than exactly zero.

# Chapter 8

## Analytical edge localization and finite-size splitting

### 8.1 Zero-mode recurrence

A zero-energy Majorana operator can be sought as a real linear combination of site Majoranas. For the code convention and $\mu=0$, one edge sector obeys a two-step recurrence whose characteristic magnitude is controlled by

$$
r=\frac{t-\Delta}{t+\Delta},\qquad 0<\Delta<t. 
$$

Because the recurrence advances by two sites, the envelope scales as

$$
\lvert\psi_j\rvert\propto\lvert r\rvert^{j/2}.
$$

Writing the envelope as $\exp(-j/\xi)$ gives

$$
e^{-2/\xi}=\left\lvert\frac{t-\Delta}{t+\Delta}\right\rvert, 
$$

so

$$
\xi=\frac{2}{\ln\left\lvert\frac{t+\Delta}{t-\Delta}\right\rvert}. 
$$

For $t=1$ and $\Delta=0.5$ gives

$$
\xi=\frac{2}{\ln3}\approx1.8205\ \text{sites}.
$$

This value will be recovered from the numerical splitting plot.

At $\Delta=t$, $r=0$ and the edge mode is confined to one site, reproducing the sweet spot. As $\Delta\to0^+$, $r\to1$ and $\xi\to\infty$: the topological gap becomes small and the Majorana delocalizes.

### 8.2 Finite-size energy splitting

Let the left and right Majorana envelopes scale as

$$
\psi_L(j)\sim e^{-j/\xi},\qquad \psi_R(j)\sim e^{-(N+1-j)/\xi}.
$$

Their overlap contains the factor $e^{-N/\xi}$. The effective low-energy Hamiltonian has the form

$$
H_{\mathrm{edge}}=\frac{i\varepsilon_M}{2}\gamma_L\gamma_R=\varepsilon_M\left(f^{\dagger} f-\frac{1}{2}\right),
$$

with

$$
\varepsilon_M\sim Ae^{-N/\xi}\cos(k_FN+\phi) 
$$

for generic parameters. The cosine represents possible oscillations from the underlying Fermi momentum. At the special $\mu=0$ parameters used in the repository, the decay is especially clean. Once $\varepsilon_M$ reaches roughly machine precision, numerical eigensolvers can no longer assign reliable signs or distinguish the two zero modes by energy alone.

> **Numerical consequence**
>
> A very small computed eigenvalue is not automatically a precisely resolved physical splitting. When $\lvert E_M\rvert$ is comparable to floating-point roundoff, the two-dimensional near-zero subspace is meaningful, but the individual eigenvectors and even the signs of the two eigenvalues are not. Majorana localization must then be extracted from the full near-zero subspace rather than by selecting “the smallest positive” and “the matching negative” eigenvalue using strict sign tests.

# Chapter 9

## Finite-chain observables used in the repository

### 9.1 BdG probability density

For a normalized eigenvector $\Phi_n=(u_n,v_n)^{T}$, the site-resolved BdG probability is

$$
p_n(j)=\lvert u_{nj}\rvert^2+\lvert v_{nj}\rvert^2,\qquad \sum_jp_n(j)=1. 
$$

This is the quantity returned by `densidad_estado`. It shows where the quasiparticle wavefunction is localized, but by itself it does not tell whether the state is topological, because ordinary Andreev bound states can also be localized.

### 9.2 Electron and hole weights

The separate components are

$$
p_n^{(e)}(j)=\lvert u_{nj}\rvert^2,\qquad p_n^{(h)}(j)=\lvert v_{nj}\rvert^2.
$$

The integrated electron and hole weights are

$$
W_e=\sum_j\lvert u_{nj}\rvert^2,\qquad W_h=\sum_j\lvert v_{nj}\rvert^2,\qquad W_e+W_h=1. 
$$

A useful BdG charge is

$$
Q_n=W_e-W_h=\Phi_n^{\dagger}\tau_z\Phi_n.
$$

For an ideal isolated Majorana zero mode, $W_e=W_h=1/2$ and $Q=0$. This equality must be computed with absolute squares. Using $u_j^2$ instead of $\lvert u_j\rvert^2$ is only accidentally harmless when the eigenvectors are real; it is mathematically incorrect for complex eigenvectors.

### 9.3 Constructing two Majorana combinations from a resolved $\pm E$ pair

Let $\Phi_E$ be a positive-energy state and $\Phi_{-E}=\mathcal{C}\Phi_E$ its particle–hole partner. Define

$$
\Phi_L=\frac{\Phi_E+\Phi_{-E}}{\sqrt2},\qquad \Phi_R=-\frac{i(\Phi_E-\Phi_{-E})}{\sqrt2}. 
$$

Applying particle–hole symmetry,

$$
\mathcal{C}\Phi_L=\frac{\Phi_{-E}+\Phi_E}{\sqrt2}=\Phi_L, 
$$

$$
\mathcal{C}\Phi_R=\frac{i(\Phi_{-E}-\Phi_E)}{\sqrt2}=\Phi_R. 
$$

Both combinations are self-conjugate. In a long topological chain, an appropriate phase convention makes one localized predominantly at the left and the other at the right.

When the splitting is unresolved, a robust alternative is:

1. take the two eigenvectors with smallest $\lvert E\rvert$;
2. form the projector onto that two-dimensional subspace;
3. diagonalize the position operator projected into the subspace to obtain maximally left- and right-localized combinations;
4. enforce particle–hole self-conjugacy within the subspace;
5. verify orthogonality and localization weights.

This procedure is stable even when both numerical eigenvalues are reported with the same sign or as exact zero.

### 9.4 Bulk bands

For each selected $\mu$, the code evaluates Eq. (6.18) on a uniform momentum grid. The plot should show a finite trivial gap at $\mu=3t$, a gap closing at one particle–hole-invariant momentum when $\mu=2t$, and a reopened topological gap at $\mu=0$.

### 9.5 Finite spectrum versus chemical potential

For every $\mu$ in a scan, the open-chain $2N\times2N$ matrix is diagonalized. Plotting all low-energy eigenvalues shows:

1. exact $\pm E$ spectral symmetry;
2. bulk-gap closing near $\mu=\pm2t$;
3. a near-zero pair in the topological interval $\lvert\mu\rvert<2t$;
4. finite-size deviations from the thermodynamic transition near the boundaries.

The near-zero pair is the fermionic combination of two overlapping edge Majoranas.

### 9.6 Minimum absolute energy

The diagnostic

$$
E_{\min}(\mu)=\min_n\lvert E_n(\mu)\rvert
$$

is easy to compute but must be interpreted carefully. In a finite topological chain it measures the edge-mode splitting, not the bulk gap. It becomes very small throughout the topological regime.

The bulk excitation gap should instead be estimated after excluding the near-zero edge pair or computed from the infinite-chain bands.

### 9.7 Local density of states

For exact discrete levels, the ideal spectral density is a sum of delta functions. Numerically, the script replaces each delta function with a Lorentzian

$$
L_\eta(E-E_n)=\frac{\eta}{\pi[(E-E_n)^2+\eta^2]},\qquad \int_{-\infty}^{\infty}L_\eta(E)\,dE=1.
$$

The plotted site-resolved LDOS is

$$
\rho_j(E)=\sum_np_n(j)L_\eta(E-E_n).
$$

The broadening $\eta$ is a visualization parameter, not a physical temperature or lifetime unless a specific experimental model justifies that interpretation.

In the topological regime, $\rho_j(E\approx0)$ is concentrated near the two ends. At the transition, the bulk gap closes and low-energy weight extends across the chain. In the trivial regime, there is no zero-energy end contribution.

### 9.8 Splitting versus length

For each chain length $N$, the matrix is rebuilt and diagonalized, and $\min\lvert E\rvert$ is recorded. In the regime where the physical splitting exceeds numerical roundoff, a semilog plot should be approximately linear:

$$
\ln E_{\min}(N)\approx\ln A-\frac{N}{\xi}. 
$$

The slope yields $\xi$. Data below the numerical precision floor must be excluded from the fit.

### 9.9 Pairing strength and localization

At fixed $t$ and $\mu$ in the topological phase, increasing $\lvert\Delta\rvert$ generally increases the superconducting gap and decreases the localization length. At $\mu=0$, Eq. (8.4) gives the dependence exactly for $0<\Delta<t$. The physically clean comparison is the density of one localized Majorana component, not an arbitrary eigenvector selected from a numerically degenerate pair.

### 9.10 Phase diagram

For the ideal homogeneous chain,

$$
\text{topological}\quad\Longleftrightarrow\quad\lvert\mu\rvert<2\lvert t\rvert\quad\text{and}\quad\Delta\neq0.
$$

The vertical phase boundaries do not depend on the magnitude of nonzero $\Delta$, but the size of the gap and the localization length do. The line $\Delta=0$ within $\lvert\mu\rvert<2\lvert t\rvert$ is gapless and should be visually distinguished from a gapped trivial phase.

# Chapter 10

## Exact mapping from equations to the current Python model

| Code object | Mathematical object | Physical meaning |
|---|---|---|
| `matriz_normal[j,j] = -mu` | $h_{jj}=-\mu$ | onsite chemical-potential energy |
| `matriz_normal[j,j+1] = -t` | $h_{j,j+1}=h_{j+1,j}=-t$ | nearest-neighbour hopping |
| `matriz_pairing[j,j+1] = Delta` | $\Delta_{j,j+1}=+\Delta$ | odd-parity pair creation amplitude |
| `matriz_pairing[j+1,j] = -Delta` | $\Delta_{j+1,j}=-\Delta$ | fermionic antisymmetry |
| `np.block(...)` | Eq. (4.6) | Nambu BdG Hamiltonian assembled from the four blocks |
| `la.eigh(H)` | $H_{\mathrm{BdG}}\Phi_n=E_n\Phi_n$ | Hermitian eigensystem |
| first $N$ vector entries | $u_{nj}$ | electron amplitudes |
| last $N$ vector entries | $v_{nj}$ | hole amplitudes |
| `abs(u)**2 + abs(v)**2` | Eq. (9.1) | site-resolved BdG probability |
| particle–hole transform | $(u,v)\mapsto(v^{*},u^{*})$ | action of $\mathcal{C}=\tau_xK$ |
| `min(abs(energias))` | Eq. (9.8) | lowest absolute finite-chain energy |
| Lorentzian sum | Eq. (9.10) | broadened LDOS |
| bulk routine | Eq. (6.18) | infinite-chain bands |
| phase test | Eq. (7.3) with $\Delta\neq0$ | ideal bulk topological criterion |

# Chapter 11

## Interpretive cautions

### 11.1 A zero eigenvalue is necessary but not sufficient in realistic devices

Within the ideal homogeneous Kitaev chain, a near-zero pair inside the topological interval has a clear interpretation. In realistic inhomogeneous nanowires, ordinary Andreev bound states can also approach zero. A convincing analysis combines bulk-gap closure and reopening, end localization, particle–hole balance, spatially separated Majorana components, length dependence, robustness to parameter variation, and an appropriate topological invariant.

### 11.2 Finite-chain $E_{\min}$ is not the bulk gap

Inside the topological phase, the smallest energy is the Majorana splitting. The bulk gap is the next positive excitation after excluding the near-zero pair. Confusing these quantities can lead to the incorrect statement that the topological phase is “gapless” because $E_{\min}\approx0$.

### 11.3 Dense diagonalization and scaling

A dense Hermitian eigendecomposition of a $2N\times2N$ matrix scales approximately as $\mathcal{O}(N^3)$ in time and $\mathcal{O}(N^2)$ in memory. Full spectra are appropriate for small pedagogical chains and LDOS calculations. For long chains or large parameter maps where only a few eigenvalues near zero are needed, sparse shift-invert methods are substantially more efficient.

### 11.4 Floating-point degeneracy

When the physical splitting is below about $10^{-14}t$-$10^{-16}t$ in double precision, the exact numerical threshold depends on matrix size, conditioning, and solver. Individual eigenvectors in the degenerate subspace are basis-dependent. Only subspace-invariant quantities, or quantities obtained after a physically motivated localization procedure, should be interpreted.

# Chapter 12

## Conclusions

The complete chain of reasoning is now explicit:

1. Fermionic second quantization introduces creation and annihilation operators obeying canonical anticommutation relations.
2. Every complex fermion can be decomposed into two Hermitian Majorana operators, and occupation is encoded in their bilinear.
3. Superconducting pairing breaks number conservation but preserves fermion parity and forces a Nambu particle–hole description.
4. The general quadratic Hamiltonian becomes the Hermitian BdG matrix $\begin{pmatrix}h&\Delta\\ -\Delta^{*}&-h^{*}\end{pmatrix}$.
5. For the Kitaev chain, nearest-neighbour odd-parity pairing makes $\Delta$ antisymmetric and produces the exact real-space matrix used in the Python script.
6. Fourier transformation yields $E_{\pm}(k)=\pm\sqrt{(-\mu-2t\cos k)^2+(2\Delta\sin k)^2}$.
7. The gap closes at $\mu=\pm2t$, and the invariant is non-trivial for $\lvert\mu\rvert<2\lvert t\rvert$ with $\Delta\neq0$.
8. Open boundaries then host two Majorana modes. Their overlap produces an exponentially small finite-size splitting and a near-zero BdG pair.
9. Every repository plot follows from a defined observable: bulk bands, finite spectrum, minimum absolute energy, LDOS, edge density, Majorana decomposition, particle–hole weights, length splitting, localization versus pairing, and the phase map.

The ideal Kitaev chain should therefore serve as the conceptual and numerical benchmark of the repository. More realistic nanowire models can be evaluated by asking which of these exact structural signatures survive when spin, spin–orbit coupling, Zeeman fields, inhomogeneity, electrostatic barriers, and induced superconductivity are added.

# Bibliography

1. A. Yu. Kitaev, “Unpaired Majorana fermions in quantum wires,” *Physics-Uspekhi* **44**, 131–136 (2001), arXiv:cond-mat/0010440.
2. J. Alicea, “New directions in the pursuit of Majorana fermions in solid state systems,” *Reports on Progress in Physics* **75**, 076501 (2012), arXiv:1202.1293.
3. M. Leijnse and K. Flensberg, “Introduction to topological superconductivity and Majorana fermions,” *Semiconductor Science and Technology* **27**, 124003 (2012), arXiv:1206.1736.
4. C. W. J. Beenakker, “Search for Majorana fermions in superconductors,” *Annual Review of Condensed Matter Physics* **4**, 113–136 (2013), arXiv:1112.1950.
5. A. L. Fetter and J. D. Walecka, *Quantum Theory of Many-Particle Systems*, McGraw–Hill (1971).
6. A. Altland and B. Simons, *Condensed Matter Field Theory*, 2nd ed., Cambridge University Press (2010).
7. P. G. de Gennes, *Superconductivity of Metals and Alloys*, Addison–Wesley (1966).
