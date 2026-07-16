# Majorana BdG Nanowire Simulations

A progressive numerical study of Majorana bound states, from the ideal Kitaev chain to a spatially inhomogeneous semiconductor-superconductor nanowire with proximity renormalization and screened Hartree electrostatics.

> **Central evidentiary principle:** a state close to zero energy is a clue, not a proof. The interpretation becomes convincing only when several independent signatures agree: a topological bulk region, closing and reopening of the excitation gap, persistent near-zero finite-wire states, localization at opposite ends, Majorana-like spatial separation, particle-hole structure, finite-size behaviour and a compatible local density of states.

## Scientific progression

| Stage | Model | Physical advance | Documentation |
|---:|---|---|---|
| 1 | [Kitaev chain](models/01-kitaev-chain/) | Minimal spinless $p$-wave model; Majorana algebra, parity, bulk topology and boundary zero modes | [Theory and implementation](models/01-kitaev-chain/theory/THEORY.md) · [Results](models/01-kitaev-chain/results/RESULTS_AND_CONCLUSIONS.md) |
| 2 | [Self-consistent BdG nanowire](models/02-self-consistent-bdg/) | Spinful BdG formalism, Rashba-Zeeman nanowire and a self-consistent superconducting order parameter | [Theory and implementation](models/02-self-consistent-bdg/theory/THEORY.md) · [Results](models/02-self-consistent-bdg/results/RESULTS_AND_CONCLUSIONS.md) |
| 3 | [Realistic proximitized nanowire](models/03-realistic-proximitized-nanowire/) | Spatial proximity profiles, quasiparticle-weight renormalization, smooth electrostatics, screened Hartree response and sparse low-energy diagonalization | [Theory and implementation](models/03-realistic-proximitized-nanowire/theory/THEORY.md) · [Results](models/03-realistic-proximitized-nanowire/results/RESULTS_AND_CONCLUSIONS.md) |

Each stage is a controlled extension of the previous one. The Kitaev chain isolates the topological mechanism. The self-consistent BdG model introduces spinful superconducting quasiparticles and derives the nonlinear gap equation. The final model translates that framework into a more realistic effective description of a proximitized InSb nanowire.

## Repository structure

```text
.
├── README.md
├── requirements.txt
├── models/
│   ├── 01-kitaev-chain/
│   │   ├── README.md
│   │   ├── theory/THEORY.md
│   │   ├── code/kitaev_chain_1d.py
│   │   └── results/
│   │       ├── RESULTS_AND_CONCLUSIONS.md
│   │       └── figures/
│   ├── 02-self-consistent-bdg/
│   │   ├── README.md
│   │   ├── theory/THEORY.md
│   │   ├── code/self_consistent_bdg_nanowire.py
│   │   └── results/
│   │       ├── RESULTS_AND_CONCLUSIONS.md
│   │       └── figures/
│   └── 03-realistic-proximitized-nanowire/
│       ├── README.md
│       ├── theory/THEORY.md
│       ├── code/nanohilosimetricofinal.py
│       └── results/
│           ├── RESULTS_AND_CONCLUSIONS.md
│           └── figures/
├── supplementary/
│   ├── README.md
│   └── parametric-study-100-figures/
├── VALIDATION_REPORT.md
└── FILE_MANIFEST.json
```

The long-form theory files reproduce the approved documentation directly in GitHub-compatible Markdown. Equations are written as native mathematical source using inline `$...$` notation and display-math blocks delimited by double dollar signs; they are not screenshots. Only the original numerical graphs are embedded as images. The Python files are preserved as supplied.

## Core models

### Kitaev chain

The ideal spinless Hamiltonian is


$$
H=-\mu\sum_j c_j^\dagger c_j-t\sum_{j=1}^{N-1}\left(c_j^\dagger c_{j+1}+c_{j+1}^\dagger c_j\right)+\Delta\sum_{j=1}^{N-1}\left(c_j^\dagger c_{j+1}^\dagger+c_{j+1}c_j\right).
$$


It provides the exact benchmark for bulk bands, the topological transition at $|\mu|=2|t|$, finite-chain zero modes, Majorana localization and finite-size splitting.

### Self-consistent BdG nanowire

The second model derives the mean-field BdG problem and closes it through the gap equation. Its continuum Hamiltonian contains kinetic energy, chemical potential, electrostatic potential, Rashba coupling, Zeeman splitting and local singlet pairing.

### Realistic proximitized nanowire

The final model uses


$$
H_{\mathrm{BdG}}=Z(x)\left[\left(\frac{p_x^2}{2m^*}-\mu+V_{\mathrm{eff}}(x)\right)\tau_z+\frac{\alpha_R}{\hbar}p_x\sigma_y\tau_z+\Gamma\sigma_x\right]+\Delta_{\mathrm{ind}}(x)\tau_x,
$$


with


$$
V_{\mathrm{eff}}(x)=V_{\mathrm{ext}}(x)+U_H(x).
$$


The numerical evidence includes proximity and electrostatic profiles, Hartree convergence, bulk bands, Zeeman and chemical-potential sweeps, end-resolved spectral weight, particle-hole content and spatial Majorana decomposition.

## Installation

```bash
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

On Windows PowerShell, activate the environment with:

```powershell
.venv\Scripts\Activate.ps1
```

## Running the models

```bash
python models/01-kitaev-chain/code/kitaev_chain_1d.py
python models/02-self-consistent-bdg/code/self_consistent_bdg_nanowire.py
python models/03-realistic-proximitized-nanowire/code/nanohilosimetricofinal.py
```

The final realistic calculation is substantially more expensive than the pedagogical models. It uses separate spatial resolutions and sparse shift-invert diagonalization for low-energy states.

## Supplementary robustness study

The [supplementary parameter study](supplementary/README.md) contains one-parameter-at-a-time variations covering wire length, spatial resolution, chemical potential, Zeeman energy, Rashba coupling, effective mass, interface coupling, parent gap, electrostatic profiles, Hartree strength, screening, barriers, superconducting coverage, temperature, LDOS broadening and asymmetric interfaces.

## Preservation and validation

No PDF files are required to read the project: the approved long-form material is written directly into the Markdown documents. The source code and archived graphs have not been silently rewritten or regenerated. [`VALIDATION_REPORT.md`](VALIDATION_REPORT.md) records the structural, image-link, mathematical-markup and source-integrity checks, while [`FILE_MANIFEST.json`](FILE_MANIFEST.json) stores file sizes and SHA-256 hashes.
