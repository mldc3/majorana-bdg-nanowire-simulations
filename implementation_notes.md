# Implementation Notes

## Code architecture

The repository contains three Python scripts corresponding to three levels of modeling complexity.

### `kitaev_chain_1d.py`

This script implements the spinless Kitaev chain. It builds the normal matrix $h$, the antisymmetric p-wave pairing matrix $\Delta$, and the BdG block matrix

$$
H_{\mathrm{BdG}}=\begin{pmatrix}h&\Delta\\-\Delta^*&-h^*\end{pmatrix}.
$$

Main diagnostics:

- bulk bands,
- finite spectrum versus chemical potential,
- minimum absolute energy,
- LDOS,
- edge-mode density,
- particle-hole content,
- Majorana decomposition,
- splitting versus length,
- ideal phase diagram.

### `self_consistent_bdg_nanowire.py`

This script implements a spinful one-dimensional BdG nanowire with Rashba spin-orbit coupling, Zeeman field, and a pairing field computed self-consistently from the eigenvectors.

Main algorithm:

1. initialize $\Delta_i$;
2. build the BdG Hamiltonian;
3. diagonalize it;
4. update $\Delta_i$ from anomalous amplitudes;
5. mix the update with the previous profile;
6. iterate until convergence.

Main diagnostics:

- convergence of the self-consistent gap,
- spatial gap profile,
- bulk bands,
- finite spectrum versus Zeeman energy,
- LDOS,
- zero-mode density,
- Majorana decomposition,
- particle-hole balance,
- finite-size splitting,
- interaction-strength dependence.

### `proximitized_bdg_nanowire.py`

This script implements the more realistic proximitized Rashba nanowire. The gap is not computed self-consistently. Instead, the semiconductor is assigned an induced gap profile $\Delta_{\mathrm{ind}}(x)$ representing proximity to an external superconductor.

Main ingredients:

- effective mass $m^*$,
- Rashba coupling $\alpha_R$,
- Zeeman energy $E_Z=g\mu_BB/2$,
- induced gap $\Delta_{\mathrm{ind}}$,
- smooth superconducting coverage profile,
- Gaussian electrostatic barrier,
- finite-difference kinetic and spin-orbit terms.

Main diagnostics:

- bulk bands,
- finite spectra versus $E_Z$ and $\mu$,
- ideal phase diagram,
- minimum-energy maps,
- LDOS,
- Majorana decomposition,
- BdG charge,
- spin polarization,
- Rashba sweep,
- barrier sweep,
- mesh convergence.

## Numerical method

All scripts use explicit dense matrices and dense Hermitian diagonalization. This makes the implementation easy to inspect and mathematically transparent. The tradeoff is computational scaling. For a spinful BdG nanowire with $N$ sites, the matrix dimension is $4N\times4N$. Dense diagonalization scales approximately as $O(N^3)$ in time and $O(N^2)$ in memory.

For learning and figure generation, dense matrices are acceptable at moderate $N$. For large production simulations, the next improvement would be to use sparse matrices and target only eigenvalues near zero with shift-invert methods.

## Validation checks

The code includes or supports the following validation checks:

1. **Hermiticity**: the BdG Hamiltonian must satisfy $H=H^\dagger$.
2. **Particle-hole symmetry**: eigenvalues must appear in $\pm E$ pairs.
3. **Bulk/finite consistency**: finite-wire low-energy states should be compared to bulk gap diagnostics.
4. **Mesh convergence**: changing the spatial step should not qualitatively change the physical conclusion.
5. **Length scaling**: Majorana splitting should decrease as wire length increases.
6. **Phase alignment**: particle-hole partner eigenvectors must be phase-aligned before Majorana decomposition.

## Practical notes

The included scripts have been adjusted to use demo-safe default site counts. The selected figures in the repository were generated from the original larger calculations. To regenerate high-resolution results, increase the number of sites carefully and monitor memory use.

Recommended workflow:

```bash
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
python kitaev_chain_1d.py
python self_consistent_bdg_nanowire.py
python proximitized_bdg_nanowire.py
```

If a script becomes slow, reduce:

- `numero_sitios`,
- `numero_sitios_mapa`,
- number of Zeeman/chemical-potential sweep points,
- number of energy points used in LDOS plots.

## Suggested future improvements

A professional extension of the project could include:

- sparse diagonalization with `scipy.sparse.linalg.eigsh`,
- automated command-line arguments for all sweeps,
- unit tests for Hermiticity and particle-hole symmetry,
- cached data files for expensive sweeps,
- a notebook version explaining the derivations interactively,
- a comparison between topological Majoranas and quasi-Majoranas using identical diagnostic metrics,
- conductance calculations from a normal lead rather than only eigenstate diagnostics,
- disorder averaging and robustness analysis.
