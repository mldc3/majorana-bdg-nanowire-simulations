# Model 3: Realistic Proximitized Nanowire

This is the final and most detailed model in the repository. It carries the Kitaev-chain and BdG foundations into an effective one-dimensional description of an InSb semiconductor nanowire coupled to a conventional superconductor.

The Hamiltonian includes kinetic energy, chemical potential, Rashba spin-orbit coupling, Zeeman splitting, spatial superconducting coverage, induced pairing, quasiparticle-weight renormalization, smooth external electrostatics, screened Hartree response and open boundaries. The documentation derives the continuum model, topological criterion, proximity profiles, Hartree iteration, finite-difference discretization, bulk lattice Hamiltonian and numerical workflow.

## Files

- [Theory and numerical implementation](theory/THEORY.md) — the complete approved contents of documents 01, 02 and 03 in sequence, with all equations written directly as GitHub mathematics. The non-graph device schematic is intentionally omitted; the original numerical graphs appearing in the approved documents are retained.
- [Python source](code/nanohilosimetricofinal.py) — the original supplied implementation.
- [Results and conclusions](results/RESULTS_AND_CONCLUSIONS.md) — all fourteen original result graphs in the exact analytical sequence, each accompanied by its calculation, purpose and physical interpretation.
- [`results/figures/`](results/figures/) — the fourteen original numerical graphs.

## Role in the progression

This model is the repository's final physical case. A candidate Majorana mode is evaluated through the agreement of bulk, finite-size, spatial, particle-hole and spectral evidence rather than through a near-zero eigenvalue alone.
