# Model 1: Kitaev Chain

The Kitaev chain is the conceptual and numerical baseline of the repository. It isolates the minimal ingredients required for one-dimensional topological superconductivity: spinless fermions, nearest-neighbour hopping, odd-parity pairing and open boundaries.

The documentation proceeds from fermionic second quantization and Majorana operators to the exact real-space BdG matrix implemented in the code. It then derives the bulk dispersion, gap-closing conditions, class-D invariant, edge localization, finite-size splitting and every observable used in the figures.

## Files

- [Theory and numerical implementation](theory/THEORY.md) — the complete approved theory document followed by the complete implementation document, with equations written directly in Markdown mathematics.
- [Python source](code/kitaev_chain_1d.py) — the original model script.
- [Results and conclusions](results/RESULTS_AND_CONCLUSIONS.md) — all ten original figures placed next to their corresponding calculation, physical objective and interpretation.
- [`results/figures/`](results/figures/) — the ten original numerical graphs.

## Role in the progression

This model supplies exact reference signatures. The following BdG nanowire models can be assessed by asking which of these signatures survive after spin, Rashba coupling, Zeeman splitting, self-consistency, inhomogeneity and electrostatics are introduced.
