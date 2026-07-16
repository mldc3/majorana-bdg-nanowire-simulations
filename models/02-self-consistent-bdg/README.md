# Model 2: Self-Consistent BdG Nanowire

This section starts from the second-quantized framework established for the Kitaev chain and develops the Bogoliubov-de Gennes formalism for a spinful superconducting nanowire. The elementary Fock-space construction is not repeated. Instead, the interacting Hamiltonian is decoupled step by step, the Bogoliubov transformation is derived, the BdG equations are obtained, and the pairing field is closed through a self-consistency equation.

The model then specializes the formalism to a one-dimensional Rashba-Zeeman nanowire and discretizes the continuum operator into the block-tridiagonal matrix used by the Python implementation.

## Files

- [Theory and numerical implementation](theory/THEORY.md) — the complete approved BdG derivation followed by the complete implementation document and its unchanged source appendix.
- [Python source](code/self_consistent_bdg_nanowire.py) — the original self-consistent model script.
- [Results and conclusions](results/RESULTS_AND_CONCLUSIONS.md) — all thirteen original figures, each followed by what is calculated, why it is calculated and how the result is interpreted.
- [`results/figures/`](results/figures/) — the thirteen original numerical graphs.

## Role in the progression

The model replaces the ideal spinless pairing of the Kitaev chain with spinful quasiparticles, Rashba spin-orbit coupling, Zeeman splitting and a pairing profile determined together with the spectrum. The final repository model adds realistic proximity renormalization and electrostatic Hartree feedback.
