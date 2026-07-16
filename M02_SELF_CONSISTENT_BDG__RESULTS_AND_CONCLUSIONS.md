# Self-Consistent BdG Nanowire

*Results, figure-by-figure interpretation, and physical conclusions using the original archived figures*

> Continuation of the Kitaev-chain documentation for the repository `majorana-bdg-nanowire-simulations`. The material starts from the second-quantized framework already established there and develops the Bogoliubov–de Gennes formalism and the self-consistent nanowire model.

Technical study document in English  
Source model: `src/self_consistent_bdg_nanowire.py`

# Result set and interpretation protocol

All figures in this document are copied from the repository’s `assets/figures/self_consistent` directory. They have not been recalculated, edited, or replaced. Each section contains two distinct explanations: first, the theoretical quantity and plotting procedure; second, the physical content visible in the archived result.

A consistent evidence chain is used throughout. A low energy alone is not treated as sufficient. The analysis compares convergence of the order parameter, bulk gap closing and reopening, finite-spectrum evolution, end localization, electron–hole balance, spatial Majorana decomposition, and finite-length splitting.

# Self-consistent order parameter

## Convergence of the gap iteration

![Convergence history of the self-consistent BdG loop.](figures/fig_selfconsistent_convergence.png)

*Convergence history of the self-consistent BdG loop.*

### What is calculated and why

The left panel plots $\langle|\Delta_i|\rangle$ after each mixed update. The right panel plots $\epsilon^{(m)}=\max_i|\Delta_i^{(m+1)}-\Delta_i^{(m)}|$ on a logarithmic scale. The first curve shows the physical fixed point approached by the pairing amplitude; the second provides the numerical stopping diagnostic. Both are necessary because a visually stable mean can conceal residual sitewise changes.

The data are produced from the histories stored by `diagonalizar_autoconsistente`. At each iteration the entire BdG matrix is rebuilt and diagonalized before a new anomalous pair amplitude is calculated. The plotted error is evaluated after linear mixing, so it measures the actual change passed to the next Hamiltonian.

### Result and interpretation

The mean gap decreases smoothly from the initial seed and approaches an approximately stationary value near $0.26$ meV. The update residual falls by several orders of magnitude and reaches the stated convergence scale after roughly eighty iterations. The absence of persistent oscillations indicates that the chosen mixing factor damps the fixed-point iteration effectively for this parameter set.

The two panels also show why the final gap must not be identified with the seed value. The self-consistent solution is determined by the spectrum and wavefunctions generated during the iteration. The seed only selects the starting point; the converged profile is the physical output within the mean-field model.

## Spatial pairing profile in three Zeeman regimes

![Magnitude of the converged pairing profile in trivial, intermediate, and topological regimes.](figures/fig_selfconsistent_gap_profile.png)

*Magnitude of the converged pairing profile in trivial, intermediate, and topological regimes.*

### What is calculated and why

For each representative Zeeman field, the figure plots $|\Delta_i|$ over the finite open chain. Because the gap equation is local, the profile can respond to boundaries and to changes in the quasiparticle spectrum. Comparing the three fields demonstrates directly that the order parameter is not a fixed external input.

The plotted quantity is the absolute value because the global condensate phase is gauge dependent and is fixed numerically. Spatial variation in the magnitude is gauge invariant. The central plateau approximates the homogeneous bulk value, whereas the end structure reflects the open boundary and finite-size quasiparticle amplitudes.

### Result and interpretation

The zero-field curve has the largest central gap, while the intermediate and topological curves are progressively suppressed. This is the expected qualitative competition between singlet pairing and Zeeman polarization: increasing $E_Z$ reduces the ability of opposite-spin electrons to form local Cooper pairs.

All three curves display enhanced and oscillatory structure near the ends before settling into a nearly uniform central value. The edge response is a genuine consequence of solving the gap locally in an open geometry. The topological case retains a nonzero central order parameter while supporting low-energy end states, showing that the topological regime is superconducting rather than simply gapless.

## Pairing profile versus interaction strength

![Self-consistent pairing profile for three effective attractive interactions.](figures/fig_selfconsistent_gap_vs_interaction.png)

*Self-consistent pairing profile for three effective attractive interactions.*

### What is calculated and why

The model is solved at the same nanowire and Zeeman parameters while $V_{\mathrm{eff}}$ is varied. Since $\Delta_i$ is proportional to $V_{\mathrm{eff}}$ times the anomalous expectation value, this sweep tests the nonlinear response of the condensate rather than imposing three chosen gaps by hand.

Every line is a separately converged fixed point. The comparison therefore includes both the direct prefactor change and the indirect change of the BdG eigenvectors caused by the altered pairing field.

### Result and interpretation

A larger attraction raises the central pairing magnitude throughout the wire. The curves remain ordered over almost the entire system, demonstrating a stable monotonic response within this parameter window. The boundary peaks also increase, so the interaction affects both the bulk plateau and the finite-size edge reconstruction.

The separation between the curves is substantial despite the relatively small change in $V_{\mathrm{eff}}$. This illustrates the nonlinear nature of the gap equation: the interaction modifies $\Delta$, the modified $\Delta$ changes the spectrum, and the new spectrum feeds back into the anomalous amplitude.

# Bulk spectral diagnostics

## Bulk bands in three regimes

![Homogeneous bulk bands obtained from the central self-consistent gap.](figures/fig_selfconsistent_bulk_bands.png)

*Homogeneous bulk bands obtained from the central self-consistent gap.*

### What is calculated and why

The central portion of each finite-wire gap profile is averaged to define a representative homogeneous $\Delta_{\rm bulk}$. The Bloch Hamiltonian $H(k)$ is then evaluated across the Brillouin zone. Solid curves are the analytical dispersion and dashed black curves are the direct eigenvalues of the $4\times4$ matrix.

This comparison serves two purposes. The analytical–numerical overlap checks the band formula and basis convention. The evolution of the lower positive band shows whether the bulk gap remains open, closes, and reopens as the Zeeman field is increased.

### Result and interpretation

In the trivial panel the positive and negative branches are separated by a clear gap. In the intermediate panel the inner branches approach one another near $k=0$, indicating proximity to the transition. In the topological panel a gap has reopened after the inversion of the low-energy branches.

The dashed numerical curves lie on top of the analytical bands over the entire momentum range. This agreement is a strong internal validation because the two calculations use different procedures. The sequence of gapped–closing–reopened spectra is the bulk signature required before finite-wire near-zero states can be interpreted topologically.

## Minimum bulk gap versus Zeeman energy

![Minimum positive bulk excitation energy versus Zeeman energy.](figures/fig_selfconsistent_bulk_gap_vs_zeeman.png)

*Minimum positive bulk excitation energy versus Zeeman energy.*

### What is calculated and why

At each $E_Z$, a self-consistent finite wire is solved, the central gap is extracted, and the homogeneous $H(k)$ is scanned to obtain $\min_kE_-(k)$. The resulting curve is the bulk excitation gap, not the energy of a finite-wire boundary state.

A topological phase transition in an infinite BdG system requires this gap to vanish and then reopen. Because $\Delta$ changes with field, the curve includes both the direct Zeeman effect on the band structure and the indirect suppression of the pairing field.

### Result and interpretation

The gap decreases from approximately $0.63$ meV at zero field toward a near-closure around the transition region. It then reopens to a finite value before decreasing more gradually at larger fields. This nonmonotonic behavior is the characteristic bulk transition pattern.

The reopened gap is smaller than the original zero-field gap. That asymmetry is physically meaningful in the self-consistent model: the Zeeman field both drives the band inversion and weakens the superconducting order. The topological protection scale is therefore set by the smaller reopened gap, not by the initial superconducting gap.

# Finite-wire spectral diagnostics

## Finite spectrum versus Zeeman energy

![Subgap finite-wire BdG levels as the Zeeman energy is swept.](figures/fig_selfconsistent_spectrum_vs_zeeman.png)

*Subgap finite-wire BdG levels as the Zeeman energy is swept.*

### What is calculated and why

For every field value, the full self-consistent $4N\times4N$ matrix is diagonalized and all levels within the displayed energy window are plotted. Particle–hole symmetry produces the visible reflection about zero. The dashed vertical line gives the continuum guide $E_Z\simeq\sqrt{\mu^2+\Delta^2}$ using a reference gap.

The finite spectrum complements the bulk gap. It reveals discrete confinement levels and the emergence of a pair of near-zero levels associated with hybridized end states. It should be interpreted together with spatial profiles because ordinary bound states can also approach zero.

### Result and interpretation

The low-energy branches move toward zero as the bulk transition is approached. Beyond the transition region, a pair remains extremely close to zero while higher excitations stay separated. This is the expected finite-wire manifestation of two end Majoranas whose overlap produces a small $\pm E$ splitting.

The spectrum is dense away from zero because the finite wire discretizes the bulk bands. The near-zero pair is distinguished not merely by its energy but by its persistence over a field interval and by the end-localized wavefunctions shown later. The plotted analytical guide is approximate because the pairing field itself changes during the sweep.

## Finite-wire map in the $(\mu,E_Z)$ plane

![Logarithm of the minimum absolute finite-wire energy over chemical potential and Zeeman energy.](figures/fig_selfconsistent_mu_zeeman_map.png)

*Logarithm of the minimum absolute finite-wire energy over chemical potential and Zeeman energy.*

### What is calculated and why

The nested parameter sweep stores $\min_n|E_n|$ at each grid point and displays its base-ten logarithm. Dark regions represent very small finite-wire energies. The logarithm makes exponentially small splittings visible over the same scale as ordinary subgap energies.

This is a computational phase-like diagnostic rather than a direct topological invariant. The expected continuum boundary grows with $|\mu|$, so the low-energy region should require a larger Zeeman field away from $\mu=0$.

### Result and interpretation

The darkest low-energy region appears at elevated Zeeman energy and is centered around small $|\mu|$. Its boundary bends upward as $|\mu|$ increases, reproducing the qualitative geometry of $E_Z^2\gtrsim\mu^2+|\Delta|^2$.

The region is structured rather than perfectly smooth because the plot comes from a finite chain, a coarse grid, and a self-consistent order parameter. Oscillatory finite-size splitting and discrete levels modulate $\min|E|$. The map is therefore most persuasive when read alongside the independently calculated bulk gap and spatial edge diagnostics.

# Local spectral weight

## LDOS in trivial, intermediate, and topological regimes

![Position- and energy-resolved local density of states.](figures/fig_selfconsistent_ldos_three_regimes.png)

*Position- and energy-resolved local density of states.*

### What is calculated and why

The electronic LDOS is built from positive-energy eigenstates with Lorentzian broadening. Electron amplitudes contribute at $+E_n$ and hole amplitudes at $-E_n$. The horizontal axis is position, the vertical axis energy, and color encodes spectral weight.

This observable links the global eigenvalue spectrum to a local probe. A bulk gap appears as an energy interval with little interior spectral weight. A topological end mode should add low-energy weight near the two boundaries rather than uniformly throughout the wire.

### Result and interpretation

The trivial panel shows a clear depletion around zero energy and strong spectral bands away from the gap. The intermediate panel shows the gap narrowing as the transition is approached. In the topological panel, low-energy spectral weight is associated with the boundaries while the interior retains a reopened gap.

The contrast between end and bulk weight is central. A purely extended gap closing would produce low-energy intensity throughout the wire. The coexistence of a gapped interior and boundary-concentrated near-zero weight supports the interpretation of localized edge excitations.

# Wavefunction and Majorana diagnostics

## State closest to zero: trivial versus topological

![Normalized density of the BdG state closest to zero in two regimes.](figures/fig_selfconsistent_mode_trivial_topological.png)

*Normalized density of the BdG state closest to zero in two regimes.*

### What is calculated and why

For each selected field, the eigenvalue with minimum absolute magnitude is identified. Its four Nambu–spin components are squared and summed site by site, then normalized. This plot tests whether the lowest state is an ordinary extended excitation or an end-localized subgap state.

The energy printed in each panel is as important as the density. A Majorana candidate requires both near-zero energy and boundary localization. Either property on its own is insufficient.

### Result and interpretation

In the trivial regime the nearest state lies at a finite energy of order the superconducting gap and extends across the wire with a standing-wave profile. It is therefore a conventional confined bulk quasiparticle.

In the topological regime the displayed energy is extremely close to zero and the density is concentrated at both ends with rapid oscillatory decay into the interior. The two-ended density is expected for one finite-energy BdG eigenstate formed from the overlap of left and right Majorana components.

## Electron–hole content of the near-zero state

![Electron and hole densities of the state closest to zero.](figures/fig_selfconsistent_particle_hole.png)

*Electron and hole densities of the state closest to zero.*

### What is calculated and why

The first two Nambu components define the electron density and the last two define the hole density. Both are normalized by the total BdG norm and plotted over position. Particle–hole symmetry predicts equal integrated weights for an ideal zero-energy Majorana state.

This test is complementary to spatial localization. An ordinary electron-like bound state may be localized but strongly unbalanced between particle and hole sectors. A self-conjugate state must contain both sectors coherently.

### Result and interpretation

The electron and hole curves nearly coincide at both ends and throughout their oscillatory tails. Their matching spatial structure indicates that the near-zero eigenstate is approximately charge neutral in the BdG sense.

The balance is consistent with Majorana self-conjugacy, but it is not by itself conclusive because particle–hole-related Andreev states can also become balanced near zero energy. Its significance comes from agreement with the bulk transition, end localization, and Majorana decomposition.

## Decomposition into left and right Majorana components

![Spatial densities of the two Majorana combinations formed from the near-zero pair.](figures/fig_selfconsistent_majoranas.png)

*Spatial densities of the two Majorana combinations formed from the near-zero pair.*

### What is calculated and why

The closest positive and negative BdG partners are phase aligned. Their symmetric and antisymmetric combinations define two self-conjugate Majorana wavefunctions. The code then compares the integrated density in the left half and labels the components by their dominant end.

A finite wire does not normally return left and right Majoranas directly as separate eigenvectors. The Hamiltonian eigenstates are the $\pm E$ fermionic combinations. Constructing $\Gamma_1$ and $\Gamma_2$ reveals whether that low-energy fermion can be resolved into spatially separated components.

### Result and interpretation

One component is sharply localized near the left boundary and the other near the right boundary. Both exhibit damped oscillations, as expected from a complex localization wave vector in a spin–orbit-coupled nanowire. Their densities are negligible over most of the central region.

The spatial separation explains the tiny but nonzero energy of the finite-wire state. The remaining overlap of the exponentially decaying tails couples the two components. Increasing the wire length reduces this coupling on average, as verified by the length sweep.

# Finite-size and interaction dependence

## Energy splitting versus wire length

![Minimum absolute energy versus the number of sites on a logarithmic scale.](figures/fig_selfconsistent_splitting_vs_length.png)

*Minimum absolute energy versus the number of sites on a logarithmic scale.*

### What is calculated and why

A separate self-consistent problem is solved for every wire length. The smallest absolute eigenvalue is plotted logarithmically. For spatially separated Majoranas, the expected asymptotic behavior is an exponentially decaying envelope multiplied by an oscillatory factor.

A logarithmic vertical scale is essential because the splitting spans many orders of magnitude. A monotonic straight line is not generally expected: the oscillatory phase $k_FL+\varphi$ produces nodes and revivals within the exponential envelope.

### Result and interpretation

The overall trend decreases strongly with length, from readily resolved subgap energies in short wires to values many orders of magnitude smaller in longer wires. Superimposed on this envelope are pronounced minima and partial revivals.

The nonmonotonic pattern is physically consistent with oscillatory Majorana hybridization rather than numerical evidence against exponential localization. Near a node of the cosine factor the splitting becomes exceptionally small; a modest increase in length can move away from the node and raise it again while the long-range envelope continues to fall.

## Interaction strength and edge localization

![Near-zero-state density for three attractive interaction strengths.](figures/fig_selfconsistent_interaction_localization.png)

*Near-zero-state density for three attractive interaction strengths.*

### What is calculated and why

The same topological-field calculation is repeated for three values of $V_{\mathrm{eff}}$. For each converged solution, the state closest to zero is selected and its normalized density is plotted. The legend reports the corresponding mean self-consistent gap.

The localization length depends inversely on the relevant topological excitation gap and on the effective velocity. Altering the pairing interaction changes both the condensate and the low-energy wavefunction, so this comparison exposes feedback between microscopic attraction and Majorana localization.

### Result and interpretation

All three states remain concentrated at the two ends, confirming that the qualitative edge-state structure survives across the tested interaction range. The mean gap increases from the weakest to the strongest attraction, in agreement with the separate gap-profile figure.

The decay and oscillatory amplitudes change with interaction strength. Stronger pairing generally confines the low-energy weight more tightly by increasing the superconducting energy scale, although the normalized peak heights also depend on the detailed oscillatory shape. The result demonstrates that localization is not controlled by Zeeman energy alone.

# Integrated conclusions

The archived calculations form a mutually consistent sequence. The gap iteration converges; the order parameter is suppressed by Zeeman field and enhanced by attraction; the homogeneous bulk gap closes and reopens; a finite-wire near-zero pair emerges after the transition; the corresponding density is concentrated at both ends; electron and hole weights are balanced; the pair can be decomposed into opposite-end Majorana components; and its splitting decreases with an oscillatory exponential envelope as the wire is lengthened.

The self-consistent treatment adds an important physical layer beyond a fixed-gap model. The same field that drives the topological band inversion also modifies the superconducting condensate. Consequently the reopened protection gap, the localization length, and the parameter boundary must be calculated from the converged order parameter rather than inferred from an unchanged input $\Delta$.

The results remain diagnostics of the specified mean-field lattice model. A finite-wire near-zero state is interpreted most reliably through the combined evidence above, not from any single plot in isolation.
