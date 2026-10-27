# Numerical methods and observables

## Dynamics

Spins occupy an `L × L` square lattice and take values `−1` or `+1`. With `J = k_B = 1`, the Hamiltonian is the negative sum of spin products over undirected nearest-neighbour bonds; each bond is counted once. The horizontal direction is periodic. The default vertical boundary consists of fixed top (`−1`) and bottom (`+1`) rows, including their bond energies. Vertical periodic boundaries are available through `--periodic-y 1`.

At each attempted move, a mobile site is selected uniformly, followed by one of four directions uniformly. Moves that encounter a fixed row or select equal spins are rejected. The remaining exchanges use `min(1, exp(−ΔE/T))`; the local energy difference includes only bonds incident to the two sites. This preserves total magnetization exactly and maintains symmetric proposals for valid neighbouring mobile pairs.

A sweep comprises `L²` attempted moves, not accepted exchanges. The first `burn` sweeps are excluded from the statistics, and one configuration is sampled every `stride` sweeps thereafter.

## Initial state

For fixed vertical boundaries, only interior rows are initialized randomly while boundary spins remain fixed. The requested `m0` is converted to a feasible count of positive mobile spins; the actual global magnetization is reported as `m0_realized`. For fully periodic boundaries, all spins may be exchanged.

## Computed quantities

- **Global magnetization:** sum of spins divided by `L²`, unchanged by all exchanges.
- **Domain magnetizations:** average spin in the top and bottom regions. For nonzero `m0_realized`, the bottom region is assigned approximately `(1 + m0_realized)/2` of the lattice rows, rounded to a valid integer. For zero magnetization, the regions have equal height.
- **Domain order parameter:** weighted average of the absolute domain magnetizations, using the corresponding row fractions.
- **Mean energy per site:** `⟨E⟩/L²`.
- **Specific heat per site:** `[⟨E²⟩ − ⟨E⟩²]/(L²T²)`.
- **Domain response statistic:** `L² Var(q)/T`, where `q` denotes the domain order parameter. This is a descriptive fluctuation statistic, **not** the uniform magnetic susceptibility.
- **Density profile:** fraction of positive spins per row averaged over retained configurations.

## Statistical interpretation

- Sampling a Markov chain produces correlated observations; the standard error of repeated measurements along one chain is not automatically the error of an independent sample.
- Nonzero lattice dimensions and fixed boundaries shift the positions and shapes of response peaks. A finite-grid maximum should be interpreted as a pseudocritical temperature for those exact settings.
- A reliable limit as `L → ∞` requires additional sizes, convergence checks, uncertainty estimates and finite-size analysis.
- The code records all run parameters and the actual random seed to facilitate exact repetition on a compatible platform and compiler implementation. Random distributions are not guaranteed bitwise identical across different standard-library implementations.
