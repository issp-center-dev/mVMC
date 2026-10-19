---
date: 2026-10-19
datetime: 2026-10-19 11:28 JST
model: OpenAI GPT-6 Codex (number-output documentation); Claude Opus 5.5 (original sample)
summary: |
  Expert-mode input for the first-step stabilized power-Lanczos estimator
  with grand-canonical sampling (NGrandCanonical=1) and anti-parallel
  pairing on a two-site Hubbard ring. The same files are run by the
  GCLanczos_Sample test. Final-state particle-number moments and block errors
  are written separately from the eleven energy columns.
---

# First-step Lanczos with grand-canonical sampling

This directory holds a complete expert-mode input that evaluates the energy
of the first Lanczos step

    |Phi> = (c0 + c1 H) |Psi_GC>

for a grand-canonical pair wave function `|Psi_GC>`. Run it with

```
vmc.out -e namelist.def initial.def
```

The run writes one result row to `output/zvo_pl_out_001.dat`. The two
header lines give the estimator version, `order=1`, the number of chains,
the samples per chain and in total, and the number of statistical blocks.
The eleven columns are

```
energy energy_se variance variance_se energy_imag variance_imag
tau_int effective_samples retained_rank condition gevp_residual
```

`energy_se` is the block standard error of the final chain for the frozen
coefficients `c0`, `c1`; it does not include the uncertainty of the
coefficients themselves. Standard output also reports the start of the final
chain (`GC Lanczos bootstrap: ...`), the energy line and `E_coeff`, the
energy of the coefficient chain for comparison. No Green functions and no
`zvo_gc.dat` are written in this mode.

The run also writes `output/zvo_gc_pl_number_001.dat`:

```
summary N_mean N_se N2_mean N2_se N_var total_samples
hist N count
block global_block_id count sum_N sum_N2
```

These are column names: the file contains one `summary` row, a `hist` row
for each N = 0, 1, 2, 3, 4, and one `block` row for each statistical block.
The header gives the frozen complex coefficients and sample/chain counts.
`N_se` and `N2_se` use the same equal-length blocks as energy and exclude
coefficient uncertainty. `N_var = N2_mean - N_mean*N_mean` has no error bar.
Only saved samples of the final chain contribute; warm-up and coefficient
samples are excluded. The histogram and block sums reproduce the moments.
For multiple chains the blocks partition all samples in chain-major order,
so a block may span a chain boundary; the histogram is not a time series.

For the exact optimal coefficients of this sample, independent Fock-space
summation gives N = 1.5206116258 and N² = 3.0413861424. A sampled run also
has uncertainty in its coefficients; its number estimates need not equal
these optimal-coefficient values. Assess mixing and block length before
interpreting the reported errors. The energy ESS is not an ESS for N.

## What the input describes

- Wave function: `OrbitalAntiParallel` with the complex 2x2 table
  `F = [[1, (1+2i)/8], [(3-i)/16, 3/4]]` (`orbitalidx.def`, values in
  `initial.def`). It has weight on N = 0, 2 and 4 with equal numbers of up
  and down electrons, so `2Sz 0` is required. No Gutzwiller or Jastrow
  factor and no projection (`NSPGaussLeg 1`, `NMPTrans 1`, the identity in
  `qptransidx.def`).
- Hamiltonian: `trans.def` holds the hopping `t=1` between the two sites
  (both directions, both spins) and the chemical potential as diagonal rows
  `i s i s 0.25`, that is `-mu n_is` with `mu=0.25`; `coulombintra.def`
  holds `U=1.5`. The reported energy is the expectation value of this H,
  including `-mu N`. The number file reports N and N² of the same final state;
  the program does not infer mu or convert the energy.
- Control: `NVMCCalMode 1`, `NLanczosMode 1`, `NLanczosStep 1`,
  `NLanczosEstimatorMode 1`, `NExUpdatePath 0`, 8192 samples per chain,
  512 warm-up sweeps. `NGCInitNelec` selects the electron number of the
  first configuration only.

For this small system the run gives an energy of about -1.35 with a
standard error of about 0.02. The exact energy of `|Phi>` for the optimal
coefficients is -1.3665, and the exact ground state of H is -1.8860.

## Limitations

The mode supports `Trans` (spin conserving and Hermitian), `CoulombIntra`,
`CoulombInter` and `Hund`, `OrbitalAntiParallel` (with spin and momentum
projection) or general pairing (momentum projection only), and complex
parameters. Green-function or anomalous measurements, other interaction
terms, `reweight`, LocSpin, BackFlow, RBM, OptTrans and Twist are rejected
when the input is read.

The run stops without writing a result row when a configuration proposed
after the start of the final chain has `Phi` that cannot be distinguished
from zero within its floating-point error (`GC Lanczos failure: ...
amplitude_state=numeric_zero` on standard error). Inputs with nodes, for
example a pair table that is nearly singular or singular only up to
rounding, can therefore stop even
when ordinary grand-canonical sampling of the same input completes. So can
momentum projection with nonzero total momentum, unless the `TransSym`
weights sum to exactly zero: otherwise the vacuum and fully occupied
amplitudes are only numerically zero. The first configuration of the final
chain is redrawn at most 1024 times.
