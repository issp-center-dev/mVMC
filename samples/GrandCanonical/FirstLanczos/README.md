---
date: 2026-10-16
datetime: 2026-10-16 23:58 JST
model: Claude Opus 5.5
summary: |
  Expert-mode input for the first-step stabilized power-Lanczos estimator
  with grand-canonical sampling (NGrandCanonical=1) and anti-parallel
  pairing on a two-site Hubbard ring. The same files are run by the
  GCLanczos_Sample test.
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
  including `-mu N`; no separate particle-number observable is reported.
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
example an exactly or nearly singular pair table, can therefore stop even
when ordinary grand-canonical sampling of the same input completes. The
first configuration of the final chain is redrawn at most 1024 times.
