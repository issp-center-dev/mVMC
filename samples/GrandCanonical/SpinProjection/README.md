---
date: 2026-10-16
datetime: 2026-10-16 11:33 JST
model: OpenAI GPT-6 Codex
summary: |
  Four-site grand-canonical S=0 spin-projection example.
  Combined spin/momentum support, boundary zeros, and remaining restrictions.
---

# Grand-canonical S=0 spin projection

This expert-mode sample is a four-site Hubbard ring in the grand-canonical
ensemble (`NGrandCanonical=1`) with anti-parallel pairing
(`OrbitalAntiParallel`), `2Sz=0`, `NSPStot=0`, and one translation
(`abs(NMPTrans)=1`).  The Hamiltonian contains hopping, a chemical potential,
on-site and nearest-neighbour Coulomb terms, and an up-down pair field.  It
conserves `S^z` but is not `SU(2)` symmetric, which is allowed: the spin
projection acts on the trial state, not on the Hamiltonian.

The pair amplitudes are not a singlet, so the unprojected state mixes `S=0`,
`S=1`, and `S=2` components.  `NSPGaussLeg=G` with `G>=2` applies the
Gauss-Legendre spin projection

    Psi(x) = P_ch(x) sum_g a_g Pf[A(beta_g)|x],   a_g = (1/2) sin(beta_g) w_g

in every even particle-number sector.  `NSPGaussLeg=1` keeps the
unprojected state.  The vacuum amplitude is `sum_g a_g`, which tends to 1
only for large `G`; the sectors are not renormalized separately.  A finite
`G` is an approximation, so check convergence in `G`.

## Generate and run

`generate.py` writes all input files (it needs only Python 3):

```sh
python3 generate.py --gauss 8 --calmode 1 --samples 10000 --outdir g8
cd g8
OMP_NUM_THREADS=1 <build>/src/mVMC/vmc.out -e namelist.def initial.def
```

Options: `--gauss` (`NSPGaussLeg`), `--calmode` (0: one SR step with
`DSROptStepDt=0.002`, 1: Green functions), `--samples`, `--seed`,
`--warmup`, and `--boundary pbc|apbc`.  For APBC the boundary bond of `Trans`
carries the phase `-1`, the orbital rows carry the sign column, and
`NMPTrans=-1`.

## Convergence in G

Run the same input with `--gauss 1`, `--gauss 8`, and `--gauss 16` and
compare the energy, the particle number, and `S^2`.  `G=1` is the
unprojected state with `S^2 > 0`.  For this four-site input, `G=8` and
`G=16` agree within the Monte Carlo error, and `S^2` vanishes within that
error.  The particle number is not fixed, so `<N>` and `<N^2>` in `zvo_gc.dat`
are part of the comparison.

## S^2 from the Green functions

`greentwo.def` lists `S+_i S-_j` (`i 0 i 1 j 1 j 0`) and `S-_i S+_j`
(`i 1 i 0 j 0 j 1`) for all site pairs.  With the corresponding values
`G2` in `output/zvo_cisajscktalt_001.dat` and `<Sz^2>` (sixth column of
`output/zvo_out_001.dat`),

    S^2 = <Sz^2> + (1/2) sum_ij [ G2(S+_i S-_j) + G2(S-_i S+_j) ].

Every sample has `M=0`, so `<Sz^2>=0`.  Operators that change `S^z` as a
whole (the last row of `greentwo.def`, the same-spin rows of
`anomalousg.def`) are exactly zero.

## Combined spin and momentum projection

Spin projection also supports `abs(NMPTrans)>1`. Supply the corresponding
`TransSym` mappings and momentum weights; use signed mappings and a negative
`NMPTrans` for APBC. The generated sample above still uses one translation.
With `R=abs(NMPTrans)`, the combined state is

    Psi(x) = P_ch(x) sum_r sum_g w_r a_g Pf[A_r(beta_g)|x].

The same finite-G weights apply to every particle-number sector. If the
sum of the input translation weights is exactly zero, the vacuum and full
configuration have exactly zero amplitude. Start with
`0 < NGCInitNelec < 2*Nsite` in this case; either boundary fails explicitly
after one initialization attempt. Small nonzero sums are preserved. Exact
cancellation at intermediate particle numbers is not guaranteed.

`GC_SMP_` tests cover PBC/APBC, momentum 0 and pi/2, G=2/8, measurement and
SR modes, both `NStoreO` settings, MPI including empty QP ranks, and OMP.
They compare all sampled amplitudes and observables with an independent
finite-G oracle. The existing convergence tests remain applicable.

## Limits

Not supported together with the spin projection: `SpinJastrow` (it is not a
spin scalar and is rejected for `NSPGaussLeg>=2`), general pairing
(`OrbitalGeneral`), nonzero `NSPStot`, OptTrans, BackFlow,
RBM, and Lanczos.

`GC_SP_Sample_m0` and `GC_SP_Sample_m1` check that the generator output
matches the test fixture and replay a shortened run against an independent
oracle.
