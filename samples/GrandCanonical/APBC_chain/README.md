# Grand-canonical Hubbard chain with anti-periodic boundary conditions

Expert-mode input for `NGrandCanonical=1` with `NMPTrans=-1` on a four-site
Hubbard chain (`t=1`, `U=4`, `mu=2`, `phase0=180`).

- `stan.in`: the StdFace input that defines the lattice and Hamiltonian.
- `make_gc_apbc_input.py`: runs `vmcdry.out stan.in` and converts the result.
- `expert/`: the converted input, ready for
  `vmc.out -e namelist.def initial.def`.

Regenerate `expert/` with

```
python3 make_gc_apbc_input.py <path-to>/vmcdry.out expert
```

StdFace already writes the anti-periodic boundary phase into `trans.def`
(the hopping across the seam has the opposite sign) and one variational
parameter with its boundary sign per pair into `orbitalidxgen.def`. The
script leaves those files unchanged and only

- keeps `OrbitalGeneral` as the single orbital entry of `namelist.def`,
- adds `NGrandCanonical=1`, `NGCInitNelec=4` and `NSPGaussLeg=1` to
  `modpara.def` (`2Sz=-1`, `NMPTrans=-1` and `NExUpdatePath=0` are checked),
- writes explicit complex initial parameters to `initial.def`.

Grand-canonical sampling uses only the first `TransSym` pattern (the identity
map); projection over several translations is not performed.
`NMPTrans=-1` enables the sign column of `orbitalidxgen.def` but does not
change the Hamiltonian, so the boundary phase must be present in `trans.def`.
The test `GC_AP_StdFaceExample` checks that this directory is reproduced by
the script and that the signs agree with the lattice geometry.
