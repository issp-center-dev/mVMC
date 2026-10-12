# Grand-canonical momentum projection

This expert-mode sample is a four-site periodic ring with anti-parallel
pairing and total momentum `K=pi`.  The four rows in `qptransidx.def` use
weights `(1, -1, 1, -1)/4`, and `NMPTrans=4` selects all four translations.
The pair function itself is not translation invariant, so the projection is
nontrivial.

The `K=pi` vacuum has exactly zero projected amplitude.  Keep
`NGCInitNelec=4` (or another nonzero even value); an initial configuration
with zero projected amplitude stops with a diagnostic.

From this directory, run with the path to a built executable:

```sh
OMP_NUM_THREADS=1 <build>/src/mVMC/vmc.out \
  -e namelist.def initial.def
```

The input uses 10,000 samples.  `GC_MP_Sample` checks that the committed
files match the test fixture generator and also runs a shortened smoke test.
