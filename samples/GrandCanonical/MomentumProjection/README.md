---
date: 2026-10-16
datetime: 2026-10-16 13:12 JST
model: OpenAI GPT-6 Codex
summary: |
  Four-site grand-canonical momentum-projection example.
  Clarify exact-zero vacuum/full boundaries for momentum-only projection.
---

# Grand-canonical momentum projection

This expert-mode sample is a four-site periodic ring with anti-parallel
pairing and total momentum `K=pi`.  The four rows in `qptransidx.def` use
weights `(1, -1, 1, -1)/4`, and `NMPTrans=4` selects all four translations.
The pair function itself is not translation invariant, so the projection is
nontrivial.

The translation weights sum exactly to zero, so both the vacuum and the
fully occupied configuration have exactly zero projected amplitude. This
also applies with `NSPGaussLeg=1` (momentum projection without spin projection).
Keep `NGCInitNelec=4`, or another even value with `0 < NGCInitNelec < 8`;
starting at either boundary stops with a diagnostic after one attempt.
Small nonzero weight sums are not treated as zero, and exact cancellation
at intermediate particle numbers is not guaranteed.

From this directory, run with the path to a built executable:

```sh
OMP_NUM_THREADS=1 <build>/src/mVMC/vmc.out \
  -e namelist.def initial.def
```

The input uses 10,000 samples.  `GC_MP_Sample` checks that the committed
files match the test fixture generator and also runs a shortened smoke test.
