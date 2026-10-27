---
date: 2026-10-27
datetime: 2026-10-27 13:43 JST
model: OpenAI GPT-6 Codex
summary: |
  Five-site S=2 first-step stabilized grand-canonical Lanczos input.
  Fixed complex anti-parallel parameters, PBC, G=16, and N=4/6 support.
  The public input is checked against an independent Fock-space replay.
---

# Finite total spin with the first Lanczos step

Copy this directory to a separate working directory and run:

```sh
export OMP_NUM_THREADS=1
vmc.out -e namelist.def initial.def
```

**This sample requires an enabled GMP/MPFR refinement backend to run.**
Install the GMP/MPFR development libraries and configure with
`-DMVMC_GC_REFINEMENT=ON`. Check that CMake prints
`GC exact-input refinement enabled (MPFR/GMP)`. The default `AUTO` enables
the backend only when those libraries are found; otherwise it leaves the
backend disabled. With `OFF`, or `AUTO` without the libraries, the input
is accepted but this sample stops in final warm-up with `numeric_zero`.
This occurred with the supplied seed and all eight alternative seeds
checked; backend-free execution is not a supported way to run this sample.

This is a five-site Hubbard ring with spin-independent nearest-neighbor
hopping t=1, U=1.5 and diagonal transfer 0.25 (the term −0.25 N in H).
The complex `OrbitalAntiParallel` pair function is projected with `2Sz=0`,
`NSPStot=2` (S, not 2S) and `NSPGaussLeg=16`. It has weight at N=4 and N=6.
The boundary is periodic and the translation list contains the identity.
There are no correlation factors or Green-function requests in this example.

The fixed parameters define the base state Ψ. Separate coefficient and final
chains evaluate Φ=(c₀+c₁H)Ψ with 8192 samples per chain and 512 warm-up sweeps.
`NGCInitNelec` is omitted, so the initial number is selected automatically
within the base support. To test initialization dependence, set it to 4 or 6
in a fresh copy. For projection parallelism, change `NSplitSize` to 2 and use
`mpiexec -np 2 vmc.out -e namelist.def initial.def`.

The run writes eleven energy/diagnostic columns to `output/zvo_pl_out_001.dat`
and number moments, histogram and block sums to `output/zvo_gc_pl_number_001.dat`.
Their format and conditional statistical errors are explained in the
[unprojected first-Lanczos example](../../FirstLanczos/README.md). The errors
do not include coefficient uncertainty; the energy ESS is not an ESS for N.
Check independent seeds, block length, initial N and longer sampling.

This H commutes with S², so it preserves an exactly projected S=2 state.
G=16 is still a finite-quadrature approximation: increase `NSPGaussLeg` to
check convergence. The implementation also accepts Sz-conserving H that
does not commute with S²; its final state need not retain the base S.
Structural zeros of Ψ are not imposed again on HΨ or Φ.

Refinement uses the existing binary64 inputs exactly and the same proposed
configuration and random value. It does not improve quadrature convergence.
It has finite memory, expression, precision and node-call limits; unresolved
or exhausted calculations stop without a result. See `NLanczosEstimatorMode`
in the expert manual. Large-system performance has not been established by
this small example.

## 日本語

5サイトHubbard環でS=2・N=4/6を持つ有限S＋1st stabilized Lanczosの入力です。
反平行複素pair、`2Sz=0`、周期境界、求積点数16を使います。`NSPStot` はSです。
別の作業場所へコピーし、上記コマンドを実行してください。
**このsampleの実行には、有効なGMP/MPFR精密backendが必要です。**
GMP/MPFR開発libraryを用意し、`-DMVMC_GC_REFINEMENT=ON` でconfigureしてください。
CMakeの出力に `GC exact-input refinement enabled (MPFR/GMP)` があることを確認します。
既定の `AUTO` はlibraryを検出した場合だけ有効化し、見つからなければ無効のままです。
`OFF`、またはlibraryのない `AUTO` では入力自体は受理されますが、このsampleは
final warm-upで `numeric_zero` により停止します。既定seedと変更した8 seedの
すべてで確認しており、backendなしで実行できるsampleではありません。

粒子数初期値は自動選択です。別runで `NGCInitNelec=4` または6を指定し、
独立seed、warm-up、採取期間、block長への依存を調べてください。粒子数の
平均・二乗・histogramはenergyと別fileへ出ます。報告誤差は係数の不確かさを
含まず、energyのESSを粒子数のESSとみなせません。

このHはS²を保存しますが、有限求積は厳密なスピン射影の近似なので点数の
収束確認が必要です。一般にはSz保存だけを要求するため、HがS²を保存しなければ
最終状態は同じSに留まりません。基底の構造零をHΨやΦに掛け直すこともしません。
精密backendは求積誤差を直さず、資源・精度上限に達すると停止します。
この小系sampleは大系での性能を保証するものではありません。
