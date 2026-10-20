---
date: 2026-10-20
datetime: 2026-10-20 14:58 JST
model: OpenAI GPT-6 Codex
summary: |
  Expert-mode ordinary GC examples for total spin S=1 and S=2, without Lanczos.
  Includes fixed-parameter measurement and dense-SR optimization inputs.
---

# Finite total spin in grand-canonical VMC

`S1` is a four-site ring with `NSPStot=1`. `S2` is a five-site ring with
`NSPStot=2`, with nonzero projected weight at both N=4 and N=6. Both use
anti-parallel complex pairing, `2Sz=0`, PBC, and `NSPGaussLeg=16`.
`NSPStot` denotes S, not 2S. Lanczos is disabled.

Copy a sample directory to a separate working directory, enter it, and run:

```sh
export OMP_NUM_THREADS=1
vmc.out -e namelist.def initial.def
```

The measurement input produces energy and particle-number moments, all
one-body Green functions, the two-body rows needed for S², and pair-add/remove
expectations. The Hamiltonian contains hopping, U=1.5, a diagonal transfer
of 0.25, and a complex up-down pair source. It conserves Sz; it need not
commute with S². The trial state is spin projected.

For a short dense-SR optimization, use a fresh copy:

```sh
vmc.out -e namelist_opt.def initial.def
```

The optimization input performs 20 steps. This is an execution example,
not a claim of optimization convergence. To measure the optimized state,
use `output/zqp_opt.dat` as the initial parameter file for the measurement
input in another working directory.

Both samples start from generic anti-parallel pair parameters with singlet
and triplet components. Do not replace them with an exactly pure triplet
state for unconstrained SR: derivatives on zero-amplitude configurations
can contribute to the energy gradient and overlap matrix but are absent
from the ordinary sampling estimator. This also affects S=0 and is not
automatically rejected. An arbitrarily small perturbation need not provide
adequate sampling. Fixed-parameter pure-triplet measurements are supported.

Increase `NSPGaussLeg` to check quadrature convergence. Check N distributions,
independent seeds, initial configurations and autocorrelation separately.
For S>0, explicit initial N must satisfy `2*S <= N <= 2*Nsite-2*S` and the
pair-function selection rules. `NGCInitNelec=-1` selects a permitted number.
The estimator reported as energy squared averages |E_local|² over nonzero
wave-function support; contributions of H Psi on nodes are not included.

For MPI, set `NSplitSize` to the QP group size and run, for example,
`mpiexec -np 2 vmc.out -e namelist.def initial.def` with `NSplitSize=2`.
G=1 omits spin projection. General pairing, SpinJastrow, BackFlow and
finite-S plus Lanczos are outside this sample's supported scope.

## 日本語

`S1` は4サイト・S=1、`S2` は5サイト・S=2の通常GC入力です。
後者はN=4とN=6の両方に重みを持ちます。`2Sz=0`、反平行複素pair、
周期境界、求積点数16を使用し、Lanczosは併用しません。
`NSPStot` はSを表します。

sampleを別の作業場所へコピーし、上記のコマンドで測定または20 stepの
dense SR最適化を実行してください。最適化の収束を保証する入力ではありません。
両sampleの初期対関数はsinglet成分とtriplet成分を持ちます。制約なしSRの初期値を
厳密な純tripletへ置き換えないでください。振幅零でも微分が非零の配置からの
勾配・重なり行列の寄与を通常SRは拾えません。S=0にもある制約で、自動拒否はしません。
任意に小さい摂動だけでは十分な標本化を保証できません。固定純triplet状態の測定は可能です。
求積点数、粒子数分布、独立seed、初期配置、自己相関をそれぞれ確認します。
`NGCInitNelec=-1` は許される初期粒子数を選びます。

HamiltonianはSzを保存する複素pair生成・消滅項を含み、S²保存は要求しません。
射影は試行状態に適用されます。出力のエネルギー二乗は波動関数の非零支持内の
|E_local|²平均で、節上のH Psiの寄与を含みません。
