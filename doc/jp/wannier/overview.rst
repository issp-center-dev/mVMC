概要
====

本資料では,
`RESPACK <https://sites.google.com/view/kazuma7k6r>`_ と
mVMC および :math:`{\mathcal H}\Phi` を用いて,
ダウンフォールディングをした格子モデルを計算する機能について説明する.

要件
----

`QuantumESPRESSO <http://www.quantum-espresso.org/>`_
もしくは
`xTAPP <http://xtapp.cp.is.s.u-tokyo.ac.jp/>`_
を用いてKohn-Sham軌道を用いたのちに,
RESPACKでWannier関数, 誘電関数, 有効相互作用を計算し,
それらを用いて構成した格子モデルを
mVMC もしくは :math:`{\mathcal H}\Phi`
で計算する.
したがってそれらのプログラムが使用可能である必要がある.
