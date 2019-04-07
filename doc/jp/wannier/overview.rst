概要
====

本資料では,
`RESPACK <https://sites.google.com/view/kazuma7k6r>`_ と
mVMC および :math:`{\mathcal H}\Phi` を用いて,
ダウンフォールディングをした格子モデルを計算する機能について説明する.

.. math::

   \begin{aligned}
   {\cal H} &=
   \sum_{R, R', i, j, \sigma}
   \left(t_{(R'-R) i j} - t_{(R'-R) i j}^{\rm DC}\right)
   c_{R' j \sigma}^{\dagger} c_{R i \sigma}
   \nonumber \\
   &+ \sum_{R, i}
   U_{0 i j} n_{R i \uparrow} n_{R i \downarrow}
   + \sum_{(R, i) < (R', j)}
   U_{(R'-R) i j} n_{R i} n_{R' j}
   - \sum_{(R, i) < (R', j)}
   J_{(R'-R) i j} (n_{R i \uparrow} n_{R' j \uparrow}
   + n_{R i \downarrow} n_{R' j \downarrow})
   \nonumber \\
   &+ \sum_{(R, i) < (R', j)}
   J_{(R'-R) i j} (
   c_{R i \uparrow}^{\dagger} c_{R' j \downarrow}^{\dagger}
   c_{R i \downarrow} c_{R' j \uparrow} +
   c_{R' j \uparrow}^{\dagger} c_{R i \downarrow}^{\dagger}
   c_{R' j \downarrow} c_{R i \uparrow} )
   \nonumber \\
   &+ \sum_{(R, i) < (R', j)}
   J_{(R'-R) i j} (
   c_{R i \uparrow}^{\dagger} c_{R i \downarrow}^{\dagger}
   c_{R' j \downarrow} c_{R' j \uparrow} +
   c_{R' j \uparrow}^{\dagger} c_{R' j \downarrow}^{\dagger}
   c_{R i \downarrow} c_{R i \uparrow} ),
   \\
   t_{0 i i}^{\rm DC} &\equiv U_{0 i i} D_{0 i i}
   + \sum_{(R, j) (\neq 0, i)} (2 U_{R i j} - J_{R i j}) D_{0 j j},
   \\
   t_{R i j}^{\rm DC} &\equiv (J_{R i j} - U_{R i j}) D_{R i j}
   + 2 J_{R i j} {\rm Re} [D_{R i j}],
   \quad (R, j) \neq (0, i),
   \\
   D_{R i j} &\equiv \left\langle c_{R' j \sigma}^{\dagger} c_{0 i \sigma}\right\rangle.
   \end{aligned}

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
