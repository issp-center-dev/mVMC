非制限Hartree-Fock近似プログラム
================================

mVMCでは補助プログラムとして、多変数変分モンテカルロ法のペア軌道 :math:`f_{ij}` の初期値を
非制限Hartree-Fock(UHF)近似から与えるためのプログラムを用意しています(対応関係は :ref:`PuffAndSlater` を参照)。
なお、本プログラムは遍歴電子系を対象としており、スピン系、近藤系では正しく動作しません。

概要
----

UHF近似では揺らぎ :math:`\delta A \equiv A-\langle A \rangle` の一次までを考慮することで、二体項を一体項へと近似します。
たとえば、サイト間クーロン相互作用

.. math::

   {\cal H}_V = \sum_{i,j, \sigma, \sigma'}V_{ij} n_ {i\sigma}n_{j\sigma'}

について考えます。簡単化のため、 :math:`i\equiv (i, \sigma)`,
:math:`j\equiv (j, \sigma')` とすると相互作用の項は揺らぎの二次を落とすことで、

.. math::

   \begin{aligned}
   n_ {i}n_{j} &=& (\langle n_{i} \rangle +\delta n_i) (\langle n_{j} \rangle +\delta n_j) - \left[ \langle c_{i}^{\dag}c_j \rangle +\delta (c_{i}^{\dag}c_j ) \right] \left[ \langle c_{j}^{\dag}c_i \rangle +\delta (c_{j}^{\dag}c_i )\right] \nonumber\\
   &\sim&\langle n_{i} \rangle n_j+\langle n_{j} \rangle  n_i - \langle c_{i}^{\dag}c_j \rangle  c_{j}^{\dag}c_i  -  \langle c_{j}^{\dag}c_i \rangle c_{i}^{\dag}c_j 
   -\langle n_{i} \rangle \langle n_j \rangle +  \langle c_{j}^{\dag}c_i \rangle \langle c_{i}^{\dag}c_j \rangle
   \end{aligned}

と近似されます。このような形式で、その他の相互作用についても近似を行うことで、一体問題に帰着させることができます。
計算では、上記の各平均値がself-consistentになるまで計算を行います。

ソースコード
~~~~~~~~~~~~

ソースコード一式は ``src/ComplexUHF/src`` 以下に入っています。

コンパイル方法
~~~~~~~~~~~~~~

コンパイルはmVMCのコンパイルと同様にmVMCのフォルダ直下で

::

    $ make mvmc

を実行することで行われます。コンパイルが終了すると、
``src/ComplexUHF/src`` に実行ファイル ``UHF`` が作成されます。

必要な入力ファイル
~~~~~~~~~~~~~~~~~~

入力ファイル指定用ファイル (namelsit.def)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

| UHFで指定するファイルは以下のファイルです。
  ``namelist.def`` は :ref:`InputFileList` で定義されているファイルと同じ様式です。

-  ``ModPara``

-  ``LocSpin``

-  ``Trans``

-  ``CoulombIntra``

-  ``CoulombInter``

-  ``Hund``

-  ``PairHop``

-  ``Exchange``

-  ``Orbital`` / ``OrbitalAntiParallel``

-  ``OrbitalParallel``

-  ``OrbitalGeneral``

-  ``Initial``

基本的にはmVMCと同じファイルとなりますが、

-  ``ModPara`` ファイルで指定されるパラメータ

-  ``Initial`` ファイルの追加

がmVMCと異なります。以下、その詳細を記載します。

ModParaファイルで指定するパラメータ
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

UHFで指定するパラメータは以下のパラメータです。

-  ``Nsite``

-  ``Ne``

-  ``Mix``

-  ``EPS``

-  ``IterationMax``

``Nsite``,
``Ne`` はmVMCと共通のパラメータで、以下の三つがUHF独特のパラメータです。

-  | ``Mix``
   | linear
     mixingをdouble型で指定します。mix=1とすると完全に新しいGreen関数に置き換えられます。

-  | ``EPS``
   | 収束判定条件をint型で指定します。新しく計算されたGreen関数と一つ前のGreen関数の残差が :math:`10^{-\verb|eps|}` の場合に、計算が打ち切られます。

-  | ``IterationMax``
   | ループの最大数をint型で指定します。

なお、mVMCで使用するその他パラメータが存在する場合はWarningが標準出力されます(計算は中断せずに実行されます)。

Initialファイル
^^^^^^^^^^^^^^^

グリーン関数 :math:`G_{ij\sigma_1\sigma_2}\equiv \langle c_{i\sigma_1}^\dag c_{j\sigma_2}\rangle` の初期値を与えます。
ファイル様式は ``Trans`` ファイルと同じで、 :math:`t_{ij\sigma_1\sigma_2}` の代わりに :math:`G_{ij\sigma_1\sigma_2}` の値を記述します。
なお、値を指定しないグリーン関数には０が入ります。

使用方法
--------

UHF自体はmVMCと同じように

::

    $ UHF namelist.def

で動きます。計算の流れは以下の通りです。

#. ファイル読み込み

#. ハミルトニアンの作成

#. グリーン関数の計算 (self-consistentになるまで)

#. :math:`f_{ij}` 、各種ファイルの出力

計算後に出力されるファイルおよび出力例は以下の通りです。

-  | zvo\_result.dat: エネルギーと粒子数が出力されます。

   ::

        energy -15.2265348135
        num    36.0000000000

-  zvo\_check.dat:
   イタレーションのステップ数、グリーン関数の残差の絶対値の平均、収束過程のエネルギー、粒子数を順に出力します。

   ::

        0  0.004925645652 -544.963484605164 36.000000
        1  0.002481594941 -278.304285708488 36.000000
        2  0.001274395448 -147.247026925130 36.000000
        3  0.000681060599 -82.973664527606 36.000000
       ...

-  | zvo\_UHF\_cisajs.dat:
     収束した一体グリーン関数 :math:`G_{ij\sigma_1\sigma_2}\equiv\langle c_{i\sigma_1}^{\dag}c_{j\sigma_2}\rangle` 。
   | 全成分について :math:`i, \sigma_1, j, \sigma_2, {\rm Re}\left[G_{ij\sigma_1\sigma_2}\right], {\rm Im}\left[G_{ij\sigma_1\sigma_2}\right]` の順に出力されます。

   ::

           0    0    0    0 0.5037555283 0.0000000000
           0    0    0    1 0.4610257618 0.0003115503
           0    1    0    0 0.4610257618 -0.0003115503
           0    1    0    1 0.4962444717 0.0000000000
        ...

-  | zvo\_eigen.dat:
     収束したハミルトニアンの固有値が低エネルギー順に出力されます。

   ::

        1  -2.9425069199
        2  -2.9425069198
        3  -1.5005359205 
        ...

-  zvo\_gap.dat:
   全電子数を :math:`N_{\rm tot}` とした場合に、 :math:`\Delta E= E(N_{\rm tot}+1)-E(N_{\rm tot})` が出力されます。

   ::

         5.2208232631

-  zvo\_orbital\_opt.dat:
   スレータ行列式から生成した :math:`f_{ij}` 。 ``InOrbital``, ``InOrbitalAntiParallel``,
   ``InOrbitalParallel``, ``InOrbitalAntiGeneral`` ファイルと同じ形式のファイルが出力されます。
   :math:`f_{ij}` が ``Orbital``, ``OrbitalAntiParallel``,
   ``OrbitalParallel``, ``OrbitalAntiGeneral`` ファイルを参照し計算され、同種のパラメータについては平均化した値が採用されます。
