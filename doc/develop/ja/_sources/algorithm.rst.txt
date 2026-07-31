.. include:: ../../bib/ref.txt

アルゴリズム
============

変分モンテカルロ法
------------------

変分モンテカルロ法では, 試行波動関数を用意して,
その試行波動関数が含むパラメータを
変分原理に従って最適化することで量子多体系の基底状態(または低励起エネルギー状態)の波動関数を
近似的に求めます。 試行波動関数に対する物理量の期待値を計算する部分で,
マルコフ連鎖モンテカルロ法を利用し,
効率よく重み付きサンプリングを行います。

本パッケージでは,
サンプリングに用いる完全系として電子の実空間配置 :math:`| x\rangle` をとっています:

.. math::

   | x\rangle = \prod_{n=1}^{N_e/2} c_{r_{n\uparrow}}^{\dagger} \prod_{n=1}^{N_e/2}
   c_{r_{n\downarrow}}^{\dagger} |0 \rangle

ここで,
:math:`r_{n\sigma}` は :math:`n` 番目の電子(スピン :math:`\sigma`)の位置,
:math:`c_{r_{n\sigma}}^{\dagger}` はその位置での
電子(スピン :math:`\sigma`)の生成演算子を表します。この基底を用いると,
演算子 :math:`A` の期待値は

.. math::

   \langle A \rangle =\frac{\langle \psi| A| \psi \rangle}{\langle \psi | \psi \rangle} 
   =\sum_x \frac{\langle \psi| A | x\rangle \langle x| \psi \rangle}{\langle \psi |\psi \rangle}

となるため, マルコフ連鎖の重みを

.. math:: \rho(x)=\frac{|\langle x| \psi \rangle|^2}{\langle \psi | \psi \rangle} \ge 0, \quad \sum_{x} \rho(x)=1

 と定義して,

.. math::

   \langle A \rangle =\sum_x \rho(x) \frac{\langle \psi| A | x\rangle }{\langle \psi |x \rangle}


と書き直した後、:math:`x` に関する和をマルコフ連鎖モンテカルロ法により
評価しています。Local Green’s function
:math:`G_{ij\sigma\sigma'}(x)` は

.. math::

   G_{ij\sigma\sigma'}(x)=\frac{\langle \psi | c_{i\sigma}^{\dagger} c_{j\sigma'}
   | \psi \rangle}{\langle \psi | x \rangle}

と定義されますが,
これも演算子 :math:`A` を :math:`c_{i\sigma}^{\dagger} c_{j\sigma'}` ととることで,
同じ方法により重み付きサンプリングを行うことができます。 なお,
サンプリングに使用する乱数生成については,
メルセンヌツイスター法を使用しています[Mutsuo2008_ ]。

.. _BogoliubovRep:

Bogoliubov表現
--------------

スピン系の計算において一体項( ``transfer``),
``InterAll`` 形式での相互作用,
相関関数のインデックスの指定にはBogoliubov表現が使われています。
一般に、スピンの演算子は次のようにフェルミオンの生成 :math:`\cdot` 消滅演算子 :math:`c_{i \sigma}^\dagger`,
:math:`c_{i \sigma}` によって書き換えることができます:

.. math::

   \begin{aligned}
   S_{i z} &= \sum_{\sigma = -S}^{S} \sigma c_{i \sigma}^\dagger c_{i \sigma}
   \\
   S_{i}^+ &= \sum_{\sigma = -S}^{S-1} 
   \sqrt{S(S+1) - \sigma(\sigma+1)} 
   c_{i \sigma+1}^\dagger c_{i \sigma}
   \\
   S_{i}^- &= \sum_{\sigma = -S}^{S-1} 
   \sqrt{S(S+1) - \sigma(\sigma+1)} 
   c_{i \sigma}^\dagger c_{i \sigma+1}
   \end{aligned}

本パッケージでは、 :math:`S=1/2` のスピン系のみ取り扱っており、上記の式で
:math:`S=1/2` と置いたものを用いています。

.. _PuffAndSlater:
      
パフィアン-スレーター行列式の性質
---------------------------------

この節では,
パフィアン-スレーター行列式のもつ性質について簡単にまとめます。
:ref:`次節 <PfaffianAP>` と :ref:`次々節 <PfaffianP>` でパフィアン-スレーター行列式と単一スレーター行列式の間の関係を導出し、
:ref:`最後 <PfaffianSingular>` に :math:`f_{ij}` の特異値分解の意味について説明します。

.. _PfaffianAP:

:math:`f_{ij}` と :math:`\Phi_{in\sigma}` の関係 (スピン反平行の場合)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

多変数変分モンテカルロ法で試行波動関数の一体部分として用いられるパフィアン-スレーター行列式は

.. math::

   |\phi_{\rm Pf}\rangle=\Big(\sum_{i,j=1}^{N_{s}}f_{ij}
   c_{i\uparrow}^{\dagger}c_{j\downarrow}^{\dagger}\Big)^{N_{\rm e}/2}|0\rangle,

のように定義されます。ここで, :math:`N_{s}` はサイト数,
:math:`N_{e}` は全電子数, :math:`f_{ij}` は変分パラメータです。
簡単化のため, 以降 :math:`f_{ij}` は実数と仮定します。また,
単一スレーター行列式として

.. math::

   \begin{aligned}
   |\phi_{\rm SL}\rangle&=\Big(\prod_{n=1}^{N_{e}/2}\psi_{n\uparrow}^{\dagger}\Big)
   \Big(\prod_{m=1}^{N_{e}/2}\psi_{m\downarrow}^{\dagger}\Big)|0\rangle, \\
   \psi_{n\sigma}^{\dagger}&=\sum_{i=1}^{N_{s}}\Phi_{in\sigma}c^{\dagger}_{i\sigma}.
   \end{aligned}

を定義します。ただし, :math:`\Phi` は正規直交基底であり,
クロネッカーのデルタ :math:`\delta_{nm}` を用い

.. math::

   \sum_{i=1}^{N_{s}}\Phi_{in\sigma}\Phi_{im\sigma}=\delta_{nm},

で表されます。この直交性の関係から, 以下の関係式

.. math::

   \begin{aligned}
   {[\psi^{\dagger}_{n\sigma}, \psi_{m\sigma}]}_{+}&=\delta_{nm},\\
   G_{ij\sigma}=\langle c_{i\sigma}^{\dagger}c_{j\sigma}\rangle 
   &=\frac{\langle \phi_{\rm SL}| c_{i\sigma}^{\dagger}c_{j\sigma} | \phi_{\rm SL}\rangle}{\langle \phi_{\rm SL}|\phi_{\rm SL}\rangle } \\
   &=\sum_{n} \Phi_{in\sigma} \Phi_{jn\sigma}.
   \end{aligned}

が導かれます。

次に, :math:`|\phi_{\rm SL}\rangle` を変形し,
:math:`f_{ij}` と :math:`\Phi_{in\sigma}` の間に成り立つ関係式を示します。
:math:`\psi^{\dagger}_{n\sigma}` の交換関係を用いると,
:math:`|\phi_{\rm SL}\rangle` は

.. math::

   \begin{aligned}
   |\phi_{\rm SL}\rangle \propto \prod_{n=1}^{N_{e}/2}
   \Big(\psi_{n\uparrow}^{\dagger}\psi_{\mu(n)\downarrow}^{\dagger}\Big)|0\rangle,
   \end{aligned}

と書き換えられます。ここで,
:math:`\mu(n)` は :math:`n= 1, 2, \cdots, N_{e}/2` の置換を表します。
議論を簡単にするため,
同一のペア :math:`n=\mu(n)` となる場合を考えましょう。 このとき,
:math:`K_{n}^{\dagger}=\psi_{n\uparrow}^{\dagger}\psi_{n\downarrow}^{\dagger}` として,
:math:`K_{n}^{\dagger}K_{m}^{\dagger}=K_{m}^{\dagger}K_{n}^{\dagger}` の関係を用いることで,

.. math::

   \begin{aligned}
   |\phi_{\rm SL}\rangle &\propto \prod_{n=1}^{N_{e}/2}\Big(\psi_{n\uparrow}^{\dagger}\psi_{n\downarrow}^{\dagger}\Big)|0\rangle
   =\prod_{n=1}^{N_{e}/2} K_{n}^{\dagger}|0\rangle \\
   &\propto\Big(\sum_{n=1}^{\frac{N_{e}}{2}}K_{n}^{\dagger}\Big)^{\frac{N_{e}}{2}} |0\rangle
   =\Big(\sum_{i,j=1}^{N_{s}}\Big[\sum_{n=1}^{\frac{N_{e}}{2}}\Phi_{in\uparrow}\Phi_{jn\downarrow}\Big]
   c_{i\uparrow}^{\dagger}c_{j\downarrow}^{\dagger}\Big)^{N_e/2}|0\rangle,
   \end{aligned}

の関係が得られます。これより :math:`f_{ij}` は単一スレーター行列式の係数により

.. math::

   \begin{aligned}
   f_{ij}=\sum_{n=1}^{\frac{N_{e}}{2}}\Phi_{in\uparrow}\Phi_{jn\downarrow}.
   \end{aligned}

として表されることが分かります。なお,
この形式は単一スレーター行列式で与えられる :math:`f_{ij}` の表式の一つであり,
実際にはペアを組む自由度(どの :math:`\mu(n)` を選ぶか)およびゲージの自由度
(すなわち :math:`\Phi_{in\sigma}` の符号の自由度)に依存します。
この自由度の多さが :math:`f_{ij}` の冗長性につながっています。

.. _PfaffianP:

:math:`F_{IJ}` と :math:`\Phi_{In}` の関係 (スピン平行も含めた場合)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

前節で考察したパフィアン-スレーター波動関数と単一スレーター波動関数の間の関係は、
同種スピンのペアリングも考えた場合に拡張することができます。
パフィアン-スレーター波動関数とスレーター波動関数をそれぞれ

.. math::

   \begin{aligned}
   |\phi_{\rm Pf}\rangle&=\Big(\sum_{I,J=1}^{2N_{s}}F_{IJ}c_{I}^{\dagger}c_{J}^{\dagger}\Big)^{N_{\rm e}/2}|0\rangle, \\
   |\phi_{\rm SL}\rangle&=\Big(\prod_{n=1}^{N_{e}}\psi_{n}^{\dagger}\Big)|0\rangle,~~\psi_{n}^{\dagger}=\sum_{I=1}^{2N_{s}}\Phi_{In}c^{\dagger}_{I}.\end{aligned}

と定義します。ここで :math:`I,J` はスピン自由度も含めたサイトのインデックスです。
スピン反平行の場合とほぼ同様の議論を用いることで,

.. math::

   \begin{aligned}
   F_{IJ}=\sum_{n=1}^{\frac{N_{e}}{2}}\Big(\Phi_{I,2n-1}\Phi_{J,2n}-\Phi_{J,2n-1}\Phi_{I,2n}\Big).
   \end{aligned}

の関係を示すことができます。これはスピン反平行のペアリングにもそのまま適用できるので,
mVMC ver1.0以降ではこの表式を使用しています。

.. _PfaffianSingular:

:math:`f_{ij}` の特異値分解
~~~~~~~~~~~~~~~~~~~~~~~~~~~

行列 :math:`F`, :math:`\Phi_{\uparrow}`, :math:`\Phi_{\downarrow}`,
:math:`\Sigma` を

.. math::

   \begin{aligned}
   &(F)_{ij}=f_{ij},~~~ 
   (\Phi_{\uparrow})_{in}=\Phi_{in\uparrow},~~~ 
   (\Phi_{\downarrow})_{in}=\Phi_{in\downarrow}, \\
   &\Sigma={\rm diag}[\underbrace{1,\cdots,1}_{N_e/2},0,0,0],\end{aligned}

として定義します。前節のように :math:`f_{ij}` （すなわち :math:`F`）が単一スレーター行列と関係づけられて
いるとき、 :math:`F` の特異値分解は

.. math::

   \begin{aligned}
   F=\Phi_{\uparrow}\Sigma\Phi_{\downarrow}^{t}.\end{aligned}

となることを示すことができます。
この結果は、一般に :math:`F` を特異値分解したとき、非ゼロの特異値が :math:`N_{e}/2` 個存在し,
かつ全ての :math:`F` の非ゼロの特異値が :math:`1` であった場合,
:math:`f_{ij}` が単一スレーター波動関数を
記述すること(つまり平均場近似解として記述できること)を表しています。
言い換えると, 特異値の非ゼロ成分の数とその値が,
シングルスレータ行列式からパフィアンスレーター行列式がどのようにしてずれるのか,
という点について定量的な基準を与えることを示しています。

Power-Lanczos法
---------------

Power-Lanczos法では、最適化済みのVMC波動関数 :math:`|\psi\rangle` に
Hamiltonianの多項式を作用させた

.. math::

   |\phi_p\rangle
   =P_p(\hat H)|\psi\rangle
   =\sum_{k=0}^{p}c_k\hat H^k|\psi\rangle

を改良された試行波動関数とします。これは
:math:`\{|\psi\rangle,\hat H|\psi\rangle,\ldots,\hat H^p|\psi\rangle\}`
が張るKrylov部分空間で、エネルギーを最小にする係数 :math:`c_k` を求める方法です。
mVMCで用いる1st stepと2nd stepの試行波動関数は、それぞれ

.. math::

   \begin{aligned}
   |\phi_1\rangle &= (1+\alpha\hat H)|\psi\rangle,\\
   |\phi_2\rangle &= (1+\alpha\hat H+\beta\hat H^2)|\psi\rangle
   \end{aligned}

です。1st stepは ``NVMCCalMode=1``、``NLanczosMode>=1``、
``NLanczosStep=1`` で実行されます。2nd stepには
``NVMCCalMode=1``、``NLanczosMode=1``、``NLanczosStep=2`` が必要です。

以下では、まず1st stepの :math:`\alpha` の決定方法を説明し、
次に2nd stepの :math:`\alpha,\beta` を決める変分原理と、
その問題を数値的に安定して解くための実装を分けて説明します。
最後に、1st stepを適用した後の物理量の計算方法を説明します。

1st power-Lanczos stepと :math:`\alpha` の決定
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

最初に, 変分モンテカルロ法のサンプリングに関して簡単に説明します。
ここで :math:`x` は基底VMC波動関数 :math:`|\psi\rangle` からサンプリングする
配置です。物理量 :math:`\hat{A}` の期待値は

.. math::

   \begin{aligned}
   &\langle \hat{A}\rangle_\psi
   = \frac{\langle \psi|\hat{A}|\psi\rangle}{\langle \psi|\psi\rangle}
   = \sum_x \rho_\psi(x)F(x,\hat{A}),\\
   &\rho_\psi(x)=\frac{|\langle x|\psi\rangle|^2}{\langle \psi|\psi\rangle},
   \qquad
   F(x,\hat{A})=\frac{\langle x|\hat{A}|\psi\rangle}
                      {\langle x|\psi\rangle}
   \end{aligned}

と書けます。:math:`F(x,\hat A)` は配置 :math:`x` における
:math:`\hat A` のlocal estimatorです。演算子の積
:math:`\hat{A}\hat{B}` の期待値には、以下の二通りの表現があります。

.. math::

   \begin{aligned}
   &\langle \hat{A}\hat{B}\rangle_\psi
   = \sum_x \rho_\psi(x)F(x,\hat{A}\hat{B}),\\
   &\langle \hat{A}\hat{B}\rangle_\psi
   = \sum_x \rho_\psi(x)F^\dagger(x,\hat{A})F(x,\hat{B})
   \end{aligned}

前者は積 :math:`\hat{A}\hat{B}` 全体を1つのlocal estimatorとして
評価する方法、後者は :math:`\hat{A}` と :math:`\hat{B}` の間に完全系
:math:`\sum_{x} |x\rangle \langle x|` を挿入し、:math:`\hat{A}` を
ブラ側、:math:`\hat{B}` をケット側に分割して評価する方法です。
後者の式では :math:`\hat A` がHermite演算子であることを仮定しています。
一般の :math:`\hat A` に対しては
:math:`F^\dagger(x,\hat A)` を :math:`F^\dagger(x,\hat A^\dagger)` で
置き換える必要があります。
両者はサンプル数無限大の極限では同じ期待値を与えますが、
有限のサンプル数では一致せず、一般に後者の方が数値的に安定です。
以下では表記を簡単にするため :math:`\rho_\psi(x)` を
:math:`\rho(x)` と略記します。

例えば、エネルギーの分散
:math:`\sigma^2=\langle (\hat{H}-\langle \hat{H}\rangle)^2\rangle` を
それぞれの方法で計算すると以下のようになります。

.. math::

   \begin{aligned}
   \sigma^2 &=\sum_{x} \rho(x) F(x,  (\hat{H}-\langle \hat{H}\rangle)^2) = \sum_{x} \rho(x) F(x,  \hat{H}^2) - \left[ \sum_{x} \rho(x) F(x,  \hat{H})\right]^2 ,\\
   \sigma^2 &=\sum_{x} \rho(x) F^{\dagger}(x,  \hat{H}-\langle \hat{H}\rangle)F(x,  \hat{H}-\langle \hat{H}\rangle) \nonumber \\
   &= \sum_{x} \rho(x) F^{\dagger}(x,  \hat{H}) F(x, \hat{H})- \left[ \sum_{x} \rho(x) F(x,  \hat{H})\right]^2
   \end{aligned}

後者の方法では有限サンプリングでも分散が非負となることが
保証されているのに対して、前者の方法ではそれが保証されません。

次に、1st stepの波動関数
:math:`|\phi_1\rangle =(1+\alpha \hat{H}) |\psi \rangle` の
エネルギー期待値を考えます。:math:`\alpha` について展開すると
(分子・分母を :math:`\langle \psi | \psi \rangle` で割っています)、

.. math::

   \begin{aligned}
   E_{LS}(\alpha) =\frac{\langle \phi_1| \hat{H} |\phi_1\rangle}{\langle \phi_1|\phi_1\rangle}
   =\frac{\langle \hat{H} \rangle_\psi + 2\alpha \langle \hat{H}^2 \rangle_\psi + \alpha^2 \langle \hat{H}^3 \rangle_\psi}
   {1 + 2\alpha \langle \hat{H} \rangle_\psi + \alpha^2 \langle \hat{H}^2 \rangle_\psi}
   \end{aligned}

となり、:math:`\hat{H}` の3次までのモーメントが必要になります。
これらのモーメントを上述の「分割」によって推定するため、
以下の記法を導入します:

.. math::

   \begin{aligned}
   h_{n(ij)}=\sum_{x} \rho(x) F^{\dagger}(x,  \hat{H}^i) F(x, \hat{H}^j), ~~~~ i+j=n.
   \end{aligned}

すなわち、:math:`h_{n(ij)}` は :math:`n` 次のモーメント
:math:`\langle \hat{H}^n \rangle_\psi` を、ブラ側に :math:`\hat{H}^i`、
ケット側に :math:`\hat{H}^j` を割り振って推定した量です
(:math:`F(x, \hat{H}^0)=1` とし、:math:`h_{1(10)}` は :math:`h_1` と
略記します)。計算に必要となるのは以下の量です:

.. math::

   \begin{aligned}
   &h_1 =\sum_{x} \rho(x) F^{\dagger}(x,  \hat{H}),\\
   &h_{2(11)}=\sum_{x} \rho(x) F^{\dagger}(x,  \hat{H}) F(x, \hat{H}),\\
   &h_{2(20)}=\sum_{x} \rho(x) F^{\dagger}(x,  \hat{H}^2),\\
   &h_{3(12)}=\sum_{x} \rho(x) F^{\dagger}(x,  \hat{H})F(x,  \hat{H}^2),\\
   &h_{4(22)}=\sum_{x} \rho(x) F^{\dagger}(x,  \hat{H}^2)F(x,  \hat{H}^2).\end{aligned}

:math:`h_{2(11)}` と :math:`h_{2(20)}` はどちらも
:math:`\langle \hat{H}^2 \rangle_\psi` の推定量ですが、
分割の仕方が異なるため有限のサンプル数では一致しません。
また、次数をブラ側とケット側へできるだけ均等に割り振ることで、
local estimatorとしては :math:`F(x, \hat{H}^2)` までしか現れず、
:math:`\hat{H}^3` や :math:`\hat{H}^4` のlocal estimatorを
直接評価する必要がなくなっています。

これらの推定量を用いると、エネルギー期待値は

.. math::

   \begin{aligned}
   E_{LS}(\alpha) =\frac{h_1 + \alpha(h_{2(20)} + h_{2(11)}) + \alpha^2 h_{3(12)}}{1 + 2\alpha h_1 + \alpha^2 h_{2(11)}}
   \end{aligned}

と評価されます。分子の :math:`\alpha` の1次の項に現れる2つの
:math:`\langle \hat{H}^2 \rangle_\psi` は、それぞれブラ側・ケット側の
:math:`(1+\alpha \hat{H})` に由来する別々の項であるため、
一方を :math:`h_{2(11)}`、もう一方を :math:`h_{2(20)}` で推定しています。

停留条件 :math:`\partial E_{LS}(\alpha)/\partial \alpha=0` は、
:math:`\alpha` に関する以下の二次方程式に帰着します:

.. math::

   \begin{aligned}
   \left[2 h_1 h_{3(12)} - h_{2(11)}\left(h_{2(11)}+h_{2(20)}\right)\right] \alpha^2
   + 2\left[h_{3(12)} - h_1 h_{2(11)}\right] \alpha
   + \left(h_{2(11)}+h_{2(20)}\right) - 2h_1^2 = 0.
   \end{aligned}

この二次方程式の2つの実数解 :math:`\alpha_{\pm}` のそれぞれについて
:math:`E_{LS}(\alpha_{\pm})` を評価し、エネルギーが低くなる方を
最適な :math:`\alpha` として採用します
(判別式が負となり実数解が存在しない場合には、
Lanczos計算はエラーとして終了します)。
また、最適化された :math:`|\phi_1\rangle` に対するエネルギーの分散も

.. math::

   \begin{aligned}
   \sigma^2_{LS}(\alpha) = \frac{\langle \phi_1| \hat{H}^2 |\phi_1\rangle}{\langle \phi_1|\phi_1\rangle} - E_{LS}(\alpha)^2,~~~~
   \frac{\langle \phi_1| \hat{H}^2 |\phi_1\rangle}{\langle \phi_1|\phi_1\rangle}
   = \frac{h_{2(11)} + 2\alpha h_{3(12)} + \alpha^2 h_{4(22)}}{1 + 2\alpha h_1 + \alpha^2 h_{2(11)}}
   \end{aligned}

から同様に計算できます。

2nd power-Lanczos step
~~~~~~~~~~~~~~~~~~~~~~

変分原理による :math:`\alpha,\beta` の決定
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

``NLanczosStep=2`` では、
:math:`|\phi_2\rangle=(1+\alpha\hat H+\beta\hat H^2)|\psi\rangle`
の変分エネルギー

.. math::

   E_2(\alpha,\beta)
   =\frac{\langle\phi_2|\hat H|\phi_2\rangle}
          {\langle\phi_2|\phi_2\rangle}

を最小にする :math:`\alpha,\beta` を求めます。これは1st stepで
:math:`E_{LS}(\alpha)` の停留条件から :math:`\alpha` を決めることの
直接の拡張です。

この最小化は、Krylov基底

.. math::

   |v_a\rangle=\hat H^a|\psi\rangle,\qquad a=0,1,2

における一般化固有値問題として表せます。overlap行列とHamiltonian行列を

.. math::

   S_{ab}=\frac{\langle v_a|v_b\rangle}{\langle\psi|\psi\rangle},
   \qquad
   {\cal H}_{ab}
   =\frac{\langle v_a|\hat H|v_b\rangle}{\langle\psi|\psi\rangle}

と定義します。:math:`|v_a\rangle` の定義から、これらの行列要素は
基底波動関数に対するHamiltonianのモーメントそのものです:

.. math::

   S_{ab}=\langle \hat H^{a+b}\rangle_\psi,
   \qquad
   {\cal H}_{ab}=\langle \hat H^{a+b+1}\rangle_\psi.

したがって2nd stepでは :math:`\langle \hat H^5\rangle_\psi` までの
モーメント(後述する分散の評価まで含めると
:math:`\langle \hat H^6\rangle_\psi` まで)が必要になります。
:math:`\boldsymbol{c}=(c_0,c_1,c_2)^{\mathsf T}` とすると、

.. math::

   |\phi_2\rangle=\sum_{a=0}^{2}c_a|v_a\rangle,\qquad
   E_2(\boldsymbol{c})
   =\frac{\boldsymbol{c}^\dagger{\cal H}\boldsymbol{c}}
          {\boldsymbol{c}^\dagger S\boldsymbol{c}}.

したがって停留条件は

.. math::

   {\cal H}\boldsymbol{c}=E_2 S\boldsymbol{c}

となり、最低固有値に対応する固有ベクトルが最適な係数を与えます。
固有ベクトルの全体のスケールは任意なので、
:math:`\alpha=c_1/c_0`、:math:`\beta=c_2/c_0` とすれば、
冒頭の :math:`(1+\alpha\hat H+\beta\hat H^2)|\psi\rangle`
という表式に戻ります。以上が2nd power-Lanczos stepの原理です。
なお、:math:`c_2=0` に固定して :math:`c_1/c_0` のみを最適化すると、
この一般化固有値問題は1st stepの二次方程式と等価になります。
つまり2nd stepは、1st stepと同じ変分原理を、Krylov部分空間の次元を
2から3へ広げて適用したものです。

数値的安定性のための実装
^^^^^^^^^^^^^^^^^^^^^^^^^^

ここからは上の変分原理を変えるものではなく、同じ行列要素をMonte Carlo
sampleから安定に評価し、一般化固有値問題を安定に解くための手順です。

まず、上で定義したlocal estimator :math:`F(x,\hat A)` へ
:math:`\hat A=\hat H^a` を代入したlocal power

.. math::

   F_a(x)\equiv F(x,\hat H^a)
   =\frac{\langle x|\hat H^a|\psi\rangle}
          {\langle x|\psi\rangle},
   \qquad a=0,1,2,3

を各配置 :math:`x` で計算し、Monte Carlo sampleごとに保持します。
:math:`F_0(x)=1`、:math:`F_1(x)` は通常のlocal energyであり、
:math:`F_2(x)` と :math:`F_3(x)` はそれぞれ :math:`\hat H^2` と
:math:`\hat H^3` のlocal estimatorです。これらは試行波動関数そのものではなく、
:math:`S`、:math:`{\cal H}` および分散計算用の
Hamiltonian二乗行列をMonte Carlo平均から構成するための量です。

sample loopの終了後に、基底VMC波動関数のglobalなweighted平均エネルギー
:math:`E_{\rm shift}=\langle\hat H\rangle_\psi` を求め、
:math:`X=\hat H-E_{\rm shift}` とします。この置換は

.. math::

   {\rm span}\{|\psi\rangle,\hat H|\psi\rangle,\hat H^2|\psi\rangle\}
   =
   {\rm span}\{|\psi\rangle,X|\psi\rangle,X^2|\psi\rangle\}

を満たすため、変分空間と最適解を変えません。これは数値的安定性のための
基底変換です。momentの積を蓄積する前にlocal powerを

.. math::

   \widetilde F_a(x)
   =\sum_{j=0}^{a}{a\choose j}(-E_{\rm shift})^{a-j}F_j(x)

へ変換します。:math:`\widetilde F_a` は :math:`X^a` のlocal estimator
です。mVMCはcentered moment行列

.. math::

   \widetilde M_{ab}
   =\sum_x \rho(x)\widetilde F_a(x)^\dagger\widetilde F_b(x)

を構成します。実装は :math:`a,b=0,\ldots,3` の全組合せを保持します。
各 :math:`\widetilde M_{ab}` は、1st stepの :math:`h_{n(ij)}` と同様に、
:math:`n=a+b` 次のモーメント :math:`\langle X^{a+b}\rangle_\psi` を
ブラ側 :math:`\widetilde F_a^\dagger` とケット側
:math:`\widetilde F_b` に分割した推定量です。各側の次数が最大3なので、
:math:`\widetilde M_{33}` として :math:`X^6` のモーメントまで表しながら、
local powerとして必要なのは :math:`F_3` までです。
このaccumulation前の中心化により、extensiveなエネルギーが大きい場合に
:math:`\langle H^6\rangle` などのraw momentを形成して生じる桁落ちを避けます。
必要な追加memoryはMonte Carlo sampleあたり4個のlocal powerに比例します。
ただし、中心化が縮小するのは主に平均エネルギーによるoffsetです。
:math:`|\langle H\rangle|` がエネルギー分布の広がりに比べて小さい場合は
高次local powerの動的範囲が十分に縮まらないため、中心化だけで任意の系の
数値精度が保証されるわけではありません。

同じsampleから構成する :math:`\widetilde M_{ab}` と
:math:`\widetilde M_{ba}^\dagger` は各sampleの寄与ごとに等しいため、
有限sampleであること自体は両者の差を生みません。ただし、浮動小数点演算や
並列reductionの丸め誤差によって小さな非Hermite成分が生じ得るため、
mVMCはmoment行列をHermite化します。一方、
:math:`a+b=c+d` を満たす異なる分割
:math:`\widetilde M_{ab}` と :math:`\widetilde M_{cd}` の差には
有限sampleの統計誤差が含まれ、``hankel_residual`` として出力されます。
その後、centered Krylov基底のoverlap行列、
Hamiltonian行列、Hamiltonian二乗行列を

.. math::

   \widetilde S_{ab}=\widetilde M_{ab},\qquad
   \widetilde{\cal H}_{ab}
   =\frac{\widetilde M_{a,b+1}+\widetilde M_{a+1,b}}{2},\qquad
   \widetilde G_{ab}=\widetilde M_{a+1,b+1}

と構成します。ここで :math:`a,b=0,1,2` です。
:math:`\widetilde{\cal H}` の二つのmomentを平均することも、
有限sampleでのHermite性を保つための処理です。

さらに各基底ベクトルのノルムの差を抑えるため、

.. math::

   D_{aa}=\frac{1}{\sqrt{\widetilde S_{aa}}}

による対角scale変換を行い、

.. math::

   (D\widetilde{\cal H}D)\boldsymbol{u}
   =E_X(D\widetilde S D)\boldsymbol{u}

を解きます。overlap行列の分解に失敗した場合だけ、
:math:`D\widetilde S D` の対角へ小さな正則化項を加えて再試行します。
中心化、Hermite化、対角scale変換およびこのフォールバックの正則化が、
原理とは別に導入されている数値安定化です。

最低固有値から元のHamiltonianのエネルギー
:math:`E=E_X+E_{\rm shift}` を得ます。centered basisの係数を
:math:`\boldsymbol{d}=(d_0,d_1,d_2)^{\mathsf T}` とすると、

.. math::

   \begin{aligned}
   c_0 &= d_0-E_{\rm shift}d_1+E_{\rm shift}^2d_2,\\
   c_1 &= d_1-2E_{\rm shift}d_2,\\
   c_2 &= d_2
   \end{aligned}

によってunshifted basisへ戻します。mVMCは :math:`c_0,c_1,c_2` と、
上で定義した :math:`\alpha,\beta` をそれぞれ
:math:`\alpha_1,\alpha_2` として出力します。分散は
:math:`\boldsymbol{d}^\dagger\widetilde S\boldsymbol{d}=1` の規格化のもとで

.. math::

   \sigma_E^2
   =\boldsymbol{d}^\dagger\widetilde G\boldsymbol{d}
    -(\boldsymbol{d}^\dagger\widetilde{\cal H}\boldsymbol{d})^2

としてcentered basisで評価します。

計算量
^^^^^^

2nd stepで新たに必要になるのは、sampleごとのlocal power
:math:`F_2(x)` と :math:`F_3(x)` の評価であり、このうち
:math:`F_3(x)=\langle x|\hat H^3|\psi\rangle/\langle x|\psi\rangle`
が計算量を支配します。

現在の2nd step実装には2種類のHamiltonian classがあります。
electronic classはspinを保存する ``Transfer`` 項とnumber-operator型の
対角相互作用を含みます。pure spin-1/2 classは対角相互作用と
``Exchange`` 項を含み、全siteの1電子占有を保存し、現時点では
``NQPFull=1`` に限定されます。

electronic classの :math:`F_3` 評価では、非対角な ``Transfer`` 項に対する
outer/innerの二重loopを回します。2つの ``Transfer`` 項の積は二体演算子

.. math::

   CACA \equiv
   c_{i\sigma}^\dagger c_{j\sigma} c_{k\tau}^\dagger c_{l\tau}

の形にまとまるため、二重loopで生成される
:math:`O(N_{\mathrm{Transfer}}^2)` 個の ``CACA`` のそれぞれについて、
残りの :math:`\hat H` を含む行列要素
:math:`\langle \psi|\hat H\,CACA|x\rangle/\langle \psi|x\rangle`
を評価することになります。

量子射影がidentity（``NQPFull=1``）の場合は、``CACA`` を
サンプリングされた配置 :math:`x` に作用させて一時的な配置
:math:`x'` を作り、:math:`x'` における :math:`\hat H` の
local estimatorに重なり比
:math:`\langle \psi|CACA|x\rangle/\langle \psi|x\rangle` を掛ける、
という分解で評価できます。二重loopが生成する ``CACA`` の個数は
:math:`O(N_{\mathrm{Transfer}}^2)` ですが、各 ``CACA`` に対する
Hamiltonianのlocal estimatorがさらに ``Transfer`` loopを含むため、
この経路の最悪計算量も
:math:`O(N_{\mathrm{Transfer}}^3)` です。ただし、通常のlocal
estimatorとrank-two updateを利用でき、``NQPFull`` に比例する
``GreenFuncN`` の評価を避けられるため、非自明な射影の経路より
計算量の係数は小さくなります。

非自明な量子射影（``NQPFull>1``）では、この分解は複数の量子射影成分の
和を正しく保てません。そのため2nd power-Lanczosでは、
:math:`\hat H=\sum_\mu h_\mu\hat O_\mu` の各項について
:math:`\hat O_\mu\,CACA` を一つの多体演算子として扱い、
``GreenFuncN`` で対応する多体Green関数を元の配置 :math:`x` から
量子射影成分を含めて評価し、係数 :math:`h_\mu` 付きで加算します。
この直接評価では、各 :math:`\hat O_\mu\,CACA` の多体Green関数を
求める際に内部へもう1段の ``Transfer`` loopが加わるため、
この経路の支配的なoperator数は :math:`O(N_{\mathrm{Transfer}}^3)`
となり、``GreenFuncN`` の計算量は ``NQPFull`` にも比例します。
したがって実時間は波動関数と計算環境だけでなく、
射影設定にも強く依存します。

pure-spin classでは各exchange bondを2つの向き付き演算子
:math:`X_a` の和として表します。直接縮約により

.. math::

   \begin{aligned}
   F_3(x)={}&V_0^3
   +\sum_a j_a(V_0^2+V_0V_a+V_a^2)G_a\\
   &+\sum_{a,b}j_bj_a(V_0+V_a+V_{ba})G_{ba}
   +\sum_{a,b,c}j_cj_bj_aG_{cba}
   \end{aligned}

を評価します。ここで :math:`V_a` と :math:`V_{ba}` は各exchange作用後の
対角energy、:math:`G_a`、:math:`G_{ba}`、:math:`G_{cba}` は元のsample配置から
それぞれ2、4、6次の ``GreenFuncN`` で求める行列要素です。
途中配置の波動関数振幅で割らないため、その振幅が0でも有効です。
向き付きexchangeが :math:`A` 個activeなら、深さ3のcall数は
:math:`A^3` 以下です。

1st power-Lanczos stepでの物理量の計算
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``NLanczosMode=2`` の場合には、最適化されたパラメータ :math:`\alpha` を
用いて、1st step波動関数 :math:`|\phi_1\rangle` に対する物理量
(一体・二体グリーン関数)を計算します
(2nd stepでのGreen関数計算は未対応です)。
演算子 :math:`\hat{A}` の期待値の分子を :math:`\alpha` について
展開すると

.. math::

   \begin{aligned}
   \langle \phi_1| \hat{A} |\phi_1\rangle = \langle \psi| \hat{A} |\psi\rangle
   + \alpha \langle \psi| \hat{H}\hat{A} |\psi\rangle
   + \alpha \langle \psi| \hat{A}\hat{H} |\psi\rangle
   + \alpha^2 \langle \psi| \hat{H}\hat{A}\hat{H} |\psi\rangle
   \end{aligned}

となります。エネルギーの場合と同様に、各項をブラ側の :math:`\hat{H}^i` と
ケット側の :math:`\hat{A}\hat{H}^j` に分割して推定します:

.. math::

   \begin{aligned}
   &A_0 =\sum_{x} \rho(x) F(x,  \hat{A}),\\
   &A_{1(10)}=\sum_{x} \rho(x) F^{\dagger}(x,  \hat{H}) F(x, \hat{A}),\\
   &A_{1(01)}=\sum_{x} \rho(x) F(x, \hat{A}\hat{H}),\\
   &A_{2(11)}=\sum_{x} \rho(x) F^{\dagger}(x,  \hat{H})F(x,  \hat{A}\hat{H}).\end{aligned}

ここで添字 :math:`n(ij)` は、エネルギーの場合と同様に、
:math:`\hat{H}` の合計次数 :math:`n=i+j` と、
ブラ側(:math:`\hat{H}^i`)・ケット側(:math:`\hat{A}\hat{H}^j`)への
割り振りを表します。これらの推定量を用いると、期待値は

.. math::

   \begin{aligned}
   A_{LS}(\alpha) =\frac{\langle \phi_1| \hat{A} |\phi_1\rangle}{\langle \phi_1|\phi_1\rangle}=\frac{A_0 + \alpha(A_{1(10)} + A_{1(01)}) + \alpha^2 A_{2(11)}}{1 + 2\alpha h_1 + \alpha^2 h_{2(11)}}
   \end{aligned}

と計算されます。分母はエネルギーの場合と同じ規格化因子
:math:`\langle \phi_1 | \phi_1 \rangle` です。プログラムでは、
この表式に基づき一体グリーン関数および二体グリーン関数の計算を
行っています。
