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

   | x\rangle = \prod_{n=1}^{N/2} c_{r_{n\uparrow}}^{\dagger} \prod_{n=1}^{N/2}
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
一般に、スピンの演算子は次のようにフェルミオンの生成 :math:`\cdot` 消滅演算子 :math:`c_{i \sigma}`,
:math:`c_{i \sigma}^\dagger` によって書き換えることができます:

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
   [\psi^{\dagger}_{n\sigma},\psi_{m\sigma}]_{+}&=\delta_{nm},\\
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

:math:`F_{IJ}` と :math:`\Phi_{In\sigma}` の関係 (スピン平行も含めた場合)
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
   F_{IJ}=\sum_{n=1}^{\frac{N_{e}}{2}}\Big(\Phi_{In}\Phi_{Jn+1}-\Phi_{Jn}\Phi_{In+1}\Big).
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

として定義します。前節のように :math:`f_{ij}` (すなわち:math:`F`)が単一スレーター行列と関係づけられて
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

このセクションでは,
Power-Lanczos法での最適化パラメータ :math:`\alpha` の決定方法について述べます。
また,
ここではシングルステップのLanczos法を適用した後の物理量の計算についても説明します。

:math:`\alpha` の決定
~~~~~~~~~~~~~~~~~~~~~

最初に, 変分モンテカルロ法のサンプリングに関して簡単に説明します。
物理量 :math:`\hat{A}` は以下の手順で計算されます:

.. math::

   \begin{aligned}
   &\langle \hat{A}\rangle = \frac{\langle \phi| \hat{A}|\phi \rangle}{\langle \phi| \phi \rangle} = \sum_{x} \rho(x) F(x, {\hat{A}}),\\
   & \rho(x)=\frac{|\langle \phi|x\rangle|^2}{\langle \phi | \phi \rangle}, ~~~~F(x,  {\hat{A}}) =  \frac{\langle x| \hat{A}|\phi \rangle}{\langle x| \phi \rangle}.\end{aligned}

演算子の積 :math:`\hat{A}\hat{B}` を計算する場合には,
以下の二通りの方法があります。

.. math::

   \begin{aligned}
   &\langle \hat{A} \hat{B}\rangle = \sum_{x} \rho(x) F(x, {\hat{A}\hat{B}}),\\
   &\langle \hat{A} \hat{B}\rangle = \sum_{x} \rho(x) F^{\dagger}(x, {\hat{A})F(x, \hat{B}}).\end{aligned}

後述するように, 後者の表記の方が数値的には安定します。 例えば,
エネルギーの期待値の分散 :math:`\sigma^2=\langle (\hat{H}-\langle \hat{H}\rangle)^2\rangle` を考えてみると,
分散は以下の2通りの方法で計算できます。

.. math::

   \begin{aligned}
   \sigma^2 &=\sum_{x} \rho(x) F(x,  (\hat{H}-\langle \hat{H}\rangle)^2) = \sum_{x} \rho(x) F(x,  \hat{H}^2) - \left[ \sum_{x} \rho(x) F(x,  \hat{H})\right]^2 ,\\
   \sigma^2 &=\sum_{x} \rho(x) F^{\dagger}(x,  \hat{H}-\langle \hat{H}\rangle)F(x,  \hat{H}-\langle \hat{H}\rangle) \nonumber \\
   &= \sum_{x} \rho(x) F^{\dagger}(x,  \hat{H}) F(x, \hat{H})- \left[ \sum_{x} \rho(x) F(x,  \hat{H})\right]^2
   \end{aligned}

この定義から,
後者の方法では常に正の値となることが保証されているのに対して,
前者の方法では分散が正の値になることが必ずしも保証されないことが分かります。次に,
シングルステップでのpower-Lanczos波動関数 :math:`|\phi\rangle =(1+\alpha \hat{H}) |\psi \rangle`
に対するエネルギーの期待値とその分散を考えます。エネルギーは以下の式で計算されます：

.. math::

   \begin{aligned}
   E_{LS}(\alpha) =\frac{\langle \phi| \hat{H} |\phi\rangle}{\langle \phi|\phi\rangle}=\frac{h_1 + \alpha(h_{2(20)} + h_{2(11)}) + \alpha^2 h_{3(12)}}{1 + 2\alpha h_1 + \alpha^2 h_{2(11)}},
   \end{aligned}

ここで, :math:`h_1`, :math:`h_{2(11)},~h_{2(20)},` and
:math:`h_{3(12)}` を以下のように定義しました：

.. math::

   \begin{aligned}
   &h_1 =\sum_{x} \rho(x) F^{\dagger}(x,  \hat{H}),\\
   &h_{2(11)}=\sum_{x} \rho(x) F^{\dagger}(x,  \hat{H}) F(x, \hat{H}),\\
   &h_{2(20)}=\sum_{x} \rho(x) F^{\dagger}(x,  \hat{H}^2),\\
   &h_{3(12)}=\sum_{x} \rho(x) F^{\dagger}(x,  \hat{H})F(x,  \hat{H}^2).\end{aligned}

:math:`\frac{\partial E_{LS}(\alpha)}{\partial \alpha}=0` の条件から
:math:`\alpha` の二次方程式が導出され,
それを解くことで :math:`\alpha` が決定されます。
分散に関しても同様の手法で計算することが可能です。

物理量の計算
~~~~~~~~~~~~

最適化パラメータ :math:`\alpha` を用いることで,
演算子 :math:`\hat{A}` の期待値を以下の式から計算することが出来ます：

.. math::

   \begin{aligned}
   A_{LS}(\alpha) =\frac{\langle \phi| \hat{A} |\phi\rangle}{\langle \phi|\phi\rangle}=\frac{A_0 + \alpha(A_{1(10)} + A_{1(01)}) + \alpha^2 A_{2(11)}}{1 + 2\alpha h_1 + \alpha^2 h_{2(11)}},
   \end{aligned}

ここで, :math:`A_0`, :math:`A_{1(10)},~A_{1(01)},` and
:math:`A_{2(11)}` は

.. math::

   \begin{aligned}
   &A_0 =\sum_{x} \rho(x) F(x,  \hat{A}),\\
   &A_{1(10)}=\sum_{x} \rho(x) F^{\dagger}(x,  \hat{H}) F(x, \hat{A}),\\
   &A_{1(01)}=\sum_{x} \rho(x) F(x, \hat{A}\hat{H}),\\
   &A_{2(11)}=\sum_{x} \rho(x) F^{\dagger}(x,  \hat{H})F(x,  \hat{A}\hat{H}).\end{aligned}

として定義される変数を表します。プログラムでは,
この表式に基づき一体グリーン関数および二体グリーン関数の計算を行っています。
