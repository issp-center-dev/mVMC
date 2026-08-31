.. include:: ../../bib/ref.txt

.. _HowToExpert:

エキスパートモード入力ファイル書式
==================================

ここではmVMCで使用する詳細入力ファイル(\*def)のフォーマットに関して説明します。
入力ファイルの種別は以下の6つで分類されます。
なお、キーワードの後にある括弧内に記載されているファイル名はvmcdry.outにより作成されるファイル名を表します。

(1) リスト:

    **キーワード指定なし (namelist.def)**: 使用するinput
    fileの名前のリストを書きます。なお、ファイル名は任意に指定することができます。

(2) 基本パラメータ:

    **ModPara (modpara.def)**:
    計算時に必要な基本的なパラメーター(サイトの数、電子数など)を設定します。

    **LocSpin (locspn.def)**: 局在スピンの位置を設定します。

(3) ハミルトニアン:

    電子系の表式で記載されるハミルトニアン

    .. math::

       \begin{aligned}
       {\cal H}&={\cal H}_T+{\cal H}_U+{\cal H}_V+{\cal H}_H+{\cal H}_E+{\cal H}_P+{\cal H}_I,\\
       {\cal H}_T&={-}\sum_{i, j}\sum_{\sigma_1, \sigma_2}t_{ij\sigma_1\sigma_2} c_{i\sigma_1}^{\dagger}c_{j\sigma_2},\\
       {\cal H}_U&=\sum_{i} U_i n_ {i \uparrow}n_{i \downarrow},\\
       {\cal H}_V&=\sum_{i,j} V_{ij}n_ {i}n_{j},\\
       {\cal H}_H&={-}\sum_{i,j}J_{ij}^{\rm Hund} (n_{i\uparrow}n_{j\uparrow}+n_{i\downarrow}n_{j\downarrow}),\\
       {\cal H}_E&=\sum_{i,j}J_{ij}^{\rm Ex} (c_ {i \uparrow}^{\dagger}c_{j\uparrow}c_{j \downarrow}^{\dagger}c_{i  \downarrow}+c_ {i \downarrow}^{\dagger}c_{j\downarrow}c_{j \uparrow}^{\dagger}c_{i  \uparrow}),\\
       {\cal H}_P&=\sum_{i,j}J_{ij}^{\rm Pair} c_ {i \uparrow}^{\dagger}c_{j\uparrow}c_{i \downarrow}^{\dagger}c_{j  \downarrow},\\
       {\cal H}_I&=\sum_{i,j,k,l}\sum_{\sigma_1,\sigma_2, \sigma_3, \sigma_4}
       I_{ijkl\sigma_1\sigma_2\sigma_3\sigma_4}c_{i\sigma_1}^{\dagger}c_{j\sigma_2}c_{k\sigma_3}^{\dagger}c_{l\sigma_4},
       \end{aligned}


    について設定します。ここで、 :math:`n_{i \sigma}=c_{i\sigma}^{\dagger}c_{i\sigma}` は
    スピン :math:`\sigma` を持つサイト :math:`i` の電子密度演算子を、
    :math:`n_i=n_{i\uparrow}+n_{i\downarrow}` はサイト :math:`i` の電子密度演算子を
    それぞれ表します。
    各ハミルトニアンのパラメータは以下のファイルで指定します。

    **Trans (trans.def)**:
    :math:`{\cal H}_T` 内の :math:`t_{ij\sigma_1\sigma_2}` を指定します。

    **CoulombIntra (coulombintra.def)**:
    :math:`{\cal H}_U` 内の :math:`U_i` を指定します。

    **CoulombInter (coulombinter.def)**:
    :math:`{\cal H}_V` 内の :math:`V_{ij}` を指定します。

    **Hund (hund.def)**:
    :math:`{\cal H}_H` 内の :math:`J_{ij}^{\rm Hund}` を指定します。

    **Exchange (exchange.def)**:
    :math:`{\cal H}_E` 内の :math:`J_{ij}^{\rm Ex}` を指定します。

    **PairHop**:
    :math:`{\cal H}_P` 内の :math:`J_{ij}^{\rm Pair}` を指定します。

    **InterAll**:
    :math:`{\cal H}_I` 内の :math:`I_{ijkl\sigma_1\sigma_2\sigma_3\sigma_4}` を指定します。

    **NBodyInterAll (nbodyinterall.def)**:
    ハミルトニアンに含める可変次数の :math:`N` 体相互作用を指定します。

    .. rubric:: t-J模型の指定

    物理的なt-J模型を

    .. math::

       {\cal H}_{tJ} = -\sum_{i,j,\sigma} t_{ij}
       {\tilde c}_{i\sigma}^{\dagger}{\tilde c}_{j\sigma}
       + \sum_{i,j} J_{ij}\left({\boldsymbol S}_i\cdot{\boldsymbol S}_j
       - \frac{1}{4}n_i n_j\right)

    と書く場合、ここで :math:`{\tilde c}` は二重占有を除いた空間での
    演算子を表します。このt-J模型は ``InterAll`` ではなく、 ``Trans``、
    ``CoulombInter``、 ``Hund``、 ``Exchange`` の組み合わせで指定します。
    上記ハミルトニアンの符号規約では、正の物理的な結合 :math:`J_{ij}` に対して

    .. math::

       V_{ij}=-\frac{J_{ij}}{4}, \qquad
       J_{ij}^{\rm Hund}=-\frac{J_{ij}}{2}, \qquad
       J_{ij}^{\rm Ex}=-\frac{J_{ij}}{2}

    を指定します。また、t-J用の更新経路は ``modpara.def`` の
    ``NExUpdatePath`` で指定します。現状のt-J更新経路では
    ``BackFlow`` と ``LocSpin`` は非対応です。また、二重占有を許さないため、
    電子数がサイト数を超える入力は使用できません。 ``NExUpdatePath=4`` では
    spin hopping の行き先として少なくとも1つの空サイトが必要なため、
    ``Ncond < Nsite``、または ``2*Nelectron < Nsite`` が必要です。
    ``NExUpdatePath=5`` では ``Ncond <= Nsite``、または
    ``2*Nelectron <= Nsite`` が許されます。

(4) 最適化対象変分パラメータ:

    最適化する変分パラメータを指定します。変分波動関数は

    .. math::

       \begin{aligned}
       |\psi \rangle &= {\cal N}_{General RBM} {\cal P}_G{\cal P}_J{\cal P}_{SJ}{\cal P}_{d-h}^{(2)}{\cal P}_{d-h}^{(4)}{\cal L}^S{\cal L}^K{\cal L}^P |\phi_{\rm pair} \rangle,\\
       {\cal P}_G&=\exp\left[ \sum_i g_i n_{i\uparrow} n_{i\downarrow} \right],\\
       {\cal P}_J&=\exp\left[\frac{1}{2} \sum_{i\neq j} v_{ij} (n_i-1)(n_j-1)\right],\\
       {\cal P}_{SJ}&=\exp\left[\sum_{i<j} v^s_{ij} m_i m_j\right],\\
       {\cal P}_{d-h}^{(2)}&= \exp \left[ \sum_t \sum_{n=0}^2 (\alpha_{2nt}^d \sum_{i}\xi_{i2nt}^d+\alpha_{2nt}^h \sum_{i}\xi_{i2nt}^h)\right],\\
       {\cal P}_{d-h}^{(4)}&= \exp \left[ \sum_t \sum_{n=0}^4 (\alpha_{4nt}^d \sum_{i}\xi_{i4nt}^d+\alpha_{4nt}^h \sum_{i}\xi_{i4nt}^h)\right],\\
       {\cal N}_{\rm General RBM}&= \exp \left[ \sum_i a_{i\sigma} n_{i\sigma} \right] \prod_k^{N_h} \cosh \left[ b_k + \sum_{i\sigma} W_{i\sigma k} n_{i\sigma} \right],\\
       {\cal L}_S&=\frac{2S+1}{8 \pi^2}\int d\Omega P_s(\cos \beta) \hat{R}(\Omega),\\
       {\cal L}_K&=\frac{1}{N_s}\sum_{{\boldsymbol R}}e^{i {\boldsymbol K} \cdot{\boldsymbol R} } \hat{T}_{\boldsymbol R},\\
       {\cal L}_P&=\sum_{\alpha}p_{\alpha} \hat{G}_{\alpha},
       \end{aligned}

    で与えられます。ここで、
    :math:`\Omega=(\alpha, \beta, \gamma)` はオイラー角、
    :math:`\hat{R}(\Omega)` は回転演算子、
    :math:`P_S(x)` は :math:`S` 次のルジャンドル多項式、
    :math:`{\boldsymbol K}` は全運動量、
    :math:`\hat{T}_{\boldsymbol R}` は並進ベクトル
    :math:`{\boldsymbol R}` に対応する並進演算子、
    :math:`\hat{G}_{\alpha}` は格子の点群演算子、
    :math:`p_\alpha` はパリティをそれぞれ表します。
    Spin Jastrow因子では :math:`m_i=n_{i\uparrow}-n_{i\downarrow}` とします。
    ダブロン・ホロン相関因子に関する詳細は文献 [Tahara2008_ ]の説明を参照してください。
    また、一体部分は実空間のペア関数

    .. math::

       |\phi_{\rm pair} \rangle =
       \left[\sum_{i, j=1}^{N_s} \sum_{\sigma_1, \sigma_2}f_{i\sigma_1j\sigma_2}
       c_{i\sigma_1}^{\dagger}c_{j\sigma_2}^{\dagger} \right]^{N/2}|0 \rangle,

    を用いた波動関数で表されます。ここで :math:`N` は全電子数、
    :math:`N_s` は全サイト数です。
    最適化する変分パラメータは以下のファイルを用いて指定します
    ( :math:`{\cal L}_S` は **ModPara** ファイルでパラメータの指定をします)。

    **Gutzwiller (gutzwilleridx.def)**:
    :math:`{\cal P}_G` のうち、最適化の対象とする変分パラメータ :math:`g_i` を指定します。

    **Jastrow (jastrowidx.def)**:
    :math:`{\cal P}_J` のうち、最適化の対象とする変分パラメータ :math:`v_{ij}` を指定します。

    **SpinJastrow (spinjastrow.def)**:
    :math:`{\cal P}_{SJ}` のうち、最適化の対象とする変分パラメータ :math:`v^s_{ij}` を指定します。

    **BFRange (rangebf.def)**:
    BackFlow 相関で参照する各サイトの近傍リストと距離シェルを指定します。

    **BF (bf.def)**:
    BackFlow 相関のパラメータ group と最適化フラグを指定します。

    **DH2**:
    :math:`{\cal P}_{d-h}^{(2)}` で表される2サイトのダブロン・ホロン相関因子を指定します。

    **DH4**:
    :math:`{\cal P}_{d-h}^{(4)}` で表される4サイトのダブロン・ホロン相関因子を指定します。

    **GeneralRBM_PhysLayer**:
    :math:`{\cal N}_{\rm General RBM}` で表されるRBM相関因子のうち、最適化の対象とする変分パラメータ :math:`a_{i\sigma}` を指定します。

    **GeneralRBM_HiddenLayer**:
    :math:`{\cal N}_{\rm General RBM}` で表されるRBM相関因子のうち、最適化の対象とする変分パラメータ :math:`h_{k}` を指定します。

    **GeneralRBM_PhysHidden**:
    :math:`{\cal N}_{\rm General RBM}` で表されるRBM相関因子のうち、最適化の対象とする変分パラメータ :math:`W_{i\sigma k}` を指定します。

    **Orbital/OrbitalAntiParallel (orbitalidx.def)**:
    スピンが反平行のペア軌道 :math:`|\phi_{\rm pair} \rangle` を設定します。

    **OrbitalParallel**:
    スピンが平行のペア軌道 :math:`|\phi_{\rm pair} \rangle` を設定します。

    **OrbitalGeneral**:
    ペア軌道 :math:`|\phi_{\rm pair} \rangle` を設定します。

    **TransSym (qptransidx.def)**:
    運動量射影 :math:`{\cal L}_K` と格子対称性射影 :math:`{\cal L}_P` に関する指定を行います。

(5) 変分パラメータ初期値:

    変分パラメータに関する初期値を与えます。
    キーワード指定されない場合には :math:`0` が初期値として設定されます。

    **InGutzwiller**:
    :math:`{\cal P}_G` 内の変分パラメータ :math:`g_i` の初期値を設定します。

    **InJastrow**:
    :math:`{\cal P}_J` 内の変分パラメータ :math:`v_{ij}` の初期値を設定します。

    **InSpinJastrow**:
    :math:`{\cal P}_{SJ}` 内の変分パラメータ :math:`v^s_{ij}` の初期値を設定します。

    **InDH2**:
    :math:`{\cal P}_{d-h}^{(2)}` 内の2サイトのダブロン・ホロン相関因子
    :math:`\alpha_{2nt}^{d(h)}` の初期値を設定します。

    **InDH4**:
    :math:`{\cal P}_{d-h}^{(4)}` 内の4サイトのダブロン・ホロン相関因子
    :math:`\alpha_{4nt}^{d(h)}` の初期値を設定します。

    **InGeneralRBM_PhysLayer**:
    :math:`{\cal N}_{\rm General RBM}` で表されるRBM相関因子のうち、最適化の対象とする変分パラメータ :math:`a_{i\sigma}` の初期値を設定します。

    **InGeneralRBM_HiddenLayer**:
    :math:`{\cal N}_{\rm General RBM}` で表されるRBM相関因子のうち、最適化の対象とする変分パラメータ :math:`h_{k}` の初期値を設定します。

    **InGeneralRBM_PhysHidden**:
    :math:`{\cal N}_{\rm General RBM}` で表されるRBM相関因子のうち、最適化の対象とする変分パラメータ :math:`W_{i\sigma k}` の初期値を設定します。

    **InOrbital/InOrbitalAntiParallel**:
    ペア軌道 :math:`|\phi_{\rm pair} \rangle` の :math:`f_{i\uparrow j\downarrow}`
    に関する初期値を設定します。

    **InOrbitalParallel**:
    ペア軌道 :math:`|\phi_{\rm pair} \rangle` の :math:`f_{i\sigma j\sigma}`
    に関する初期値を設定します。

    **InOrbitalGeneral**:
    ペア軌道 :math:`|\phi_{\rm pair} \rangle` の
    :math:`f_{i\sigma j\sigma_1}` に関する初期値を設定します。

(6) 出力:

    **OneBodyG (greenone.def)**:出力する一体Green関数を指定します。

    **TwoBodyG (greentwo.def)**:出力する二体Green関数を指定します。

    **NBodyG (nbodyg.def)**:出力する :math:`N` 体相関関数を指定します。

    **Twist (twist.def)**:出力するTwist演算子を指定します。

(7) その他:

    **Lattice (lattice.def)**:サイト番号に対するユニットセル内部座標を含めた格子ベクトルを指定します。Twist演算子計算時に使用されます。

.. _InputFileList:

入力ファイル指定用ファイル(namelist.def)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

計算で使用する入力ファイル一式を指定します。ファイル形式に関しては、
以下のようなフォーマットをしています。

::

    ModPara  modpara.def
    LocSpin  zlocspn.def
    Trans    ztransfer.def
    InterAll zinterall.def
    LsTrans  zlstrans.def
    LsInterAll zlsinterall.def
    NBodyInterAll nbodyinterall.def
    Orbital orbitalidx.def
    OneBodyG zcisajs.def
    TwoBodyG	zcisajscktaltdc.def
    NBodyG nbodyg.def
    InUpdateWeight updateweight.def

ファイル形式
^^^^^^^^^^^^

[string01] [string02]

パラメータ
^^^^^^^^^^

-  [ string01 ]

   **形式 :** string型 (固定)

   **説明 :** キーワードを指定します。

-  [ string02 ]

   **形式 :** string型

   **説明 :** キーワードにひも付けられるファイル名を指定します(任意)。

使用ルール
^^^^^^^^^^

本ファイルを使用するにあたってのルールは以下の通りです。

-  キーワードを記載後、半角空白(複数可)を開けた後にファイル名を書きます。ファイル名は自由に設定できます。

-  ファイル読込用キーワードはTable[Table:Defs]により指定します。

-  必ず指定しなければいけないキーワードはModPara, LocSpin, Orbital,
   TransSymです。それ以外のキーワードについては、指定がない場合はデフォルト値が採用されます(変分パラメータについては最適化されず、固定する設定となります)。詳細は各ファイルの説明を参照してください。

-  各キーワードは順不同に記述できます。

-  指定したキーワード、ファイルが存在しない場合はエラー終了します。

-  :math:`\#` で始まる行は読み飛ばされます。


.. list-table::
   :widths: 30 70
   :header-rows: 1

   * - Keywords
     - 対応するファイルの概要
   * - ModPara :math:`^*`
     - 計算用のパラメータを指定します。
   * - LocSpin :math:`^*`
     - 局在・遍歴スピンを指定します。
   * - Trans
     - 一般的な一体相互作用を指定します。
   * - InterAll
     - 一般的な二体相互作用を指定します。
   * - LsTrans
     - 独立Lanczos演算子 :math:`H'` の一体部分を指定します。
   * - LsInterAll
     - 独立Lanczos演算子 :math:`H'` の二体部分を指定します。
   * - NBodyInterAll
     - ハミルトニアンに含める可変次数の :math:`N` 体相互作用を指定します。
   * - CoulombIntra
     - 内部クーロン相互作用を指定します。
   * - CoulombInter
     - サイト間クーロン相互作用を指定します。
   * - Hund
     - フント結合を指定します。
   * - PairHop
     - ペアホッピング相互作用を指定します。
   * - Exchange
     - 交換相互作用を指定します。
   * - Gutzwiller
     - 最適化するGutzwiller因子を設定します。
   * - Jastrow
     - 最適化する電荷Jastrow因子を指定します。
   * - SpinJastrow
     - 最適化するスピンJastrow因子を指定します。
   * - BFRange
     - BackFlow相関で参照する近傍リストと距離シェルを指定します。
   * - BF
     - BackFlow相関のパラメータgroupと最適化フラグを指定します。
   * - DH2
     - 最適化する2サイトダブロン・ホロン相関因子を指定します。
   * - DH4
     - 最適化する4サイトダブロン・ホロン相関因子を指定します。
   * - GeneralRBM_PhysLayer
     - 一般的なRBM相関因子のうち、最適化する物理層での変分パラメータを指定します。
   * - GeneralRBM_HiddenLayer
     - 一般的なRBM相関因子のうち、最適化する隠れ層での変分パラメータを指定します。
   * - GeneralRBM_PhysHidden
     - 一般的なRBM相関因子のうち、最適化する物理層と隠れ層を繋ぐ変分パラメータを指定します。
   * - Orbital :math:`^*`
     - 反平行のスピンを持つペア軌道因子を指定します。
   * - OrbitalAntiParallel
     - 反平行のスピンを持つペア軌道因子を指定します。
   * - OrbitalParallel
     - 平行のスピンを持つペア軌道因子を指定します。
   * - OrbitalGeneral
     - ペア軌道因子を指定します。
   * - TransSym :math:`^*`
     - 並進・格子対称演算子を設定します。
   * - InUpdateWeight
     - ローカル更新kernelの相対重みを任意指定します。
   * - InGutzwiller
     - Gutzwiller因子の初期値を設定します。
   * - InJastrow
     - 電荷Jastrow因子の初期値を設定します。
   * - InSpinJastrow
     - スピンJastrow因子の初期値を設定します。
   * - InDH2
     - 2サイトダブロン・ホロン相関因子の初期値を設定します。
   * - InDH4
     - 4サイトダブロン・ホロン相関因子の初期値を設定します。
   * - InGeneralRBM_PhysLayer
     - 一般的なRBM相関因子のうち、最適化する物理層での変分パラメータの初期値を設定します。
   * - InGeneralRBM_HiddenLayer
     - 一般的なRBM相関因子のうち、最適化する隠れ層での変分パラメータの初期値を設定します。
   * - InGeneralRBM_PhysHidden
     - 一般的なRBM相関因子のうち、最適化する物理層と隠れ層を繋ぐ変分パラメータの初期値を設定します。
   * - InOrbital
     - ペア軌道因子 :math:`f_{i\uparrow j\downarrow}` の初期値を設定します。
   * - InOrbitalAntiParallel
     - ペア軌道因子 :math:`f_{i\uparrow j\downarrow}` の初期値を設定します。
   * - InOrbitalParallel
     - ペア軌道因子 :math:`f_{i\sigma j\sigma}` の初期値を設定します。
   * - InOrbitalGeneral
     - ペア軌道因子 :math:`f_{i\sigma j\sigma'}` の初期値を設定します。
   * - OneBodyG
     - 出力する一体グリーン関数を指定します。
   * - TwoBodyG
     - 出力する二体グリーン関数を指定します。
   * - NBodyG
     - 出力する :math:`N` 体相関関数 :math:`\langle \prod_{a=1}^{N} c_{i_a\sigma_a}^{\dagger} c_{j_a\tau_a} \rangle` を指定します。
   * - Twist
     - 出力するTwist演算子を指定します。
   * - Lattice
     - サイト番号に対するユニットセル内部座標を含めた格子ベクトルを指定します。Twist演算子計算時に使用されます。

ModParaファイル (modpara.def)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

| 計算で使用するパラメータを指定します。以下のようなフォーマットをしています。

::

    --------------------
    Model_Parameters
    --------------------
    VMC_Cal_Parameters
    --------------------
    CDataFileHead  zvo
    CParaFileHead  zqp
    --------------------
    NVMCCalMode    0
    NLanczosMode   0
    NLanczosStep   1
    NLanczosEstimatorMode 0
    NLanczosSupportMode 0
    --------------------
    NDataIdxStart  1
    NDataQtySmp    1
    --------------------
    Nsite          16
    Nelectron      8
    NSPGaussLeg    1
    NSPStot        0
    NMPTrans       1
    NSROptItrStep  1200
    NSROptItrSmp   100
    DSROptRedCut   0.001
    DSROptStaDel   0.02
    DSROptStepDt   0.02
    NVMCWarmUp     10
    NVMCInterval   1
    NVMCSample     1000
    NExUpdatePath  0
    RndSeed        11272
    NSplitSize     1
    NStore         1
    NneuronGeneral 32

ファイル形式
^^^^^^^^^^^^

以下のように行数に応じ異なる形式をとります。

-  1 - 5行: ヘッダ(何が書かれても問題ありません)。

-  6行: [string01] [string02]

-  7行: [string03] [string04]

-  8行: ヘッダ(何が書かれても問題ありません)

-  9行以降: [string05] [int01] (もしくは[double01])

各項目の対応関係は以下の通りです。

-  [ string01 ]

   **形式 :** string型 (空白不可)

   **説明 :**
   アウトプットファイルのヘッダを表すためのキーワード。何を指定しても問題ありません。

-  [ string02 ]

   **形式 :** string型 (空白不可)

   **説明 :**
   アウトプットファイルのヘッダ。例えば、一体のGreen関数の出力ファイル名が **xxx\_cisajs.dat** として出力されます(xxxに指定した文字が記載)。

-  [ string03 ]

   **形式 :** string型 (空白不可)

   **説明 :**
   最適化された変分パラメータの出力ファイル名のヘッダを表すためのキーワード。何を指定しても問題ありません。

-  [ string04 ]

   **形式 :** string型 (空白不可)

   **説明 :**
   最適化された変分パラメータの出力ファイル名のヘッダ。最適化された変分パラメータが **xxx\_opt.dat** ファイルとして出力されます(xxxに指定した文字が記載)。

-  [ string05 ]

   **形式 :** string型 (固定)

   **説明 :** キーワードの指定を行います。

-  [ int01 ] ([ double01 ] )

   **形式 :** int (double)型 (空白不可)

   **説明 :** キーワードでひも付けられるパラメータを指定します。

使用ルール
^^^^^^^^^^

本ファイルを使用するにあたってのルールは以下の通りです。

-  9行目以降ではキーワードを記載後、半角空白(複数可)を開けた後に整数値を書きます。

-  9行目以降では“-”で始まる行は読み込まれません。



キーワード
^^^^^^^^^^

-  ``NVMCCalMode``

   **形式 :** int型 (デフォルト値 = 0)

   **説明 :** [0] 変分パラメータの最適化、[1] 1 体・2
   体のグリーン関数の計算。

-  ``NLanczosMode``

   **形式 :** int型 (デフォルト値 = 0)

   **説明 :** [0] 何もしない、[1] Single Lanczos Step
   でエネルギーまで計算、[2] Single Lanczos Step でエネルギー, 1 体・2
   体のグリーン関数まで計算(条件: 1, 2 は ``NVMCCalMode`` =
   1のみ使用可能。)

   ``OrbitalGeneral`` / FSZ、および BackFlow では ``NLanczosMode=1``
   のみ対応します。この場合、Hamiltonian は ``Transfer`` と
   number-operator 型の対角相互作用（``CoulombIntra``,
   ``CoulombInter``, ``Hund``）に限定されます。``NLanczosMode=2``、
   ``NBodyG`` / ``NBodyInterAll`` の Lanczos 補正は未対応です。

-  ``NLanczosStep``

   **形式 :** int型 (デフォルト値 = 1)

   **説明 :** Power-Lanczos計算の次数を指定します。[1] は従来の
   single-step計算および出力を変更せずに使用します。[2] は3次元Krylov
   部分空間で
   :math:`(c_0+c_1\hat{H}+c_2\hat{H}^2)|\psi\rangle`
   を最適化します。

   2nd stepには ``NVMCCalMode=1`` と ``NLanczosMode=1`` が必要です。
   実数または複素数のnon-FSZペア軌道について、次の2つのmodel classに
   対応します。

   * electronic classは ``NExUpdatePath=0`` を使用し、spinを保存する
     ``Transfer`` 項（:math:`\sigma_1=\sigma_2`）とnumber-operator型の
     対角相互作用（``CoulombIntra``, ``CoulombInter``, ``Hund``）に
     対応します。既存のスピン射影・運動量射影を使用できます。
   * pure spin-1/2 classは ``NExUpdatePath=2`` を使用し、
     ``NLocalSpin=Nsite=2*Ne``、``NTransfer=0``、``NQPFull=1`` を
     必須とします。Hamiltonianにはnumber-operator型の対角相互作用と
     ``Exchange`` 項を指定できます。固定 ``2Sz=0`` sectorの
     Heisenberg/XXZ、非一様結合、Ising極限を含みます。

   ``OrbitalGeneral`` / FSZ、BackFlow、RBM、``PairHop``、
   ``InterAll``、``NBodyInterAll``、``NBodyG`` はsampling開始前に
   エラー終了します。``Exchange`` はelectronic classでは未対応、
   ``Transfer`` はpure-spin classでは未対応です。局在spinと遍歴電子の
   mixed系、および ``NQPFull>1`` のpure-spin量子射影も現scope外です。
   2nd stepでは ``reweight=1`` もエラー終了します。reweightは基底
   波動関数のsupport上のweightだけを変え、基底振幅が0の配置を回復できません。
   2nd stepのGreen関数は計算しません。
   Gutzwillerの ``ProjRatio`` 経路は独立ED oracleで値を検証しています。
   非自明な :math:`k=\pi` 運動量射影は、real/complexの独立projected ED
   value oracleで検証しています。real oracleはserial、MPI、OpenMP経路を
   対象とし、complex oracleは現時点でserial経路を対象とします。
   非自明なスピン射影はruntimeおよび構造のsmoke testまでで、
   projected EDによる独立な値oracleは現時点ではありません。

   electronic classで ``NQPFull=1`` の場合、:math:`H^3` のlocal powerは
   ``Transfer`` 項に対する
   outer/inner二重loopで評価するため、支配的なoperator数は
   :math:`O(N_{\mathrm{Transfer}}^2)` で増加します。非自明な量子射影
   （``NQPFull>1``）では、射影された :math:`H\,CACA` 行列要素の直接縮約経路に
   もう1段の ``Transfer`` loopが加わります。そのため支配的なoperator数は
   :math:`O(N_{\mathrm{Transfer}}^3)` となり、``GreenFuncN`` の縮約量は
   ``NQPFull`` にも比例します。

   pure-spin classでは、各 ``Exchange`` bondを2つの向き付きspin-exchange
   演算子へ展開します。:math:`H^3` の直接縮約は、深さ1、2、3のactive pathを
   それぞれ2、4、6次の ``GreenFuncN`` で評価します。sampleあたり
   :math:`A` 個の向き付きexchangeがactiveなら、各深さのcall数上限は
   :math:`A`、:math:`A^2`、:math:`A^3` で、深さ3が支配します。
   このclassは現時点で ``NQPFull=1`` のみ対応します。

   したがって、射影により2nd stepの測定時間が大幅に増える場合があります。
   実サイズの計算前に、実計算と同じ ``NSPGaussLeg`` と ``NMPTrans`` を使った
   小さい系で測定し、性能判断には対象Linux/HPC環境のbaselineを使用してください。
   source treeにある非CIの手動probe
   ``test/python/lanczos2_cost_probe.py`` は ``--nspgaussleg`` と
   ``--nmptrans`` のgridを受け取り、JSON出力へ ``NQPFull``、Timer 41、
   Timer 95を記録します。

   local powerまたはmomentの非有限値、不正なoverlap行列、一般化固有値
   solverの失敗、出力エラーを検出すると診断を表示して終了します。
   非有限値のエラー経路を検証する環境変数
   ``MVMC_LANCZOS2_TEST_NONFINITE_SAMPLE`` はテスト専用です。
   ``Testing=ON`` のbuildだけにコンパイルされ、通常の
   ``Testing=OFF`` buildには故障注入コードが含まれません。
   overlap行列が正定値でない場合は、相対対角regularization
   :math:`10^{-12}` を加えて1回だけ再試行し、結果を ``solve_flag`` に
   記録します。

-  ``NLanczosEstimatorMode``

   **形式 :** int型 (デフォルト値 = 0)

   **説明 :** [0] 従来のPower-Lanczos計算、[1] 数値安定化した
   energy/variance estimatorを使用します。[1] はexpert modeで明示的に指定し、
   ``NVMCCalMode=1``、``NLanczosMode=1``、``NLanczosStep=1`` または ``2``
   を必要とします。係数決定用chainと最終energy用chainは独立です。
   FSZ、BackFlow、RBM、reweight、およびGreen関数は対象外です。
   出力は ``xxx_pl_out_yyy.dat`` のみです。``NVMCSample`` は8以上かつ、
   各blockに2 sample以上が入る4から16までのblock数で割り切れる値とし、
   16の倍数を推奨します。

   ``LsTrans`` または ``LsInterAll`` を指定すると、基底
   :math:`\{\Psi,H'\Psi\}` を使う独立演算子経路を選択します。この経路はさらに
   実変分パラメータ、``NLanczosStep=1``、``NExUpdatePath=0``、physical modeを
   必要とします。物理 :math:`H` と独立 :math:`H'` は ``Trans`` と
   ``InterAll`` の行だけで構成してください。FSZ、BackFlow、RBM、reweight、
   native model coupling、special update、Green関数測定は対象外です。
   ``NQPFull`` はruntime入力値であり、本機能では固定しません。

-  ``NLanczosSupportMode``

   **形式 :** int型 (デフォルト値 = 0)

   **説明 :** 独立なpower-Lanczos support監査の動作を指定します。[0] は
   production用strict modeです。``M02`` と ``M11``、2nd stepではさらに
   ``M03`` と ``M12`` が不一致の場合、または有効sampleが32未満で判定不能の
   場合にLanczos solverの前で停止します。relative differenceが0.5以上で
   scoreが4.5未満の未解決な大差も判定不能とします。reweighting weightが厳密に
   0のsampleは有効sample数から除外します。[1] は明示的なexperimental互換
   modeです。従来のsupport制限estimatorを続行し、監査失敗を
   ``xxx_ls_support_yyy.dat`` に ``biased-diagnostic-only`` と記録します。
   監査の通過は必要条件ですが、全ての高次Krylov supportが完全であることの
   十分条件ではありません。

-  ``NDataIdxStart``

   **形式 :** int型 (デフォルト値 = 0)

   **説明 :** 出力ファイルの付加番号。 ``NVMCCalMode`` = 0
   の場合は ``NDataIdxStart`` が出力され、 ``NVMCCalMode`` = 1
   の場合は、 ``NDataIdxStart`` から連番で ``NDataQtySmp`` 個のファイルを出力します。

-  ``NDataQtySmp``

   **形式 :** int型 (デフォルト値 = 1)

   **説明 :** 出力ファイルのセット数。 ``NVMCCalMode`` = 1
   の場合に使用します。

-  ``Nsite``

   **形式 :** int型 (1以上、必須)

   **説明 :** サイト数を指定する整数。

-  ``Nelectron``

   **形式 :** int型 (1以上、必須)

   **説明 :** 電子のペア数(電子数は2 ``Nelectron`` で与えられる)。

-  ``NCond``

   **形式 :** int型 (0以上)

   **説明 :** 伝導電子の数。

-  ``2Sz``

   **形式 :** int型

   **説明 :** 2 :math:`S_z` の値。電子がペアを組むため,
   2 :math:`S_z` は偶数で指定する必要がある。

-  ``NSPGaussLeg``

   **形式 :** int型 (1以上、デフォルト値 = 8)

   **説明 :**
   スピン量子数射影の :math:`\beta` 積分( :math:`S^y` 回転)のGauss-Legendre求積法の分点数。

-  ``NSPStot``

   **形式 :** int型 (0以上、デフォルト値 = 0)

   **説明 :** スピン量子数。

-  ``NMPTrans``

   **形式 :** int型 (デフォルト値 = 1)

   **説明 :**
   ``NMPTrans`` の絶対値で並進・格子対称性の量子数射影の個数を指定する。負の場合は反周期境界条件を与える。
   TransSymファイルで指定した重みで上から ``abs(NMPTrans)`` 個まで使用する。
   0 は指定できず、``abs(NMPTrans)`` は TransSym の ``NQPTrans`` 以下でなければならない。
   射影を行わない場合は1（反周期の恒等patternでは -1）に設定する必要があります。

-  ``NSROptItrStep``

   **形式 :** int型 (1以上、デフォルト値 = 1000)

   **説明 :** SR
   法で最適化する場合の全ステップ数。 ``NVMCCalMode`` =0の場合のみ使用されます。

-  ``NSROptItrSmp``

   **形式 :** int型 (1以上数、デフォルト値 = ``NSROptItrStep``/10)

   **説明 :**
   ``NSROptItrStep`` ステップ中、最後の ``NSROptItrSmp`` ステップでの各変分パラメータの平均値を最適値とする。 ``NVMCCalMode`` =0の場合のみ使用されます。

-  ``DSROptRedCut``

   **形式 :** double型 (デフォルト値 = 0.001)

   **説明 :** SR
   法安定化因子。手法論文[Tahara2008_ ]の :math:`\varepsilon_{\rm wf}` に対応。

-  ``DSROptStaDel``

   **形式 :** double型 (デフォルト値 = 0.02)

   **説明 :** SR
   法安定化因子。手法論文[Tahara2008_ ]の :math:`\varepsilon` に対応。

-  ``DSROptStepDt``

   **形式 :** double型

   **説明 :**
   SR法で使用する刻み幅。手法論文[Tahara2008_ ]の :math:`\Delta t` に対応。

-  ``NSROptCGMaxIter``

   **形式 :** int型 (デフォルト値 = 0)

   **説明 :** SR-CG法での、CG法の繰り返し回数の上限。
   0以下を指定した場合、最大で :math:`S`
   行列のサイズの数だけ実行するようになります。 ``NSRCG``!=0
   の場合のみ使用されます。

-  ``DSROptCGTol``

   **形式 :** double型 (デフォルト値 = 1.0e-10)

   **説明 :**
   SR-CG法での、CG法の収束判定条件。残差ベクトルの要素の自乗平均平方根がこの値以下になったらCG
   法を終了します。 ``NSRCG``!=0 の場合のみ使用されます。

-  ``NVMCWarmUp``

   **形式 :** int型 (1以上、デフォルト値=10)

   **説明 :** マルコフ連鎖の空回し回数。

-  ``NVMCInterval``

   **形式 :** int型 (1以上、デフォルト値=1)

   **説明 :** サンプル間のステップ間隔。ローカル更新を ``Nsite`` ×
   ``NVMCInterval`` 回行います。

-  ``NVMCSample``

   **形式 :** int型 (1以上、デフォルト値=1000)

   **説明 :** 期待値計算に使用するサンプル数。

-  ``NExUpdatePath``

   **形式 :** int型 (0以上)

   **説明 :** ローカル更新の種類を指定します。
   0: HOPPING、1: EXCHANGE または HOPPING、2: EXCHANGE、
   3: KondoGC用（HOPPING または EXCHANGE/LOCALSPINFLIP）、
   4: tJ用 SPINHOPPING、5: tJ用（EXCHANGE または SPINHOPPING）、
   6: pair hoppingによるdoublon-onlyサンプリング。
   4と5の違いはt-J空間でのローカル更新経路の違いであり、
   :math:`S_z` 保存・非保存の切り替えではありません。スピン量子数の指定は、
   固定 :math:`S_z` 計算では ``2Sz`` などで別途行います。
   t-J用の4と5では ``BackFlow`` と ``LocSpin`` は非対応です。また、
   二重占有を許さないため、電子数が ``Nsite`` を超える入力は使用できません。
   さらに、 ``NExUpdatePath=4`` では spin hopping の行き先として
   少なくとも1つの空サイトが必要なため、 ``Ncond < Nsite``、または
   ``2*Nelectron < Nsite`` が必要です。 ``NExUpdatePath=5`` では
   ``Ncond <= Nsite``、または ``2*Nelectron <= Nsite`` が許されます。
   ``NExUpdatePath=6`` では各サイトの状態を空状態またはdoublon状態、
   すなわち :math:`(n_{\uparrow}, n_{\downarrow})=(0,0)` または
   :math:`(1,1)` に制限します。このモードでは ``0 < Ne < Nsite`` と
   反平行スピンの ``Orbital`` 入力が必要です。現状では ``LocSpin``、
   ``BackFlow``、RBM、``OrbitalGeneral``/FSZ入力には対応していません。

   ``InUpdateWeight``を指定しない場合、``NExUpdatePath=2``の従来selectorは
   変更されません。固定 :math:`S_z` または非FSZ計算ではEXCHANGEのみ、
   ``OrbitalGeneral``かつ``2Sz=-1``ではEXCHANGEとLOCALSPINFLIPを等確率で
   選びます。重み付きselectorとPAIRSPINFLIPについては後述の
   ``updateweight.def``を参照してください。

-  ``RndSeed``

   **形式 :** int型

   **説明 :** 乱数の初期seed。MPI 並列では各計算機に ``RndSeed`` +my
   rank+1 で初期seed が与えられます。

-  ``NSplitSize``

   **形式 :** int型 (1以上、デフォルト値=1)

   **説明 :** 1つの内部MPI並列グループに含まれるMPIプロセス数。
   並列化の対象は計算段階で異なります。

   サンプル生成段階 (``VMCMakeSample*``) では、パラメータ最適化計算と
   物理量計算のどちらでも、量子数射影の添字を ``NQPFull`` 個の分点に
   わたってこのグループ内のプロセスで分割します。ここで
   ``NQPFull = NSPGaussLeg * NMPTrans * NQPOptTrans`` であり、
   ``NQPOptTrans`` は ``OptTrans`` mode を使わない通常時は1です。

   サンプル生成後の主評価段階 (``VMCMainCal*``) では、射影分点ではなく
   モンテカルロサンプル (``NVMCSample``) を同じグループ内のプロセスで
   分割します。各プロセスは担当サンプルについて全ての射影分点を評価
   します。この段階は、パラメータ最適化計算ではSR量の集計、物理量計算
   では物理量の集計に使われます。

   ``NSplitSize`` は1以上である必要があります。全MPIプロセス数が
   ``NSplitSize`` で割り切れない場合、最後のグループが小さくなり
   load imbalance が起こり得ます。直接SRでは ``NStore=0`` と
   ``NStore=1`` のどちらも ``NSplitSize > 1`` に対応しています。
   SR-CG の stored :math:`O` 行列ベクトル積は内部MPI分割に対応していないため、
   ``NSplitSize > 1`` の場合は ``NSRCG=0`` を指定してください。

-  ``NStore``

   **形式 :** int型 (0もしくは1、デフォルト値=1)

   **説明 :**
   期待値 :math:`\langle O_k O_l \rangle` を計算するとき行列-行列積にして高速化するオプション
   (1で機能On、モンテカルロサンプリング数に応じてメモリの消費が増大します [3]_)。
   直接SR (``NSRCG=0``) では ``NSplitSize > 1`` と併用できます。

-  ``NSRCG``

   **形式 :** int型 (0もしくは1、デフォルト値=0)

   **説明 :** SR法で連立一次方程式 :math:`Sx=g`
   を解くときに、 :math:`S`
   を陽に構築せずに解くことでメモリを削減する [4]_ オプション[NeuscammanUmrigarChan_ ](1で機能On,
   ``NStore`` は1に固定されます)。 ``NSRCG`` は ``NSplitSize > 1`` と併用できません。
   この制限は内部MPI分割に固有のものです。 ``NSplitSize=1`` の場合、SR-CG は従来通り
   stored :math:`O` 経路を内部で使用します。

-  ``useDiagScale``

   **形式 :** int型 (0もしくは1、デフォルト値=0)

   **説明 :** SR法での連立一次方程式 :math:`Sx=g` をCG法により解く際に、 Point Jacobi法 (:math:`S` 行列の対角スケーリング)による前処理付きCG法を使用するオプション(1で機能ON, ``NSRCG=1`` である必要がある)。

   ``NSRCG=2`` を指定した場合も受理され、内部では
   ``NSRCG=1`` かつ ``useDiagScale=1`` として扱われます。

-  ``NSRCGFallback``

   **形式 :** int型 (0もしくは1、デフォルト値=0)

   **説明 :** 選択されたSR-CGソルバが収束しない、または数値不安定になった場合に、
   もう一方のSR-CGソルバで再試行するオプション(0で無効、1で有効)。
   通常CGはDiagScale-CGへ、DiagScale-CGは通常CGへフォールバックします。

-  ``NSRCGAbortOnFail``

   **形式 :** int型 (0もしくは1、デフォルト値=1)

   **説明 :** SR-CGが、必要ならフォールバックも試した後で失敗した場合に
   計算を終了するオプション(0で警告を出して近似解で継続、1で終了)。

   デフォルトは ``NSRCGFallback=0`` かつ ``NSRCGAbortOnFail=1`` です。
   この設定では、SR-CGが収束しない、または数値不安定になった場合に、
   入力または数値条件の失敗として計算を終了します。
   overlap行列の統計ノイズが大きい、または条件が悪い場合は、
   ``NVMCSample`` を増やす、SR-CGの収束判定値や反復回数上限を調整する、
   あるいは ``useDiagScale`` を有効にすることを検討してください。
   ``useDiagScale=1`` はPoint Jacobi前処理付きCGを選択します。
   収束性の改善に有効な場合があります。

-  ``RescaleSmat``

   **形式 :** int型 (0もしくは1、デフォルト値=0)

   **説明 :** SR-CGで :math:`Sx=g` を解く前に、Slater関連ブロックを
   リスケーリングするオプション(1で機能ON)。
   ``RescaleSmat=1`` を使う場合は ``NSRCG=1`` が必要です
   (リスケーリングはSR-CG経路でのみ適用されます)。
   ``ModPara`` でのSR-CGの典型設定は次の通りです。

   ::

       NSRCG = 1
       useDiagScale = 1
       RescaleSmat = 1

-  ``NneuronGeneral``

   **形式 :** int型 (デフォルト値=0)

   **説明 :** RBMの隠れ層にあるニューロン数 :math:`N_{\rm General RBM}` を指定する整数。

UpdateWeight指定ファイル(updateweight.def)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``NExUpdatePath=2``で用いるローカル更新kernelの相対重みを指定する任意入力です。
``namelist.def``へ``InUpdateWeight updateweight.def``を追加すると有効になります。
ファイルを指定しない場合、従来のselector、乱数消費順、``zvo_time``の列構成を
そのまま維持します。

::

    ================================
    NUpdateWeight 3
    ================================
    ========UpdateType Weight=======
    ================================
    Exchange          1.0
    LocalSpinFlip     1.0
    PairSpinFlip      1.0

2行目にdata行数を書き、6行目以降にkernel名と有限・非負の重みを指定します。
行順は任意で、省略した対応済みkernelの重みは0です。重複名、未知の名前、余分な列、
正でない重み合計、またはpath 2以外での使用はエラー終了します。重みは内部で正規化し、
入力読込時に表示します。

``LocalSpinFlip``には``OrbitalGeneral``と``2Sz=-1``が必要です。
この重みを0にすると、up-spin数のparityを意図的に保存します。宣言した固定parity
sectorだけをsampleする場合には正しい設定です。target分布が両parity sectorに
nonzero weightを持つ場合、``LocalSpinFlip``の正の重みは必要条件であり、実際に
acceptされたlocal-spin-flipも確認する必要があります。proposal重みが正であること
だけではergodicityを保証しません。

このファイルを有効にすると、正の重みを持つkernelが1個だけでもkernel選択ごとに
乱数を1個消費します。したがって、選ばれるkernel列が同じでもweighted runの乱数列は
legacy runと一致しません。乱数消費順の互換性を保証するのは``InUpdateWeight``を
省略した場合だけです。

``PairSpinFlip``は、初期spinが同じ2個の異なる局在spinを同時に反転します。さらに
``Ncond=0``、``NLocalSpin=Nsite``、全siteで``LocSpin=1``を要求します。初期spinが
逆向きのpairは棄却proposalとして数え、現在stateに応じた重みの再正規化は行いません。
これによりproposal確率をstate非依存かつ逆更新と対称に保ちます。このファイルが
有効な場合、PAIRSPINFLIPのproposal数とacceptanceを``zvo_time``末尾へ追記します。

LocSpin指定ファイル(locspn.def)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

局在スピンを指定します。以下のようなフォーマットをしています。

::

    ================================
    NlocalSpin     6
    ================================
    ========i_0LocSpn_1IteElc ======
    ================================
        0      1
        1      0
        2      1
        3      0
        4      1
        5      0
        6      1
        7      0
        8      1
        9      0
       10      1
       11      0

ファイル形式
^^^^^^^^^^^^

以下のように行数に応じ異なる形式をとります。

-  1行: ヘッダ(何が書かれても問題ありません)。

-  2行: [string01] [int01]

-  3-5行: ヘッダ(何が書かれても問題ありません)。

-  6行以降: [int02] [int03]

パラメータ
^^^^^^^^^^

-  [ string01 ]

   **形式 :** string型 (空白不可)

   **説明 :** 局在スピンの総数を示すキーワード(任意)。

-  [ int01 ]

   **形式 :** int型 (空白不可)

   **説明 :** 局在スピンの総数を指定する整数。

-  [ int02 ]

   **形式 :** int型 (空白不可)

   **説明 :**
   サイト番号を指定する整数。0以上 ``Nsite`` 未満で指定します。

-  [ int03 ]

   **形式 :** int型 (空白不可)

   | **説明 :** 遍歴電子か局在スピンかを指定する整数(0: 遍歴電子, 1:
     局在スピン)。

使用ルール
^^^^^^^^^^

本ファイルを使用するにあたってのルールは以下の通りです。

-  行数固定で読み込みを行う為、ヘッダの省略はできません。

-  [ int01 ] と [ int03 ] で指定される局在電子数の総数が異なる場合はエラー終了します。

-  [ int02 ] の総数が全サイト数と異なる場合はエラー終了します。

-  [ int02 ] が全サイト数以上もしくは負の値をとる場合はエラー終了します。

Trans指定ファイル(trans.def)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

一般的な一体相互作用をハミルトニアンに付け加え、ハミルトニアン中の相互作用パラメータ :math:`t_{ij\sigma_1\sigma2}` を指定します。付け加える項は以下で与えられます。

.. math::

   \begin{aligned}
   {\cal H}_T =-\sum_{ij\sigma_1\sigma_2} t_{ij\sigma_1\sigma_2}c_{i\sigma_1}^{\dagger}c_{j\sigma_2}\end{aligned}

以下にファイル例を記載します。

::

    ========================
    NTransfer      24
    ========================
    ========i_j_s_tijs======
    ========================
        0     0     2     0   1.000000  0.000000
        2     0     0     0   1.000000  0.000000
        0     1     2     1   1.000000  0.000000
        2     1     0     1   1.000000  0.000000
        2     0     4     0   1.000000  0.000000
        4     0     2     0   1.000000  0.000000
        2     1     4     1   1.000000  0.000000
        4     1     2     1   1.000000  0.000000
        4     0     6     0   1.000000  0.000000
        6     0     4     0   1.000000  0.000000
        4     1     6     1   1.000000  0.000000
        6     1     4     1   1.000000  0.000000
        6     0     8     0   1.000000  0.000000
        8     0     6     0   1.000000  0.000000
    …

ファイル形式
^^^^^^^^^^^^

以下のように行数に応じ異なる形式をとります。

-  1行: ヘッダ(何が書かれても問題ありません)。

-  2行: [string01] [int01]

-  3-5行: ヘッダ(何が書かれても問題ありません)。

-  6行以降: [int02]  [int03]  [int04]  [int05]  [double01]  [double02]

パラメータ
^^^^^^^^^^

-  [ string01 ]

   **形式 :** string型 (空白不可)

   **説明 :** 定義するパラメータの総数のキーワード名を指定します(任意)。

-  [ int01 ]

   **形式 :** int型 (空白不可)

   **説明 :** 定義するパラメータの総数を指定します。

-  [ int02 ], [ int04 ]

   **形式 :** int型 (空白不可)

   **説明 :**
   サイト番号を指定する整数。0以上 ``Nsite`` 未満で指定します。

-  [ int03 ], [ int05 ]

   **形式 :** int型 (空白不可)

   | **説明 :** スピンを指定する整数。
   | 0: アップスピン
   | 1: ダウンスピン
   | を選択することが出来ます。

-  [ double01 ]

   **形式 :** double型 (空白不可)

   **説明 :** :math:`t_{ij\sigma_1\sigma_2}` の実部を指定します。

-  [ double02 ]

   **形式 :** double型 (空白不可)

   **説明 :** :math:`t_{ij\sigma_1\sigma_2}` の虚部を指定します。

使用ルール
^^^^^^^^^^

本ファイルを使用するにあたってのルールは以下の通りです。

-  行数固定で読み込みを行う為、ヘッダの省略はできません。

-  空行は許されません。

-  [ int01 ] と定義されているTrasferの総数が異なる場合はエラー終了します。

-  [ int02 ]-[ int05 ] を指定する際、範囲外の整数を指定した場合はエラー終了します。

-  Hamiltonianがエルミートという制限から :math:`t_{ij\sigma_1\sigma_2}=t_{ji\sigma_2\sigma_1}^{\dagger}` の関係を満たす必要があります。

InterAll指定ファイル
~~~~~~~~~~~~~~~~~~~~

一般的な二体相互作用をハミルトニアンに付け加え、ハミルトニアン中の相互作用パラメータを指定します。
付け加える項は以下で与えられます。

.. math::

   {\cal H}_{I}=\sum_{i,j,k,l}\sum_{\sigma_1,\sigma_2, \sigma_3, \sigma_4}
   I_{ijkl\sigma_1\sigma_2\sigma_3\sigma_4}c_{i\sigma_1}^{\dagger}c_{j\sigma_2}c_{k\sigma_3}^{\dagger}c_{l\sigma_4}

以下にファイル例を記載します。

::

    ======================
    NInterAll      36
    ======================
    ========zInterAll=====
    ======================
    0    0    0    1    1    1    1    0   0.50  0.0
    0    1    0    0    1    0    1    1   0.50  0.0
    0    0    0    0    1    0    1    0   0.25  0.0
    0    0    0    0    1    1    1    1  -0.25  0.0
    0    1    0    1    1    0    1    0  -0.25  0.0
    0    1    0    1    1    1    1    1   0.25  0.0
    2    0    2    1    3    1    3    0   0.50  0.0
    2    1    2    0    3    0    3    1   0.50  0.0
    2    0    2    0    3    0    3    0   0.25  0.0
    2    0    2    0    3    1    3    1  -0.25  0.0
    2    1    2    1    3    0    3    0  -0.25  0.0
    2    1    2    1    3    1    3    1   0.25  0.0
    4    0    4    1    5    1    5    0   0.50  0.0
    4    1    4    0    5    0    5    1   0.50  0.0
    4    0    4    0    5    0    5    0   0.25  0.0
    4    0    4    0    5    1    5    1  -0.25  0.0
    4    1    4    1    5    0    5    0  -0.25  0.0
    4    1    4    1    5    1    5    1   0.25  0.0
    …

ファイル形式
^^^^^^^^^^^^

以下のように行数に応じ異なる形式をとります。

-  1行: ヘッダ(何が書かれても問題ありません)。

-  2行: [string01] [int01]

-  3-5行: ヘッダ(何が書かれても問題ありません)。

-  6行以降:
   [int02] [int03] [int04] [int05] [int06] [int07] [int08] [int09] [double01] [double02]

パラメータ
^^^^^^^^^^

-  [ string01 ]

   **形式 :** string型 (空白不可)

   **説明 :** 二体相互作用の総数のキーワード名を指定します(任意)。

-  [ int01 ]

   **形式 :** int型 (空白不可)

   **説明 :** 二体相互作用の総数を指定します。

-  [ int02 ], [ int04 ],
   [ int06 ], [ int08 ]

   **形式 :** int型 (空白不可)

   **説明 :**
   サイト番号を指定する整数。0以上 ``Nsite`` 未満で指定します。

-  [ int03 ], [ int05 ],
   [ int07 ], [ int09 ]

   **形式 :** int型 (空白不可)

   **説明 :** スピンを指定する整数。

   0: アップスピン

   1: ダウンスピン

   を選択することが出来ます。

-  [ double01 ]

   **形式 :** double型 (空白不可)

   **説明 :**
   :math:`I_{ijkl\sigma_1\sigma_2\sigma_3\sigma_4}` の実部を指定します。

-  [ double02 ]

   **形式 :** double型 (空白不可)

   **説明 :**
   :math:`I_{ijkl\sigma_1\sigma_2\sigma_3\sigma_4}` の虚部を指定します。

使用ルール
^^^^^^^^^^

本ファイルを使用するにあたってのルールは以下の通りです。

-  行数固定で読み込みを行う為、ヘッダの省略はできません。

-  ハミルトニアンがエルミートという制限から :math:`I_{ijkl\sigma_1\sigma_2\sigma_3\sigma_4}=I_{lkji\sigma_4\sigma_3\sigma_2\sigma_1}^{\dagger}` の関係を満たす必要があります。

-  [ int01 ] と定義されているInterAllの総数が異なる場合はエラー終了します。

-  [ int02 ]-[ int09 ] を指定する際、範囲外の整数を指定した場合はエラー終了します。

独立Lanczos演算子指定ファイル(lstrans.def, lsinterall.def)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``LsTrans`` と ``LsInterAll`` は、物理Hamiltonianとは独立な演算子 :math:`H'` を
定義します。ファイル形式と符号規約は ``Trans`` と ``InterAll`` と同じです。
``LsTrans`` の1行は :math:`-t c^\dagger c`、``LsInterAll`` の1行は記載した
CACA順で作用します。どちらか一方は省略できますが、両者を合わせた演算子は非空、
実係数、Hermitianで、粒子数と :math:`S_z` を保存する必要があります。mVMCは
行ごとの有限性・index・sector条件を検査し、Hermitian closureは入力生成側の責務です。

NBodyInterAll指定ファイル(nbodyinterall.def)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

可変次数の :math:`N` 体相互作用をハミルトニアンに付け加えます。付け加える項は以下で与えられます。

.. math::

   {\cal H}_{N} =
   \sum_{\alpha} V_{\alpha}
   \prod_{a=1}^{N_{\alpha}}
   c_{i_{\alpha a}\sigma_{\alpha a}}^{\dagger}
   c_{j_{\alpha a}\tau_{\alpha a}}

各データ行に書かれた因子の順序が演算子の順序として使われます。エルミート共役項は自動生成されないため、
目的のハミルトニアンに必要な項はすべて入力ファイルで指定してください。

::

    =============================================
    NNBodyInterAll   3
    =============================================
    ======== NBodyInterAll interactions =========
    =============================================
        1     0     0     2     0    0.37    0.21
        2     0     0     2     0     5     1     1     1   -0.28    0.16
        3     0     0     2     0     1     0     3     0     5     1     1     1    0.22   -0.17

ファイル形式
^^^^^^^^^^^^

以下のように行数に応じ異なる形式をとります。

-  1行: ヘッダ(何が書かれても問題ありません)。

-  2行: [string01] [int01]

-  3-5行: ヘッダ(何が書かれても問題ありません)。

-  6行以降:
   [N] [i_1] [sigma_1] [j_1] [tau_1] ... [i_N] [sigma_N] [j_N] [tau_N] [double01] [double02]

パラメータ
^^^^^^^^^^

-  [ string01 ]

   **形式 :** string型 (空白不可)

   **説明 :** :math:`N` 体ハミルトニアン項の総数のキーワード名を指定します(任意)。

-  [ int01 ]

   **形式 :** int型 (空白不可)

   **説明 :** :math:`N` 体ハミルトニアン項の総数を指定します。

-  [ N ]

   **形式 :** int型 (空白不可)

   **説明 :** 相互作用項の次数を指定します。正の整数である必要があります。

-  [ i_a ], [ j_a ]

   **形式 :** int型 (空白不可)

   **説明 :** サイト番号を指定する整数。0以上 ``Nsite`` 未満で指定します。

-  [ sigma_a ], [ tau_a ]

   **形式 :** int型 (空白不可)

   | **説明 :** スピンを指定する整数。
   | 0: アップスピン
   | 1: ダウンスピン。

-  [ double01 ], [ double02 ]

   **形式 :** double型 (空白不可)

   **説明 :** 相互作用係数の実部と虚部を指定します。

使用ルール
^^^^^^^^^^

本ファイルを使用するにあたってのルールは以下の通りです。

-  行数固定で読み込みを行う為、ヘッダの省略はできません。

-  各成分行は :math:`1 + 4N + 2` 個のフィールドを持つ必要があります。

-  [ int01 ] と定義されている :math:`N` 体ハミルトニアン項の総数が異なる場合はエラー終了します。

-  サイト番号またはスピン番号に範囲外の整数を指定した場合、または係数が非有限値の場合はエラー終了します。

-  スピンを変える因子 :math:`\sigma_a \neq \tau_a` は orbital-general モードでのみ指定できます。
   それ以外のモードでは各因子が :math:`\sigma_a = \tau_a` を満たす必要があります。

-  ``NBodyInterAll`` は通常のエネルギー出力へ寄与し、相互作用専用の出力ファイルは生成しません。

-  real local-energy kernel と :math:`N` 体 Lanczos 補正は未実装です。
   このため ``NBodyInterAll`` を使う場合は複素変分パラメータが必要で、
   ``NLanczosMode`` は 0 とする必要があります。

-  BackFlow では任意の正の入力次数を使用できます。通常の ``Orbital`` /
   ``OrbitalAntiParallel`` 形式では全因子がスピンを保存する必要があります。
   ``OrbitalGeneral`` / FSZ でスピンを変える因子を使用できるのは
   ``2Sz=-1`` の場合だけです。

-  代数的な縮約後、non-FSZ BackFlow は実効次数1をone-body kernelへdispatchし、
   実効次数2以上ではcandidateのBackFlow Slater/Pfaffian状態を完全に再構築します。
   BF-FSZは検証済みの次数1/2 dispatchを使用し、それより高い次数を再構築します。
   したがって真正な実効次数3以上では、残った各項ごとにBackFlow
   Slater/Pfaffianを完全構築する計算量が必要です。

-  native の ``PairHop``、``Exchange``、``InterAll`` 入力は、通常の
   ``Orbital`` / ``OrbitalAntiParallel`` 形式でも ``OrbitalGeneral`` / FSZ でも
   BackFlow 用の二体 Green 関数で評価されます。``NBodyInterAll`` が必要なのは
   3 次以上の項だけです。BackFlowと ``reweight=1`` の併用は
   エラー終了します。

CoulombIntra指定ファイル(coulombintra.def)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

オンサイトクーロン相互作用をハミルトニアンに付け加えます。付け加える項は以下で与えられます。

.. math::

   {\cal H}_U =\sum_{i}U_i n_ {i \uparrow}n_{i \downarrow}

以下にファイル例を記載します。

::

    ======================
    NCoulombIntra 6
    ======================
    ========i_0LocSpn_1IteElc ======
    ======================
       0  4.000000
       1  4.000000
       2  4.000000
       3  4.000000
       4  4.000000
       5  4.000000

ファイル形式
^^^^^^^^^^^^

以下のように行数に応じ異なる形式をとります。

-  1行: ヘッダ(何が書かれても問題ありません)。

-  2行: [string01] [int01]

-  3-5行: ヘッダ(何が書かれても問題ありません)。

-  6行以降: [int02] [double01]

パラメータ
^^^^^^^^^^

-  [ string01 ]

   **形式 :** string型 (空白不可)

   **説明 :**
   オンサイトクーロン相互作用の総数のキーワード名を指定します(任意)。

-  [ int01 ]

   **形式 :** int型 (空白不可)

   **説明 :** オンサイトクーロン相互作用の総数を指定します。

-  [ int02 ]

   **形式 :** int型 (空白不可)

   **説明 :**
   サイト番号を指定する整数。0以上 ``Nsite`` 未満で指定します。

-  [ double01 ]

   **形式 :** double型 (空白不可)

   **説明 :** :math:`U_i` を指定します。

使用ルール
^^^^^^^^^^

本ファイルを使用するにあたってのルールは以下の通りです。

-  行数固定で読み込みを行う為、ヘッダの省略はできません。

-  [ int01 ] と定義されているオンサイトクーロン相互作用の総数が異なる場合はエラー終了します。

-  [ int02 ] を指定する際、範囲外の整数を指定した場合はエラー終了します。

CoulombInter指定ファイル(coulombiter.def)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

サイト間クーロン相互作用をハミルトニアンに付け加えます。付け加える項は以下で与えられます。

.. math::

   {\cal H}_V = \sum_{i,j}V_{ij} n_ {i}n_{j}

以下にファイル例を記載します。

::

    ======================
    NCoulombInter 12
    ======================
    ========CoulombInter ======
    ======================
        0     1         1.500000000000000
        0     3         1.500000000000000
        1     2         1.500000000000000
        1     4         1.500000000000000
        2     0         1.500000000000000
        2     5         1.500000000000000
        3     4         1.500000000000000
        3     0         1.500000000000000
        4     5         1.500000000000000
        4     1         1.500000000000000
        5     3         1.500000000000000
        5     2         1.500000000000000

ファイル形式
^^^^^^^^^^^^

以下のように行数に応じ異なる形式をとります。

-  1行: ヘッダ(何が書かれても問題ありません)。

-  2行: [string01] [int01]

-  3-5行: ヘッダ(何が書かれても問題ありません)。

-  6行以降: [int02] [int03] [double01]

パラメータ
^^^^^^^^^^

-  [ string01 ]

   **形式 :** string型 (空白不可)

   **説明 :**
   サイト間クーロン相互作用の総数のキーワード名を指定します(任意)。

-  [ int01 ]

   **形式 :** int型 (空白不可)

   **説明 :** サイト間クーロン相互作用の総数を指定します。

-  [ int02 ], [ int03 ]

   **形式 :** int型 (空白不可)

   **説明 :**
   サイト番号を指定する整数。0以上 ``Nsite`` 未満で指定します。

-  [ double01 ]

   **形式 :** double型 (空白不可)

   **説明 :** :math:`V_{ij}` を指定します。

使用ルール
^^^^^^^^^^

本ファイルを使用するにあたってのルールは以下の通りです。

-  行数固定で読み込みを行う為、ヘッダの省略はできません。

-  [ int01 ] と定義されているオフサイトクーロン相互作用の総数が異なる場合はエラー終了します。

-  [ int02 ]-[ int03 ] を指定する際、範囲外の整数を指定した場合はエラー終了します。

Hund指定ファイル(hund.def)
~~~~~~~~~~~~~~~~~~~~~~~~~~

Hundカップリングをハミルトニアンに付け加えます。付け加える項は以下で与えられます。

.. math::

   {\cal H}_H =-\sum_{i,j}J_{ij}^{\rm Hund} (n_{i\uparrow}n_{j\uparrow}+n_{i\downarrow}n_{j\downarrow})

以下にファイル例を記載します。

::

    ======================
    NHund 6
    ======================
    ========Hund ======
    ======================
       0     1 -0.250000
       1     2 -0.250000
       2     3 -0.250000
       3     4 -0.250000
       4     5 -0.250000
       5     0 -0.250000

ファイル形式
^^^^^^^^^^^^

以下のように行数に応じ異なる形式をとります。

-  1行: ヘッダ(何が書かれても問題ありません)。

-  2行: [string01] [int01]

-  3-5行: ヘッダ(何が書かれても問題ありません)。

-  6行以降: [int02] [int03] [double01]

パラメータ
^^^^^^^^^^

-  [ string01 ]

   **形式 :** string型 (空白不可)

   **説明 :** Hundカップリングの総数のキーワード名を指定します(任意)。

-  [ int01 ]

   **形式 :** int型 (空白不可)

   **説明 :** Hundカップリングの総数を指定します。

-  [ int02 ], [ int03 ]

   **形式 :** int型 (空白不可)

   **説明 :**
   サイト番号を指定する整数。0以上 ``Nsite`` 未満で指定します。

-  [ double01 ]

   **形式 :** double型 (空白不可)

   **説明 :** :math:`J_{ij}^{\rm Hund}` を指定します。

使用ルール
^^^^^^^^^^

本ファイルを使用するにあたってのルールは以下の通りです。

-  行数固定で読み込みを行う為、ヘッダの省略はできません。

-  [ int01 ] と定義されているHundカップリングの総数が異なる場合はエラー終了します。

-  [ int02 ]-[ int03 ] を指定する際、範囲外の整数を指定した場合はエラー終了します。

PairHop指定ファイル
~~~~~~~~~~~~~~~~~~~

PairHopカップリングをハミルトニアンに付け加えます。付け加える項は以下で与えられます。

.. math::

   {\cal H}_P=\sum_{i,j}J_{ij}^{\rm Pair}
   (c_ {i \uparrow}^{\dagger}c_{j\uparrow}c_{i \downarrow}^{\dagger}c_{j  \downarrow}
   +c_{j \downarrow}^{\dagger}c_{i  \downarrow}c_ {j \uparrow}^{\dagger}c_{i\uparrow})

以下にファイル例を記載します。

::

    ======================
    NPairhop 6
    ======================
    ========Pairhop ======
    ======================
       0     1  0.50000
       1     2  0.50000
       2     3  0.50000
       3     4  0.50000
       4     5  0.50000
       5     0  0.50000

ファイル形式
^^^^^^^^^^^^

以下のように行数に応じ異なる形式をとります。

-  1行: ヘッダ(何が書かれても問題ありません)。

-  2行: [string01] [int01]

-  3-5行: ヘッダ(何が書かれても問題ありません)。

-  6行以降: [int02] [int03] [double01]

パラメータ
^^^^^^^^^^

-  [ string01 ]

   **形式 :** string型 (空白不可)

   **説明 :**
   PairHopカップリングの総数のキーワード名を指定します(任意)。

-  [ int01 ]

   **形式 :** int型 (空白不可)

   **説明 :** PairHopカップリングの総数を指定します。

-  [ int02 ], [ int03 ]

   **形式 :** int型 (空白不可)

   **説明 :**
   サイト番号を指定する整数。0以上 ``Nsite`` 未満で指定します。

-  [ double01 ]

   **形式 :** double型 (空白不可)

   **説明 :** :math:`J_{ij}^{\rm Pair}` を指定します。

使用ルール
^^^^^^^^^^

本ファイルを使用するにあたってのルールは以下の通りです。

-  行数固定で読み込みを行う為、ヘッダの省略はできません。

-  [ int01 ] と定義されているPairHopカップリングの総数が異なる場合はエラー終了します。

-  [ int02 ]-[ int03 ] を指定する際、範囲外の整数を指定した場合はエラー終了します。

Exchange指定ファイル (exchange.def)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Exchangeカップリングをハミルトニアンに付け加えます。 電子系,スピン系の両方の場合に

.. math::

   {\cal H}_E  =\sum_{i,j}J_{ij}^{\rm Ex}
   (c_ {i \uparrow}^{\dagger}c_{j\uparrow}c_{j \downarrow}^{\dagger}c_{i  \downarrow}
   +c_ {i \downarrow}^{\dagger}c_{j\downarrow}c_{j \uparrow}^{\dagger}c_{i  \uparrow})

が付け加えられます(スピン系の場合にHPhiとは定義が異なりますので注意してください)。

::

    ======================
    NExchange 6
    ======================
    ========Exchange ======
    ======================
       0     1  0.50000
       1     2  0.50000
       2     3  0.50000
       3     4  0.50000
       4     5  0.50000
       5     0  0.50000

ファイル形式
^^^^^^^^^^^^

以下のように行数に応じ異なる形式をとります。

-  1行: ヘッダ(何が書かれても問題ありません)。

-  2行: [string01] [int01]

-  3-5行: ヘッダ(何が書かれても問題ありません)。

-  6行以降: [int02] [int03] [double01]

パラメータ
^^^^^^^^^^

-  [ string01 ]

   **形式 :** string型 (空白不可)

   **説明 :**
   Exchangeカップリングの総数のキーワード名を指定します(任意)。

-  [ int01 ]

   **形式 :** int型 (空白不可)

   **説明 :** Exchangeカップリングの総数を指定します。

-  [ int02 ], [ int03 ]

   **形式 :** int型 (空白不可)

   **説明 :**
   サイト番号を指定する整数。0以上 ``Nsite`` 未満で指定します。

-  [ double01 ]

   **形式 :** double型 (空白不可)

   **説明 :** :math:`J_{ij}^{\rm Ex}` を指定します。

使用ルール
^^^^^^^^^^

本ファイルを使用するにあたってのルールは以下の通りです。

-  行数固定で読み込みを行う為、ヘッダの省略はできません。

-  [ int01 ] と定義されているExchangeカップリングの総数が異なる場合はエラー終了します。

-  [ int02 ]-[ int03 ] を指定する際、範囲外の整数を指定した場合はエラー終了します。

Gutzwiller指定ファイル(gutzwiller.def)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Gutzwiller因子

.. math::

   {\cal P}_G=\exp\left[ \sum_i g_i n_{i\uparrow} n_{i\downarrow} \right]

の設定を行います。指定するパラメータはサイト番号 :math:`i` と :math:`g_i` の変分パラメータの番号です。
以下にファイル例を記載します。

::

    ======================
    NGutzwillerIdx 2
    ComplexType 0
    ======================
    ======================
       0     0
       1     0
       2     0
       3     1
    (continue...)
      12     1
      13     0
      14     0
      15     0
       0     1
       1     0

ファイル形式
^^^^^^^^^^^^

以下のように行数に応じ異なる形式をとります( :math:`N_s` はサイト数、 :math:`N_g` は変分パラメータの種類の数)。

-  1行: ヘッダ(何が書かれても問題ありません)。

-  2行: [string01] [int01]

-  3行: [string02] [int02]

-  4-5行: ヘッダ(何が書かれても問題ありません)。

-  6 - (5+ :math:`N_s`)行: [int03] [int04]

-  (6+ :math:`N_s`) - (5+ :math:`N_s` + :math:`N_g`)行：[int05] [int06]

パラメータ
^^^^^^^^^^

-  [ string01 ]

   **形式 :** string型 (空白不可)

   **説明 :**
   :math:`g_i` の変分パラメータの種類の総数のキーワード名を指定します(任意)。

-  [ int01 ]

   **形式 :** int型 (空白不可)

   **説明 :** :math:`g_i` の変分パラメータの種類の総数を指定します。

-  [ string02 ]

   **形式 :** string型 (空白不可)

   **説明 :**
   :math:`g_i` の変分パラメータの型を指定するためのキーワード名を指定します(任意)。

-  [ int02 ]

   **形式 :** int型 (空白不可)

   **説明 :**
   変分パラメータの型を指定する整数。0が実数、1が複素数に対応します。

-  [ int03 ]

   **形式 :** int型 (空白不可)

   **説明 :**
   サイト番号を指定する整数。0以上 ``Nsite`` 未満で指定します。

-  [ int04 ]

   **形式 :** int型 (空白不可)

   **説明 :**
   :math:`g_i` の変分パラメータの種類を表します。0以上[int01]未満で指定します。

-  [ int05 ]

   **形式 :** int型 (空白不可)

   **説明 :**
   :math:`g_i` の変分パラメータの種類を表します(最適化有無の設定用)。0以上[int01]未満で指定します。

-  [ int06 ]

   **形式 :** int型 (空白不可)

   **説明 :**
   [int05]で指定した :math:`g_i` の変分パラメータの最適化有無を設定します。最適化する場合は1、最適化しない場合は0とします。

使用ルール
^^^^^^^^^^

本ファイルを使用するにあたってのルールは以下の通りです。

-  行数固定で読み込みを行う為、ヘッダの省略はできません。

-  [ int01 ] と定義されている変分パラメータの種類の総数が異なる場合はエラー終了します。

-  [ int02 ]-[ int06 ] を指定する際、範囲外の整数を指定した場合はエラー終了します。

Jastrow指定ファイル(jastrow.def)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Jastrow因子

.. math:: {\cal P}_J=\exp\left[\frac{1}{2} \sum_{i\neq j} v_{ij} (n_i-1) (n_j-1)\right]


の設定を行います。指定するパラメータはサイト番号 :math:`i, j` と :math:`v_{ij}` の変分パラメータの番号です。以下にファイル例を記載します。

::

    ======================
    NJastrowIdx 5
    ComplexType 0
    ======================
    ======================
       0     1     0
       0     2     1
       0     3     0
     (continue...)
       0    1
       1    1
       2    1
       3    1
       4    1

ファイル形式
^^^^^^^^^^^^

以下のように行数に応じ異なる形式をとります( :math:`N_s` はサイト数、 :math:`N_j` は変分パラメータの種類の数)。

-  1行: ヘッダ(何が書かれても問題ありません)。

-  2行: [string01] [int01]

-  3行: [string02] [int02]

-  4-5行: ヘッダ(何が書かれても問題ありません)。

-  6 - (5+ :math:`N_s\times (N_s-1)`) 行: [int03] [int04] [int05]

-  (6+ :math:`N_s\times (N_s-1)` ) -
   (5+ :math:`N_s\times (N_s-1)` + :math:`N_j`)行：[int06] [int07]

パラメータ
^^^^^^^^^^

-  [ string01 ]

   **形式 :** string型 (空白不可)

   **説明 :**
   :math:`v_{ij}` の変分パラメータの種類の総数のキーワード名を指定します(任意)。

-  [ int01 ]

   **形式 :** int型 (空白不可)

   **説明 :** :math:`v_{ij}` の変分パラメータの種類の総数を指定します。

-  [ string02 ]

   **形式 :** string型 (空白不可)

   **説明 :**
   :math:`v_{ij}` の変分パラメータの型を指定するためのキーワード名を指定します(任意)。

-  [ int02 ]

   **形式 :** int型 (空白不可)

   **説明 :**
   :math:`v_{ij}` の変分パラメータの型を指定します。0が実数、1が複素数に対応します。

-  [ int03 ], [ int04 ]

   **形式 :** int型 (空白不可)

   **説明 :**
   サイト番号を指定する整数。0以上 ``Nsite`` 未満で指定します。

-  [ int05 ]

   **形式 :** int型 (空白不可)

   **説明 :**
   :math:`v_{ij}` の変分パラメータの種類を表します。0以上[int01]未満で指定します。

-  [ int06 ]

   **形式 :** int型 (空白不可)

   **説明 :**
   :math:`v_{ij}` の変分パラメータの種類を表します(最適化有無の設定用)。0以上[int01]未満で指定します。

-  [ int07 ]

   **形式 :** int型 (空白不可)

   **説明 :**
   [int06]で指定した :math:`v_{ij}` の変分パラメータの最適化有無を設定します。最適化する場合は1、最適化しない場合は0とします。

使用ルール
^^^^^^^^^^

本ファイルを使用するにあたってのルールは以下の通りです。

-  行数固定で読み込みを行う為、ヘッダの省略はできません。

-  [ int01 ] と定義されている変分パラメータの種類の総数が異なる場合はエラー終了します。

-  [ int02 ]-[ int07 ] を指定する際、範囲外の整数を指定した場合はエラー終了します。

Spin Jastrow指定ファイル(spinjastrow.def)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Spin Jastrow因子

.. math::

   {\cal P}_{SJ}=\exp\left[\sum_{i<j} v^s_{ij} m_i m_j\right],
   \quad m_i=n_{i\uparrow}-n_{i\downarrow}

の設定を行います。指定するパラメータはサイト番号 :math:`i, j` と
:math:`v^s_{ij}` の変分パラメータの番号です。以下にファイル例を記載します。

::

    ======================
    NSpinJastrowIdx 4
    ComplexType 0
    =====================
    =====================
       0     1     0
       0     2     1
       0     3     2
       1     0     0
     (continue...)
       0     1
       1     1
       2     1
       3     1

ファイル形式
^^^^^^^^^^^^

以下のように行数に応じ異なる形式をとります( :math:`N_s` はサイト数、 :math:`N_{sj}` は変分パラメータの種類の数)。

-  1行: ヘッダ(何が書かれても問題ありません)。

-  2行: [string01] [int01]

-  3行: [string02] [int02]

-  4-5行: ヘッダ(何が書かれても問題ありません)。

-  6 - (5+ :math:`N_s\times (N_s-1)`) 行: [int03] [int04] [int05]

-  (6+ :math:`N_s\times (N_s-1)` ) -
   (5+ :math:`N_s\times (N_s-1)` + :math:`N_{sj}`)行：[int06] [int07]

パラメータ
^^^^^^^^^^

-  [ string01 ]

   **形式 :** string型 (空白不可)

   **説明 :**
   :math:`v^s_{ij}` の変分パラメータの種類の総数のキーワード名を指定します(任意)。

-  [ int01 ]

   **形式 :** int型 (空白不可)

   **説明 :** :math:`v^s_{ij}` の変分パラメータの種類の総数を指定します。

-  [ string02 ]

   **形式 :** string型 (空白不可)

   **説明 :**
   :math:`v^s_{ij}` の変分パラメータの型を指定するためのキーワード名を指定します(任意)。

-  [ int02 ]

   **形式 :** int型 (空白不可)

   **説明 :**
   :math:`v^s_{ij}` の変分パラメータの型を指定します。0が実数、1が複素数に対応します。

-  [ int03 ], [ int04 ]

   **形式 :** int型 (空白不可)

   **説明 :**
   サイト番号を指定する整数。0以上 ``Nsite`` 未満で指定します。

-  [ int05 ]

   **形式 :** int型 (空白不可)

   **説明 :**
   :math:`v^s_{ij}` の変分パラメータの種類を表します。0以上[int01]未満で指定します。

-  [ int06 ]

   **形式 :** int型 (空白不可)

   **説明 :**
   :math:`v^s_{ij}` の変分パラメータの種類を表します(最適化有無の設定用)。0以上[int01]未満で指定します。

-  [ int07 ]

   **形式 :** int型 (空白不可)

   **説明 :**
   [int06]で指定した :math:`v^s_{ij}` の変分パラメータの最適化有無を設定します。最適化する場合は1、最適化しない場合は0とします。

使用ルール
^^^^^^^^^^

本ファイルを使用するにあたってのルールは以下の通りです。

-  行数固定で読み込みを行う為、ヘッダの省略はできません。

-  [ int01 ] と定義されている変分パラメータの種類の総数が異なる場合はエラー終了します。

-  [ int02 ]-[ int07 ] を指定する際、範囲外の整数を指定した場合はエラー終了します。

-  [ int03 ] と [ int04 ] が同じサイトを指す場合はエラー終了します。

-  各サイトペアについて :math:`(i,j)` と :math:`(j,i)` の両方を指定する必要があります。また、両者の変分パラメータ番号は同一でなければなりません。

BackFlow 指定ファイル (rangebf.def, bf.def)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

BackFlow 相関は、電子配置に依存してペア軌道
:math:`f_{ij}` を有効なペア軌道 :math:`f^b_{ij}` に置き換える相関因子です。
エキスパートモードでは ``namelist.def`` に ``BFRange`` と ``BF`` の
2 つのキーワードを両方指定して有効化します。片方だけを指定した場合は
エラー終了します。

::

    BFRange  rangebf.def
    BF       bf.def

``BFRange`` は各中心サイトから BackFlow が参照する近傍サイトと距離シェルを
指定し、``BF`` はそれらの近傍ペアに対する BackFlow group と最適化フラグを
指定します。現時点の実装はエキスパートモード入力のみを対象としており、
StdFace/Standard mode はこれらのファイルを生成しません。

波動関数の形
^^^^^^^^^^^^

BackFlow 付きペア積波動関数の定義は、Ido ら [Ido2015_] の格子 BackFlow
拡張に従っています。Slater 行列式に対する格子 BackFlow については
Tocchio ら [Tocchio2008_], [Tocchio2011_] も参照してください。
現行実装では BackFlow とスピン射影は併用できないため、BackFlow を使う場合の
波動関数は、BackFlow と併用可能な対角相関因子をまとめて
:math:`{\cal P}_{\rm corr}` と書くと、概念的には

.. math::

   |\Psi^b\rangle
     = {\cal P}_{\rm corr} {\cal L}_{K} |\phi^b\rangle

です。:math:`{\cal L}_{K}` は運動量射影です。運動量射影を用いない場合は
恒等演算子とみなします。実空間配置 :math:`x` で展開すると、

.. math::

   {\cal L}_{K}|\phi^b\rangle
     = \sum_x \left(\frac{N}{2}\right)!
       \sum_{\boldsymbol R} e^{i{\boldsymbol K}\cdot{\boldsymbol R}}
       {\rm Pf}\left[X^b(x_{-{\boldsymbol R}})\right] |x\rangle .

ここで :math:`N` は電子数、:math:`{\boldsymbol R}` は並進ベクトル、
:math:`{\rm Pf}` は Pfaffian です。配置 :math:`x_{\boldsymbol R}` 中の
:math:`n` 番目・ :math:`m` 番目の電子サイトを :math:`i_n, i_m` とすると、
歪対称行列 :math:`X^b` の成分は

.. math::

   X^b_{nm}(x_{\boldsymbol R})
     =
     f^b_{T_{\boldsymbol R}(i_n),T_{\boldsymbol R}(i_m)}
       (x_{\boldsymbol R})
     -
     f^b_{T_{\boldsymbol R}(i_m),T_{\boldsymbol R}(i_n)}
       (x_{\boldsymbol R})

で与えられます。BackFlow は、元のペア軌道 :math:`f_{ij}` を次の
配置依存ペア軌道 :math:`f^b_{ij}(x)` に置き換えます。

.. math::

   f^b_{i_n i_m}(x)
     =
     \sum_{\mu,\nu=0}^{3}\sum_{\tau,\tau'}
     \eta^{\mu\nu}_{\tau\tau'}
     \Theta^{\mu\uparrow}_{i_n,i_n+\tau}(x)
     \Theta^{\nu\downarrow}_{i_m,i_m+\tau'}(x)
     f_{i_n+\tau,i_m+\tau'} .

:math:`\tau,\tau'` は ``BFRange`` で列挙した近傍サイトを表し、
:math:`\eta^{\mu\nu}_{\tau\tau'}` が BackFlow 変分パラメータです。
mVMC の ``ProjBF`` はこの :math:`\eta` を平坦化して格納したものです。
``BF`` ファイルは、各 :math:`(i,j,\tau,\tau')` に対する BackFlow group と
最適化フラグを指定します。現在の実装では group は 0 のみ指定できます。

:math:`\Theta` は配置 :math:`x` 上の局所的な doublon/holon パターンを
表す対角量です。:math:`h_{i\sigma}=1-n_{i\sigma}`,
:math:`D_i=n_{i\uparrow}n_{i\downarrow}`,
:math:`H_i=(1-n_{i\uparrow})(1-n_{i\downarrow})` とすると、

.. math::

   \begin{aligned}
   \Theta^{0\sigma}_{i,i+\tau}(x)
     &= \delta_{i,i+\tau},\\
   \Theta^{1\sigma}_{i,i+\tau}(x)
     &= \langle D_i H_{i+\tau}\rangle_x,\\
   \Theta^{2\sigma}_{i,i+\tau}(x)
     &=
     \langle
       n_{i\sigma}h_{i,-\sigma}
       n_{i+\tau,-\sigma}h_{i+\tau,\sigma}
     \rangle_x,\\
   \Theta^{3\sigma}_{i,i+\tau}(x)
     &=
     \langle
       D_i n_{i+\tau,-\sigma}h_{i+\tau,\sigma}
       + n_{i\sigma}h_{i,-\sigma}H_{i+\tau}
     \rangle_x .
   \end{aligned}

BackFlow 無しの極限は
:math:`\eta^{00}_{0,0}=1` かつその他の :math:`\eta=0` です。
mVMC の入力ではこれが ``ProjBF[0]=1``、``ProjBF[k>0]=0`` に対応します。

パラメータ数
^^^^^^^^^^^^

``rangebf.def`` の ``Nrange`` と ``NzBF``、``bf.def`` の
``NBackFlowIdx`` から、BackFlow の変分パラメータ数は以下のように決まります。

.. math::

   \begin{aligned}
   N_{\rm rangeIdx}
     &= 3\frac{N_{\rm range}-1}{N_z^{\rm BF}} + 1,\\
   N_{\rm BFIdxTotal}
     &= \frac{N_{\rm rangeIdx}(N_{\rm rangeIdx}+1)}{2},\\
   N_{\rm ProjBF}
     &= N_{\rm BFIdxTotal} N_{\rm BackFlowIdx}.
   \end{aligned}

``Nrange`` は各中心サイトに対して列挙するサイト数で、自サイトを含みます。
``NzBF`` は 1 つの非ゼロ距離シェルに含まれるサイト数です。
``Nrange-1`` は ``NzBF`` で割り切れる必要があります。
BackFlow が無い極限は ``ProjBF[0]=1``、``ProjBF[k>0]=0`` です。
すべての ``ProjBF`` を 0 にすることは BackFlow 無しの極限ではありません。
初期値ファイルを使わない場合、mVMC はこの identity 初期値を設定します。

BFRange ファイル (rangebf.def)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

以下に 1 次元 4 サイト周期鎖で、自サイトと左右最近接を BackFlow 範囲にする例を示します。
この例では ``Nrange=3``、``NzBF=2`` です。
``rangebf.def`` は 10 行のヘッダを持ち、本文は 11 行目から始まります。
1-5 行目は mVMC の ``*.def`` ファイル共通ヘッダ、6-10 行目は
BackFlow 用の追加ヘッダです。7 行目の ``Nrange`` 行は読み飛ばされますが、
可読性と既存入力との整合のため、2 行目と同じ内容を書きます。

::

    ====================
    Nrange 3 2
    ====================
    ====================
    ====================
    ====================
    Nrange 3 2
    ====================
    ====================
    ====================
    0 0 0
    0 3 1
    0 1 1
    1 1 0
    1 0 1
    1 2 1
    2 2 0
    2 1 1
    2 3 1
    3 3 0
    3 2 1
    3 0 1

ファイル形式
^^^^^^^^^^^^

-  1-5行: 共通ヘッダ。2行目は ``Nrange`` [int01] [int02] とします。

-  6-10行: BackFlow 用の追加ヘッダ。読み込み時には読み飛ばされますが省略できません。
   7 行目にも ``Nrange`` [int01] [int02] を書くことを推奨します。

-  11 - (10 + :math:`N_s \times` [int01]) 行:
   [int03] [int04] [int05]

パラメータ
^^^^^^^^^^

-  [ int01 ]

   **形式 :** int型 (1以上)

   **説明 :** 各中心サイトから列挙する BackFlow 範囲のサイト数。
   自サイトを必ず含めます。

-  [ int02 ]

   **形式 :** int型 (1以上)

   **説明 :** 1 つの非ゼロ距離シェルに含まれるサイト数。
   [int01] - 1 は [int02] で割り切れる必要があります。

-  [ int03 ]

   **形式 :** int型

   **説明 :** 中心サイト番号。0以上 ``Nsite`` 未満で指定します。

-  [ int04 ]

   **形式 :** int型

   **説明 :** [int03] から見た BackFlow 範囲内のサイト番号。
   0以上 ``Nsite`` 未満で指定します。

-  [ int05 ]

   **形式 :** int型

   **説明 :** 距離シェル番号。自サイトは必ず 0 とし、非ゼロ距離シェルは
   1 以上 :math:`(N_{\rm range}-1)/N_z^{\rm BF}` 以下で指定します。

使用ルール
^^^^^^^^^^

-  各中心サイト [int03] について、ちょうど ``Nrange`` 行を指定する必要があります。

-  各中心サイトには自サイト行 ``i i 0`` が必須です。
   ``i==j`` でシェル番号が 0 以外、または ``i!=j`` でシェル番号が 0 の行は
   エラーになります。

-  同じ中心サイトと範囲サイトの組 ``(i,j)`` を重複して指定することはできません。

-  範囲の関係は相互でなければなりません。自サイト以外の行 ``i j ...`` で
   中心サイト ``i`` の範囲にサイト ``j`` を含めた場合、逆向きの行 ``j i ...`` も
   必要です (両行のシェル番号は一致していなくても構いません)。
   BackFlow の増分更新は、Theta の個数が変化し得る anchor サイトを集めるときに
   この相互性を仮定しているため、非対称な ``BFRange`` は入力時に
   ``BFRange must be mutual`` というエラーで拒否されます。

BF ファイル (bf.def)
^^^^^^^^^^^^^^^^^^^^

``BF`` ファイルでは、``BFRange`` で定義した近傍の組に BackFlow group を
割り当てます。現在の実装では ``NBackFlowIdx==1`` のみをサポートするため、
group 番号 [int06] は常に 0 です。
``bf.def`` も 10 行のヘッダを持ち、本文は 11 行目から始まります。
1-5 行目は mVMC の ``*.def`` ファイル共通ヘッダ、6-10 行目は
BackFlow 用の追加ヘッダです。7 行目の ``NBackFlowIdx`` 行は読み飛ばされますが、
可読性と既存入力との整合のため、2 行目と同じ内容を書きます。

例:

::

    ====================
    NBackFlowIdx 1
    ====================
    ====================
    ====================
    ====================
    NBackFlowIdx 1
    ====================
    ====================
    ====================
    0 0 0 0 0
    0 0 0 3 0
    0 0 0 1 0
    0 0 3 0 0
    0 0 3 3 0
    0 0 3 1 0
    (continue...)
    0 0
    1 0
    2 0
    (continue...)
    9 0

compact 形式:

::

    ====================
    NBackFlowIdx 1
    ====================
    ====================
    ====================
    ====================
    NBackFlowIdx 1
    ====================
    ====================
    ====================
    BFCompact 1
    0 0
    1 0
    2 0
    (continue...)
    9 0

``NBackFlowIdx==1`` の場合は compact 形式を使い、
:math:`N_s^2 N_{\rm range}^2` 行の BackFlow group 行を省略できます。
この形式では、全ての近傍ペアの BackFlow group 番号は暗黙に 0 として扱われます。
上に示した従来の full 形式も引き続き使用できます。

ファイル形式
^^^^^^^^^^^^

-  1-5行: 共通ヘッダ。2行目は ``NBackFlowIdx`` [int01] とします。

-  6-10行: BackFlow 用の追加ヘッダ。読み込み時には読み飛ばされますが省略できません。
   7 行目にも ``NBackFlowIdx`` [int01] を書くことを推奨します。

-  従来の full 形式では、11 - (10 + :math:`N_s^2 N_{\rm range}^2`) 行:
   [int02] [int03] [int04] [int05] [int06]

-  compact 形式では、11 行目に ``BFCompact 1`` を書きます。

-  どちらの形式でも、続く :math:`N_{\rm ProjBF}` 行:
   [int07] [int08]

パラメータ
^^^^^^^^^^

-  [ int01 ]

   **形式 :** int型

   **説明 :** BackFlow group 数。現在の実装では 1 のみ指定できます。

-  [ int02 ], [ int03 ]

   **形式 :** int型

   **説明 :** BackFlow で補正する元のサイト対 :math:`(i,j)`。
   どちらも 0 以上 ``Nsite`` 未満で指定します。

-  [ int04 ], [ int05 ]

   **形式 :** int型

   **説明 :** ``BFRange`` で指定した近傍サイト。
   [int04] は [int02] の BackFlow 範囲に、[int05] は [int03] の
   BackFlow 範囲に含まれている必要があります。

-  [ int06 ]

   **形式 :** int型

   **説明 :** BackFlow group 番号。現在の実装では 0 のみ指定できます。

-  [ int07 ]

   **形式 :** int型

   **説明 :** BackFlow 変分パラメータ ``ProjBF`` の番号。
   0 以上 ``NProjBF`` 未満を重複なくすべて指定します。

-  [ int08 ]

   **形式 :** int型

   **説明 :** [int07] で指定した ``ProjBF`` を最適化する場合は 1、
   固定する場合は 0 とします。``ProjBF[0]`` の虚部は常に固定されます。

使用ルール
^^^^^^^^^^

-  ``BFRange`` と ``BF`` は必ず両方を ``namelist.def`` に指定してください。

-  従来の full 形式では、``BF`` のサイト行は、全ての :math:`i,j` と
   :math:`x_0 \in {\rm BFRange}(i)`、:math:`x_1 \in {\rm BFRange}(j)`
   の組を 1 回ずつ指定する必要があります。行数は
   :math:`N_s^2 N_{\rm range}^2` です。

-  compact 形式では、11 行目に ``BFCompact 1`` を書き、サイト行を省略します。
   この形式は ``NBackFlowIdx==1`` の場合のみ使用できます。省略された
   BackFlow group 番号は 0 として扱われます。

-  ``ProjBF[0]`` は BackFlow の base 項にも使われるため、identity 初期点では
   ``ProjBF[0]=1`` とします。初期値ファイルや再開ファイルを使う場合、
   ``ProjBF[0]`` の虚部は 0 にしてください。

制限事項
^^^^^^^^

BackFlow は現時点では以下の範囲でのみ使用できます。範囲外の入力は
エラー終了します。

-  ``NBackFlowIdx==1`` のみ。複数 BackFlow group は未対応です。

-  ローカル更新は、通常の ``Orbital`` / ``OrbitalAntiParallel``（非 FSZ）では
   ``NExUpdatePath==0``（hopping）または ``NExUpdatePath==1``（hopping + exchange）、
   ``OrbitalGeneral``（FSZ）では ``NExUpdatePath==0`` のみ。
   t-J 用更新経路、Kondo update、doublon-only update、``NExUpdatePath>=2`` の
   pure-spin / weighted local update との併用は未対応です。

-  ペア軌道は ``Orbital`` / ``OrbitalAntiParallel`` の通常形式、または
   ``OrbitalGeneral`` / FSZ を使用できます。``OrbitalGeneral`` では固定
   :math:`S_z` sector と ``2Sz=-1`` の両方を使用できます。``2Sz=-1`` では
   conduction electron の hopping / spin flip sampler と、spin-changing
   ``Transfer`` / ``OneBodyG`` / ``TwoBodyG`` が有効になります。
   局在スピンと BackFlow の併用は未対応です。

-  スピン射影は未対応です。``NSPGaussLeg==1`` を指定してください。

-  RBM、Twist、``NQPOptTrans>1`` は未対応です。周期境界に加えて、負の
   ``NMPTrans`` で指定する反周期境界（``APFlag=1``）を使用できます。
   BF-FSZとBackFlow :math:`N` 体入力を含む全てのBackFlow計算で
   ``reweight=1`` はエラー終了します。
   Single Lanczos Step は ``NVMCCalMode=1`` かつ ``NLanczosMode=1``
   のみ使用でき、Hamiltonian は ``Transfer`` と number-operator 型の
   対角相互作用に限定されます。``OrbitalGeneral`` / FSZ の
   spin-changing ``Transfer`` は ``2Sz=-1`` の場合だけ使用できます。
   ``NLanczosMode>=2``、2nd Lanczos、および ``NBodyG`` /
   ``NBodyInterAll`` に対する BackFlow の Lanczos 補正は未対応です。

-  ``abs(NMPTrans)>1`` による運動量射影は ``OrbitalGeneral`` / FSZ と
   通常の ``Orbital`` / ``OrbitalAntiParallel`` の両形式で、周期・反周期の
   どちらでも使用できます。``abs(NMPTrans)`` は TransSym の ``NQPTrans``
   以下でなければなりません。通常の non-FSZ 形式では ``OptTrans`` を
   rejectし、FSZでは ``NQPOptTrans==1`` の範囲で従来の ``OptTrans`` supportを
   維持します。non-FSZ の ``abs(NMPTrans)>1`` では sampling、Green関数、
   Hamiltonian、1st Lanczos、N体評価を correctness-first の完全
   Slater/Pfaffian再構築で処理します。反周期入力、および先頭変換が非恒等な
   単一pattern入力も同じ再構築を使用し、周期・恒等・単一patternの場合だけ
   従来の増分経路を維持します。再構築経路は従来経路より高コストになる場合があります。

-  ``Orbital`` / ``OrbitalAntiParallel`` の通常形式では、Hamiltonian は
   ``Trans``、number-operator 型の相互作用（``CoulombIntra``,
   ``CoulombInter``, ``Hund``）、および二体項 ``PairHop``、``Exchange``、
   ``InterAll``（スピンを保存する因子対）を使用できます。二体項は
   BackFlow 用の二体 Green 関数で評価され、局所エネルギーのループは
   項について OpenMP 並列化されています。

-  ``OrbitalGeneral`` / FSZ でも、``PairHop``、``Exchange``、``InterAll`` を
   使用できます。測定用の ``OneBodyG``、``TwoBodyG``、``TwoBodyGEx`` は
   general spin label に対応します。``NBodyG`` と ``NBodyInterAll`` で
   general spin labelを使えるのは ``2Sz=-1`` の場合だけで、固定
   :math:`S_z` sectorでは全ての :math:`N` 体因子がスピンを保存する必要があります。
   BackFlow :math:`N` 体入力では複素変分パラメータと
   ``NLanczosMode=0`` が必要です。``NLanczosMode=1`` との併用時は
   ``PairHop``、``Exchange``、``InterAll``、``NBodyG``、
   ``NBodyInterAll`` が対象外です。

-  Standard mode / StdFace から BackFlow 入力は生成されません。
   BackFlow を使う場合は、エキスパートモード入力として ``BFRange`` と ``BF`` を
   手動で用意してください。

BackFlow :math:`N` 体機能のsupport matrix
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. list-table::
   :widths: 28 17 20 35
   :header-rows: 1

   * - BackFlow modeまたは設定
     - ``NBodyG``
     - ``NBodyInterAll``
     - 条件
   * - ``Orbital`` / ``OrbitalAntiParallel``
     - 対応
     - 対応
     - 複素変分パラメータ、spin-conserving因子、
       ``NLanczosMode=0``。
   * - ``OrbitalGeneral`` / FSZ、固定 :math:`S_z`
     - 対応
     - 対応
     - 複素変分パラメータ、spin-conserving因子、
       ``NLanczosMode=0``。
   * - ``OrbitalGeneral`` / FSZ、``2Sz=-1``
     - 対応
     - 対応
     - 複素変分パラメータ、spin-changing因子を使用可能、
       ``NLanczosMode=0``。
   * - 実変分パラメータ
     - エラー終了
     - エラー終了
     - real BackFlow :math:`N` 体measurement/local-energy kernelは未実装。
   * - ``reweight=1``
     - エラー終了
     - エラー終了
     - reweightはBackFlowのsupport contract外。
   * - ``NLanczosMode>0``
     - エラー終了
     - エラー終了
     - :math:`N` 体Lanczos補正は未実装。

DH2指定ファイル
~~~~~~~~~~~~~~~

.. math::

   {\cal P}_{d-h}^{(2)}= \exp \left[ \sum_t \sum_{n=0}^2
   (\alpha_{2nt}^d \sum_{i}\xi_{i2nt}^d+\alpha_{2nt}^h \sum_{i}\xi_{i2nt}^h)\right]

で表される2サイトのdoublon-holon相関因子の設定を行います。
指定するパラメータはサイト番号 :math:`i` とその周囲2サイト、
:math:`\alpha_{2nt}^{d,h}` の変分パラメータの番号で、
変分パラメータは各サイト毎に :math:`t` 種類設定します。
各パラメータ、演算子に関する詳細は文献 [Tahara2008_ ]をご覧ください。
以下にファイル例を記載します。

::

    ====================================
    NDoublonHolon2siteIdx 2
    ComplexType 0
    ====================================
    ====================================
       0     5   15    0
       0    13    7    1
       1     6   12    0
       1    14    4    1
     (continue...)
      15     0   10    0
      15     8    2    1
       0     1
       1     1
       2     1
    (continue...)
      10     1
      11     1

ファイル形式
^^^^^^^^^^^^

以下のように行数に応じ異なる形式をとります( :math:`N_s` はサイト数、 :math:`N_{\rm dh2}` は変分パラメータの種類の数)。

-  1行: ヘッダ(何が書かれても問題ありません)。

-  2行: [string01] [int01]

-  3行: [string02] [int02]

-  4-5行: ヘッダ(何が書かれても問題ありません)。

-  6 - (5+ :math:`N_s\times N_{\rm dh2}`)行:
   [int03] [int04] [int05] [int06]

-  (6+ :math:`N_s\times N_{\rm dh2}`) -
   (5+ :math:`(N_s+6) \times N_{\rm dh2}`)行：[int07] [int08]

パラメータ
^^^^^^^^^^

-  [ string01 ]

   **形式 :** string型 (空白不可)

   **説明 :**
   変分パラメータのセット総数のキーワード名を指定します(任意)。

-  [ int01 ]

   **形式 :** int型 (空白不可)

   **説明 :** 変分パラメータのセット総数を指定します。

-  [ string02 ]

   **形式 :** string型 (空白不可)

   **説明
   :** 変分パラメータの型を指定するためのキーワード名を指定します(任意)。

-  [ int02 ]

   **形式 :** int型 (空白不可)

   **説明 :**
   変分パラメータの型を指定します。0が実数、1が複素数に対応します。

-  [ int03 ], [ int04 ],
   [ int05 ]

   **形式 :** int型 (空白不可)

   **説明 :**
   サイト番号を指定する整数。0以上 ``Nsite`` 未満で指定します。

-  [ int06 ]

   **形式 :** int型 (空白不可)

   **説明 :**
   変分パラメータの種類を表します。0以上[int01]未満で指定します。

-  [ int07 ]

   **形式 :** int型 (空白不可)

   **説明 :** 変分パラメータの種類を表します(最適化有無の設定用)。値は

   -  | :math:`n`: 周囲のdoublon(holon)数 (0, 1, 2)

   -  | :math:`s`: 中心がdoublonの場合 :math:`0`,
        中心がholonの場合 :math:`1`

   -  :math:`t`: 変分パラメータのセット番号(0, :math:`\cdots` [int1]-1)

   として、 :math:`(2n+s)\times` [int01] :math:`+t` を設定します。

-  [ int08 ]

   **形式 :** int型 (空白不可)

   **説明 :**
   [int07]で指定した変分パラメータの最適化有無を設定します。最適化する場合は1、最適化しない場合は0とします。

使用ルール
^^^^^^^^^^

本ファイルを使用するにあたってのルールは以下の通りです。

-  行数固定で読み込みを行う為、ヘッダの省略はできません。

-  [ int01 ] と定義されている変分パラメータの種類の総数が異なる場合はエラー終了します。

-  [ int02 ]-[ int08 ] を指定する際、範囲外の整数を指定した場合はエラー終了します。

DH4指定ファイル
~~~~~~~~~~~~~~~

.. math::

   {\cal P}_{d-h}^{(4)}= \exp \left[ \sum_t \sum_{n=0}^4
   (\alpha_{4nt}^d \sum_{i}\xi_{i4nt}^d+\alpha_{4nt}^h \sum_{i}\xi_{i4nt}^h)\right]

で表される4サイトのdoublon-holon相関因子の設定を行います。
指定するパラメータはサイト番号 :math:`i` とその周囲4サイト、
:math:`\alpha_{4nt}^{d,h}` の変分パラメータの番号で、
変分パラメータは各サイト毎に :math:`t` 種類設定します。
各パラメータ、演算子に関する詳細は文献[Tahara2008_ ]をご覧ください。
以下にファイル例を記載します。

::

    ====================================
    NDoublonHolon4siteIdx 1
    ComplexType 0
    ====================================
    ====================================
       0     1    3    4   12    0
       1     2    0    5   13    0
       2     3    1    6   14    0
       3     0    2    7   15    0
     (continue...)
      14    15   13    2   10    0
      15    12   14    3   11    0
       0     1
       1     1
    (continue...)
       8     1
       9     1

ファイル形式
^^^^^^^^^^^^

以下のように行数に応じ異なる形式をとります( :math:`N_s` はサイト数、 :math:`N_{\rm dh4}` は変分パラメータの種類の数)。

-  1行: ヘッダ(何が書かれても問題ありません)。

-  2行: [string01] [int01]

-  3行: [string02] [int02]

-  4-5行: ヘッダ(何が書かれても問題ありません)。

-  6 - (5+ :math:`N_s\times N_{\rm dh4}`)行:
   [int03] [int04] [int05] [int06] [int07] [int08]

-  (6+ :math:`N_s\times N_{\rm dh4}`) -
   (5+ :math:`(N_s+10) \times N_{\rm dh4}`)行：[int09] [int10]

パラメータ
^^^^^^^^^^

-  [ string01 ]

   **形式 :** string型 (空白不可)

   **説明 :**
   変分パラメータのセット総数のキーワード名を指定します(任意)。

-  [ int01 ]

   **形式 :** int型 (空白不可)

   **説明 :** 変分パラメータのセット総数を指定します。

-  [ string02 ]

   **形式 :** string型 (空白不可)

   **説明 :**
   変分パラメータの型を指定するためのキーワード名を指定します(任意)。

-  [ int02 ]

   **形式 :** int型 (空白不可)

   **説明 :**
   変分パラメータの型を指定します。0が実数、1が複素数に対応します。

-  [ int03 ], [ int04 ],
   [ int05 ], [ int06 ],
   [ int07 ]

   **形式 :** int型 (空白不可)

   **説明 :**
   サイト番号を指定する整数。0以上 ``Nsite`` 未満で指定します。

-  [ int08 ]

   **形式 :** int型 (空白不可)

   **説明 :**
   変分パラメータの種類を表します。0以上[int01]未満で指定します。

-  [ int09 ]

   **形式 :** int型 (空白不可)

   **説明 :** 変分パラメータの種類を表します(最適化有無の設定用)。値は

   -  | :math:`n`: 周囲のdoublon(holon)数 (0, 1, 2, 3, 4)

   -  | :math:`s`: 中心がdoublonの場合 :math:`0`,
        中心がholonの場合 :math:`1`

   -  :math:`t`: 変分パラメータのセット番号(0, :math:`\cdots` [int1]-1)

   として、 :math:`(2n+s)\times` [int01] :math:`[2]+t` を設定します。

-  [ int10 ]

   **形式 :** int型 (空白不可)

   **説明 :**
   [int09]で指定した変分パラメータの最適化有無を設定します。最適化する場合は1、最適化しない場合は0とします。

使用ルール
^^^^^^^^^^

本ファイルを使用するにあたってのルールは以下の通りです。

-  行数固定で読み込みを行う為、ヘッダの省略はできません。

-  [ int01 ] と定義されている変分パラメータの種類の総数が異なる場合はエラー終了します。

-  [ int02 ]-[ int10 ] を指定する際、範囲外の整数を指定した場合はエラー終了します。

GeneralRBM_PhysLayer指定ファイル
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

RBM因子

.. math::

   {\cal N}_{\rm General RBM}= \exp \left[ \sum_i a_{i\sigma} n_{i\sigma} \right] \prod_k^{N_{\rm neuronGeneral}} \cosh \left[ b_k + \sum_{i\sigma} W_{i\sigma k} n_{i\sigma} \right]

のうち、 :math:`\exp \left[ \sum_i a_{i\sigma} n_{i\sigma} \right] ` の設定を行います。指定するパラメータはサイト番号 :math:`i`、スピン番号 :math:`\sigma`、および :math:`a_{i \sigma}` の変分パラメータの番号です。以下にファイル例を記載します。

::

    --------------------
    NRBM_PhysLayerIdx	1
    ComplexType	1
    i s RBM_PhysLayer_Idx
    --------------------
     0	0	 0
     0	1	 0
     1	0	 0
     1	1	 0
     (continue...)
     0    1

ファイル形式
^^^^^^^^^^^^

以下のように行数に応じ異なる形式をとります( :math:`N_s` はサイト数、 :math:`N_v` は変分パラメータの種類の数)。

-  1行: ヘッダ(何が書かれても問題ありません)。

-  2行: [string01] [int01]

-  3行: [string02] [int02]

-  4-5行: ヘッダ(何が書かれても問題ありません)。

-  6 - (5+ :math:`2N_s`) 行: [int03] [int04] [int05]

-  (6+ :math:`2N_s` ) -
   (5+ :math:`2N_s` + :math:`N_v`)行：[int06] [int07]

パラメータ
^^^^^^^^^^

-  [ string01 ]

   **形式 :** string型 (空白不可)

   **説明 :**
   :math:`a_{i \sigma }` の変分パラメータの種類の総数のキーワード名を指定します(任意)。

-  [ int01 ]

   **形式 :** int型 (空白不可)

   **説明 :** :math:`a_{i \sigma }` の変分パラメータの種類の総数を指定します。

-  [ string02 ]

   **形式 :** string型 (空白不可)

   **説明 :**
   :math:`a_{i \sigma }` の変分パラメータの型を指定するためのキーワード名を指定します(任意)。

-  [ int02 ]

   **形式 :** int型 (空白不可)

   **説明 :**
   :math:`a_{i \sigma }` の変分パラメータの型を指定します。0が実数、1が複素数に対応します。

-  [ int03 ]

   **形式 :** int型 (空白不可)

   **説明 :**
   サイト番号を指定する整数。0以上 ``Nsite`` 未満で指定します。

-  [ int04 ]

   **形式 :** int型 (空白不可)

   **説明 :**
   スピン番号を指定する整数。0もしくは1で指定します。

-  [ int05 ]

   **形式 :** int型 (空白不可)

   **説明 :**
   :math:`a_{i \sigma }` の変分パラメータの種類を表します。0以上[int01]未満で指定します。

-  [ int06 ]

   **形式 :** int型 (空白不可)

   **説明 :**
   :math:`a_{i \sigma }` の変分パラメータの種類を表します(最適化有無の設定用)。0以上[int01]未満で指定します。

-  [ int07 ]

   **形式 :** int型 (空白不可)

   **説明 :**
   [int06]で指定した :math:`a_{i \sigma }` の変分パラメータの最適化有無を設定します。最適化する場合は1、最適化しない場合は0とします。

使用ルール
^^^^^^^^^^

本ファイルを使用するにあたってのルールは以下の通りです。

-  本機能はベータ版のため、使用には十分注意してください。また、正式リリースした際に、ファイル形式や実装が変更される可能性があります。

-  ``ComplexType=0`` と ``ComplexType=1`` の両方を指定できます。RBM層を含む有効な変分パラメータ定義がすべて ``ComplexType=0`` の場合、mVMCは実数VMC経路を使用します。いずれか1つでも ``ComplexType=1`` の定義があれば複素数経路を使用します。実数RBMは ``Orbital`` （``OrbitalGeneral`` / FSZは不可）、``NExUpdatePath=0`` または ``1``、``NLanczosMode=0`` の組合せに対応します。BackFlowおよびPower Lanczosとの併用には対応していません。同じ経路選択と制約がChargeRBMおよびSpinRBMの各層定義にも適用されます。

-  行数固定で読み込みを行う為、ヘッダの省略はできません。

-  [ int01 ] と定義されている変分パラメータの種類の総数が異なる場合はエラー終了します。

-  [ int02 ]-[ int07 ] を指定する際、範囲外の整数を指定した場合はエラー終了します。


GeneralRBM_HiddenLayer指定ファイル
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

RBM因子

.. math::

   {\cal N}_{\rm General RBM}= \exp \left[ \sum_i a_{i\sigma} n_{i\sigma} \right] \prod_k^{N_{\rm neuronGeneral}} \cosh \left[ b_k + \sum_{i\sigma} W_{i\sigma k} n_{i\sigma} \right]


のうち、 :math:`\prod_k^{N_{\rm neuronGeneral}} \cosh \left[ b_k + \sum_{i\sigma} W_{i\sigma k} n_{i\sigma} \right]` の隠れ層のみが関わる箇所の設定を行います。指定するパラメータは隠れニューロンの番号 :math:`k` と :math:`h_{k}` の変分パラメータの番号です。以下にファイル例を記載します。

::

    --------------------
    NRBM_HiddenLayerIdx	2
    ComplexType	1
    k RBM_HiddenLayer_Idx
    --------------------
     0	0
     1	0
     2	0
     3	0
     (continue...)
     0    1
     1    1

ファイル形式
^^^^^^^^^^^^

以下のように行数に応じ異なる形式をとります( :math:`N_{\rm neuronGeneral}` は隠れニューロンの数、 :math:`N_v` は変分パラメータの種類の数)。

-  1行: ヘッダ(何が書かれても問題ありません)。

-  2行: [string01] [int01]

-  3行: [string02] [int02]

-  4-5行: ヘッダ(何が書かれても問題ありません)。

-  6 - (5+ :math:`N_{\rm neuronGeneral}`) 行: [int03] [int04]

-  (6+ :math:`N_{\rm neuronGeneral}` ) -
   (5+ :math:`N_{\rm neuronGeneral}` + :math:`N_v`)行：[int05] [int06]

パラメータ
^^^^^^^^^^

-  [ string01 ]

   **形式 :** string型 (空白不可)

   **説明 :**
   :math:`h_{k}` の変分パラメータの種類の総数のキーワード名を指定します(任意)。

-  [ int01 ]

   **形式 :** int型 (空白不可)

   **説明 :** :math:`h_{k}` の変分パラメータの種類の総数を指定します。

-  [ string02 ]

   **形式 :** string型 (空白不可)

   **説明 :**
   :math:`h_{k}` の変分パラメータの型を指定するためのキーワード名を指定します(任意)。

-  [ int02 ]

   **形式 :** int型 (空白不可)

   **説明 :**
   :math:`h_{k}` の変分パラメータの型を指定します。0が実数、1が複素数に対応します。

-  [ int03 ]

   **形式 :** int型 (空白不可)

   **説明 :**
   サイト番号を指定する整数。0以上 ``NneuronGeneral`` 未満で指定します。

-  [ int04 ]

   **形式 :** int型 (空白不可)

   **説明 :**
   :math:`h_{k}` の変分パラメータの種類を表します。0以上[int01]未満で指定します。

-  [ int05 ]

   **形式 :** int型 (空白不可)

   **説明 :**
   :math:`h_{k}` の変分パラメータの種類を表します(最適化有無の設定用)。0以上[int01]未満で指定します。

-  [ int06 ]

   **形式 :** int型 (空白不可)

   **説明 :**
   [int06]で指定した :math:`h_{k}` の変分パラメータの最適化有無を設定します。最適化する場合は1、最適化しない場合は0とします。

使用ルール
^^^^^^^^^^

本ファイルを使用するにあたってのルールは以下の通りです。

-  本機能はベータ版のため、使用には十分注意してください。また、正式リリースした際に、ファイル形式や実装が変更される可能性があります。

-  ``ComplexType=0`` と ``ComplexType=1`` の両方を指定できます。実数RBMの対応範囲と制約は上記のGeneralRBM_PhysLayer指定ファイルの節を参照してください。

-  行数固定で読み込みを行う為、ヘッダの省略はできません。

-  [ int01 ] と定義されている変分パラメータの種類の総数が異なる場合はエラー終了します。

-  [ int02 ]-[ int06 ] を指定する際、範囲外の整数を指定した場合はエラー終了します。

GeneralRBM_PhysHidden指定ファイル
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

RBM因子

.. math::

   {\cal N}_{\rm General RBM}= \exp \left[ \sum_i a_{i\sigma} n_{i\sigma} \right] \prod_k^{N_{\rm neuronGeneral}} \cosh \left[ b_k + \sum_{i\sigma} W_{i\sigma k} n_{i\sigma} \right]


のうち、 :math:`\prod_k^{N_{\rm neuronGeneral}} \cosh \left[ b_k + \sum_{i\sigma} W_{i\sigma k} n_{i\sigma} \right]` の隠れ層と物理層どちらも関わる箇所の設定を行います。 指定するパラメータはサイト番号、スピン番号、および隠れニューロンの番号 :math:`i \sigma k` と :math:`W_{i \sigma k}` の変分パラメータの番号です。以下にファイル例を記載します。

::

    --------------------
    NRBM_HiddenLayerIdx	32
    ComplexType	1
    i s k RBM_PhysHidden_Idx
    --------------------
     0	0   0   0
     0	1   0   1
     1	0   0   2
     1	1   0   3
     2	0   0   4
     2	1   0   5
     (continue...)
     0    1
     1    1
     (continue...)

ファイル形式
^^^^^^^^^^^^

以下のように行数に応じ異なる形式をとります( :math:`N_s` はサイト数、:math:`N_{\rm neuronGeneral}` は隠れニューロンの数、 :math:`N_v` は変分パラメータの種類の数)。

-  1行: ヘッダ(何が書かれても問題ありません)。

-  2行: [string01] [int01]

-  3行: [string02] [int02]

-  4-5行: ヘッダ(何が書かれても問題ありません)。

-  6 - (5+ :math:`2 N_s N_{\rm neuronGeneral}`) 行: [int03] [int04] [int05] [int06]

-  (6+ :math:`2 N_s N_{\rm neuronGeneral}` ) -
   (5+ :math:`2 N_s N_{\rm neuronGeneral}` + :math:`N_v`)行：[int07] [int08]

パラメータ
^^^^^^^^^^

-  [ string01 ]

   **形式 :** string型 (空白不可)

   **説明 :**
   :math:`W_{i \sigma k}` の変分パラメータの種類の総数のキーワード名を指定します(任意)。

-  [ int01 ]

   **形式 :** int型 (空白不可)

   **説明 :** :math:`W_{i \sigma k}` の変分パラメータの種類の総数を指定します。

-  [ string02 ]

   **形式 :** string型 (空白不可)

   **説明 :**
   :math:`W_{i \sigma k}` の変分パラメータの型を指定するためのキーワード名を指定します(任意)。

-  [ int02 ]

   **形式 :** int型 (空白不可)

   **説明 :**
   :math:`W_{i \sigma k}` の変分パラメータの型を指定します。0が実数、1が複素数に対応します。

-  [ int03 ]

   **形式 :** int型 (空白不可)

   **説明 :**
   サイト番号を指定する整数。0以上 ``Nsite`` 未満で指定します。

-  [ int04 ]

   **形式 :** int型 (空白不可)

   **説明 :**
   スピン番号を指定する整数。0もしくは1で指定します。

-  [ int05 ]

   **形式 :** int型 (空白不可)

   **説明 :**
   サイト番号を指定する整数。0以上 ``NneuronGeneral`` 未満で指定します。

-  [ int06 ]

   **形式 :** int型 (空白不可)

   **説明 :**
   :math:`W_{i \sigma k}` の変分パラメータの種類を表します。0以上[int01]未満で指定します。

-  [ int07 ]

   **形式 :** int型 (空白不可)

   **説明 :**
   :math:`W_{i \sigma k}` の変分パラメータの種類を表します(最適化有無の設定用)。0以上[int01]未満で指定します。

-  [ int08 ]

   **形式 :** int型 (空白不可)

   **説明 :**
   [int06]で指定した :math:`W_{i \sigma k}` の変分パラメータの最適化有無を設定します。最適化する場合は1、最適化しない場合は0とします。

使用ルール
^^^^^^^^^^

本ファイルを使用するにあたってのルールは以下の通りです。

-  本機能はベータ版のため、使用には十分注意してください。また、正式リリースした際に、ファイル形式や実装が変更される可能性があります。

-  ``ComplexType=0`` と ``ComplexType=1`` の両方を指定できます。実数RBMの対応範囲と制約は上記のGeneralRBM_PhysLayer指定ファイルの節を参照してください。

-  行数固定で読み込みを行う為、ヘッダの省略はできません。

-  [ int01 ] と定義されている変分パラメータの種類の総数が異なる場合はエラー終了します。

-  [ int02 ]-[ int08 ] を指定する際、範囲外の整数を指定した場合はエラー終了します。

Orbital/OrbitalAntiParallel指定ファイル(orbitalidx.def)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. math:: |\phi_{\rm pair} \rangle = \left[\sum_{i, j=1}^{N_s} f_{ij}c_{i\uparrow}^{\dagger}c_{j\downarrow}^{\dagger} \right]^{N/2}|0 \rangle

で表されるペア軌道の設定を行います。指定するパラメータはサイト番号 :math:`i, j` と変分パラメータの種類を設定します。以下にファイル例を記載します。

::

    ====================================
    NOrbitalIdx 64
    ComplexType 0
    ====================================
    ====================================
       0     0     0
       0     1     1
       0     2     2
       0     3     3
     (continue...)
      15     9    62
      15    10    63
       0    1
       1    1
    (continue...)
      62    1
      63    1

ファイル形式
^^^^^^^^^^^^

以下のように行数に応じ異なる形式をとります( :math:`N_s` はサイト数、 :math:`N_{\rm o}` は変分パラメータの種類の数)。

-  1行: ヘッダ(何が書かれても問題ありません)。

-  2行: [string01] [int01]

-  3行: [string02] [int02]

-  4-5行: ヘッダ(何が書かれても問題ありません)。

-  6 - (5+ :math:`N_s^2`)行: [int03] [int04] [int05] [int06]

-  (6+ :math:`N_s^2`) - (5+ :math:`N_s^2+N_{\rm o}`)行：[int07] [int08]

パラメータ
^^^^^^^^^^

-  [ string01 ]

   **形式 :** string型 (空白不可)

   **説明 :**
   変分パラメータのセット総数のキーワード名を指定します(任意)。

-  [ int01 ]

   **形式 :** int型 (空白不可)

   **説明 :** 変分パラメータのセット総数を指定します。

-  [ string02 ]

   **形式 :** string型 (空白不可)

   **説明 :**
   変分パラメータの型を指定するためのキーワード名を指定します(任意)。

-  [ int02 ]

   **形式 :** int型 (空白不可)

   **説明 :**
   変分パラメータの型を指定します。0が実数、1が複素数に対応します。

-  [ int03 ], [ int04 ]

   **形式 :** int型 (空白不可)

   **説明 :**
   サイト番号を指定する整数。0以上 ``Nsite`` 未満で指定します。

-  [ int05 ]

   **形式 :** int型 (空白不可)

   **説明 :**
   変分パラメータの種類を表します。0以上[int01]未満で指定します。

-  [ int06 ]

   **形式 :** int型

   **説明 :**
   反周期境界条件モードがON( ``ModPara`` ファイルで ``NMPTrans`` が負の場合に有効)の場合、変分パラメータ :math:`f_{ij}` の番号の他に符号を反転するか否かを直接指定する。 [ int06 ] = :math:`\pm1` により符号を指定する。反周期境界条件モードがOFFの場合は省略可能。

-  [ int07 ]

   **形式 :** int型 (空白不可)

   **説明 :**
   変分パラメータの種類を表します(最適化有無の設定用)。0以上[int01]未満で指定します。

-  [ int08 ]

   **形式 :** int型 (空白不可)

   **説明 :**
   [int07]で指定した変分パラメータの最適化有無を設定します。最適化する場合は1、最適化しない場合は0とします。

使用ルール
^^^^^^^^^^

本ファイルを使用するにあたってのルールは以下の通りです。

-  行数固定で読み込みを行う為、ヘッダの省略はできません。

-  [ int01 ] と定義されている変分パラメータの種類の総数が異なる場合はエラー終了します。

-  [ int02 ]-[ int08 ] を指定する際、範囲外の整数を指定した場合はエラー終了します。

OrbitalParallel指定ファイル
~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. math:: |\phi_{\rm pair} \rangle = \left[\sum_{i, j=1}^{N_s} f_{i\sigma j\sigma}c_{i\sigma}^{\dagger}c_{j\sigma}^{\dagger} \right]^{N/2}|0 \rangle

で表されるペア軌道の設定を行います。指定するパラメータはサイト番号 :math:`i, j` と変分パラメータの種類を設定します。以下にファイル例を記載します。

::

    ====================================
    NOrbitalIdx 120
    ComplexType 0
    ====================================
    ====================================
       0     1     0
       0     2     1
       0     3     2
     (continue...)
      15    13    118
      15    14    119
       0    1
       1    1
    (continue...)
      118    1
      119    1

ファイル形式
^^^^^^^^^^^^

以下のように行数に応じ異なる形式をとります( :math:`N_s` はサイト数、 :math:`N_{\rm o}` は変分パラメータの種類の数)。

-  1行: ヘッダ(何が書かれても問題ありません)。

-  2行: [string01] [int01]

-  3行: [string02] [int02]

-  4-5行: ヘッダ(何が書かれても問題ありません)。

-  6 - (5+ :math:`N_s*(N_s-1)/2`)行: [int03] [int04] [int05] [int06]

-  (6+ :math:`N_s*(N_s-1)/2`) -
   (5+ :math:`N_s*(N_s-1)/2+N_{\rm o}`)行：[int07] [int08]

パラメータ
^^^^^^^^^^

-  [ string01 ]

   **形式 :** string型 (空白不可)

   **説明 :**
   変分パラメータのセット総数のキーワード名を指定します(任意)。

-  [ int01 ]

   **形式 :** int型 (空白不可)

   **説明 :** 変分パラメータのセット総数を指定します。

-  [ string02 ]

   **形式 :** string型 (空白不可)

   **説明 :**
   変分パラメータの型を指定するためのキーワード名を指定します(任意)。

-  [ int02 ]

   **形式 :** int型 (空白不可)

   **説明 :**
   変分パラメータの型を指定します。0が実数、1が複素数に対応します。

-  [ int03 ], [ int04 ]

   **形式 :** int型 (空白不可)

   **説明 :**
   サイト番号を指定する整数。0以上 ``Nsite`` 未満で指定します。

-  [ int05 ]

   **形式 :** int型 (空白不可)

   **説明 :**
   変分パラメータの種類を表します。0以上[int01]未満で指定します。

-  [ int06 ]

   **形式 :** int型

   **説明 :**
   反周期境界条件モードがON( ``ModPara`` ファイルで ``NMPTrans`` が負の場合に有効)の場合、変分パラメータ :math:`f_{ij}` の番号の他に符号を反転するか否かを直接指定する。 [ int06 ] = :math:`\pm1` により符号を指定する。反周期境界条件モードがOFFの場合は省略可能。

-  [ int07 ]

   **形式 :** int型 (空白不可)

   **説明 :**
   変分パラメータの種類を表します(最適化有無の設定用)。0以上[int01]未満で指定します。

-  [ int08 ]

   **形式 :** int型 (空白不可)

   **説明 :**
   [int06]で指定した変分パラメータの最適化有無を設定します。最適化する場合は1、最適化しない場合は0とします。

使用ルール
^^^^^^^^^^

本ファイルを使用するにあたってのルールは以下の通りです。

-  行数固定で読み込みを行う為、ヘッダの省略はできません。

-  [ int01 ] と定義されている変分パラメータの種類の総数が異なる場合はエラー終了します。

-  [ int02 ]-[ int08 ] を指定する際、範囲外の整数を指定した場合はエラー終了します。

OrbitalGeneral指定ファイル
~~~~~~~~~~~~~~~~~~~~~~~~~~

.. math:: |\phi_{\rm pair} \rangle = \left[\sum_{i, j=1}^{N_s} \sum_{\sigma_1, \sigma_2}f_{i\sigma_1 j \sigma_2}c_{i\sigma_1}^{\dagger}c_{j\sigma_2}^{\dagger} \right]^{N/2}|0 \rangle

で表されるペア軌道の設定を行います。指定するパラメータはサイト番号 :math:`i, j`,
スピン :math:`\sigma_1, \sigma_2` と変分パラメータの種類を設定します。 :math:`i+\sigma_1 N_s < j+\sigma_2 N_s` (:math:`\sigma=0 ,1`)を満たすように指定する必要があります。以下にファイル例を記載します。

::

    ====================================
    NOrbitalIdx 255
    ComplexType 0
    ====================================
    ====================================
       0  0  0  1  0
       0  0  1  1  1
     (continue...)
      14  0  15 1 253
      15  0  15  1  254
       0    1
       1    1
    (continue...)
      253   1
      254    1

ファイル形式
^^^^^^^^^^^^

以下のように行数に応じ異なる形式をとります( :math:`N_s` はサイト数、 :math:`N_{\rm o}` は変分パラメータの種類の数)。
変分パラメータの総数 :math:`N_p` は :math:`i+\sigma_1 N_s < j+\sigma_2 N_s` (:math:`\sigma=0 ,1`)を満たすペアの総数に対応し、
:math:`S_z=0` の場合は :math:`N_p=N_s^2` 、 :math:`S_z` 非保存の場合は :math:`N_p=2N_s^2-N_s` 個となります。

-  1行: ヘッダ(何が書かれても問題ありません)。

-  2行: [string01] [int01]

-  3行: [string02] [int02]

-  4-5行: ヘッダ(何が書かれても問題ありません)。

-  6 - (5+ :math:`N_p`)行:
   [int03] [int04] [int05] [int06] [int07] [int08]

-  (6+ :math:`N_p`) - (5+ :math:`N_p+N_{\rm o}`)行：[int09] [int10]

パラメータ
^^^^^^^^^^

-  [ string01 ]

   **形式 :** string型 (空白不可)

   **説明 :**
   変分パラメータのセット総数のキーワード名を指定します(任意)。

-  [ int01 ]

   **形式 :** int型 (空白不可)

   **説明 :** 変分パラメータのセット総数を指定します。

-  [ string02 ]

   **形式 :** string型 (空白不可)

   **説明 :**
   変分パラメータの型を指定するためのキーワード名を指定します(任意)。

-  [ int02 ]

   **形式 :** int型 (空白不可)

   **説明 :**
   変分パラメータの型を指定します。0が実数、1が複素数に対応します。

-  [ int03 ], [ int05 ]

   **形式 :** int型 (空白不可)

   **説明 :**
   サイト番号を指定する整数。0以上 ``Nsite`` 未満で指定します。

-  [ int04 ], [ int06 ]

   **形式 :** int型 (空白不可)

   **説明 :** スピンを指定する整数。0が :math:`\uparrow` スピン,
   1が :math:`\downarrow` スピンに対応します。

-  [ int07 ]

   **形式 :** int型 (空白不可)

   **説明 :**
   変分パラメータの種類を表します。0以上[int01]未満で指定します。

-  [ int08 ]

   **形式 :** int型

   **説明 :**
   反周期境界条件モードがON( ``ModPara`` ファイルで ``NMPTrans`` が負の場合に有効)の場合、変分パラメータ :math:`f_{ij}` の番号の他に符号を反転するか否かを直接指定する。 [ int08 ] = :math:`\pm1` により符号を指定する。反周期境界条件モードがOFFの場合は省略可能。

-  [ int09 ]

   **形式 :** int型 (空白不可)

   **説明 :**
   変分パラメータの種類を表します(最適化有無の設定用)。0以上[int01]未満で指定します。

-  [ int10 ]

   **形式 :** int型 (空白不可)

   **説明 :**
   [int06]で指定した変分パラメータの最適化有無を設定します。最適化する場合は1、最適化しない場合は0とします。

使用ルール
^^^^^^^^^^

本ファイルを使用するにあたってのルールは以下の通りです。

-  行数固定で読み込みを行う為、ヘッダの省略はできません。

-  [ int01 ] と定義されている変分パラメータの種類の総数が異なる場合はエラー終了します。

-  [ int02 ]-[ int10 ] を指定する際、範囲外の整数を指定した場合はエラー終了します。

TransSym指定ファイル(qptransidx.def)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

運動量射影 :math:`{\cal L}_K=\frac{1}{N_s}\sum_{{\boldsymbol R}}e^{i {\boldsymbol K} \cdot{\boldsymbol R} } \hat{T}_{\boldsymbol R}` と格子対称性射影 :math:`{\cal L}_P=\sum_{\alpha}p_{\alpha} \hat{G}_{\alpha}` について、重みとサイト番号に関する指定を行います。射影するパターンは :math:`(\alpha, {\boldsymbol R})` で指定されます。射影を行わない場合も重み1.0
で“恒等演算”を指定してください。 以下にファイル例を記載します。

::

    ====================================
    NQPTrans 4
    ====================================
    == TrIdx_TrWeight_and_TrIdx_i_xi  ==
    ====================================
       0  1.000000  0.000000
       1  1.000000  0.000000
       2  1.000000  0.000000
       3  1.000000  0.000000
       0     0    0
     (continue...)
       3    12    1
       3    13    2

ファイル形式
^^^^^^^^^^^^

以下のように行数に応じ異なる形式をとります( :math:`N_s` はサイト数、 :math:`N_{\rm TS}` は射影演算子の種類の総数)。

-  1行: ヘッダ(何が書かれても問題ありません)。

-  2行: [string01] [int01]

-  3-5行: ヘッダ(何が書かれても問題ありません)。

-  6 - (5+ :math:`N_{\rm TS}`)行: [int02] [double01] [double02]

-  (6+ :math:`N_{\rm TS}`) -
   (5+ :math:`(N_s+1) \times N_{\rm TS}`)行：[int03] [int04] [int05] [int06]

パラメータ
^^^^^^^^^^

-  [ string01 ]

   **形式 :** string型 (空白不可)

   **説明 :** 射影パターンの総数に関するキーワード名を指定します(任意)。

-  [ int01 ]

   **形式 :** int型 (空白不可)

   **説明 :** 射影パターンの総数を指定します。
   正の値を指定し、``ModPara`` の ``abs(NMPTrans)`` 以上でなければなりません。
   実際に使うのは先頭 ``abs(NMPTrans)`` 個です。

-  [ int02 ]

   **形式 :** int型 (空白不可)

   **説明 :**
   射影パターン :math:`(\alpha, {\boldsymbol R})` を指定する整数。0以上
   [ int01 ] 未満で指定します。

-  [ double01 ]

   **形式 :** double型 (空白不可)

   **説明 :**
   射影パターン :math:`(\alpha, {\boldsymbol R})` の重み :math:`p_{\alpha}\cos ({\boldsymbol K}\cdot {\boldsymbol R})` の実部を指定します。

-  [ double02 ]

   **形式 :** double型

   **説明 :**
   射影パターン :math:`(\alpha, {\boldsymbol R})` の重み :math:`p_{\alpha}\exp (i{\boldsymbol K}\cdot {\boldsymbol R})` の虚部を指定します。

-  [ int03 ]

   **形式 :** int型 (空白不可)

   **説明 :**
   射影パターン :math:`(\alpha, {\boldsymbol R})` を指定する整数。0以上
   [ int01 ] 未満で指定します。

-  [ int04 ], [ int05 ]

   **形式 :** int型 (空白不可)

   **説明 :**
   サイト番号を指定する整数。0以上 ``Nsite`` 未満で指定します。 [ int03 ] で指定した並進・点群移動をサイト番号 [ int04 ] に作用させた場合の行き先が、サイト番号 [ int05 ] となるように設定します。

-  [ int06 ]

   **形式 :** int型 (空白不可)

   **説明 :**
   反周期境界条件モードがON( ``ModPara`` ファイルで ``NMPTrans`` が負の場合に有効)の場合、並進演算で生成消滅演算子の符号が反転するか否かを直接指定する。 [ int06 ] = :math:`\pm1` により符号を指定する。反周期境界条件モードがOFFの場合は省略可能。

使用ルール
^^^^^^^^^^

本ファイルを使用するにあたってのルールは以下の通りです。

-  行数固定で読み込みを行う為、ヘッダの省略はできません。

-  [ int01 ] と定義されている射影パターンの総数が異なる場合はエラー終了します。

-  [ int02 ]-[ int06 ] を指定する際、範囲外の整数を指定した場合はエラー終了します。

.. _InputParam:

変分パラメータ初期値指定ファイル
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

各変分パラメータの初期値を設定することが可能です。
変分パラメータの種類は :ref:`InputFileList` において
``InGutzwiller``, ``InJastrow``, ``InSpinJastrow``, ``InDH2``, ``InDH4``,  ``InGeneralRBM_PhysLayer``, ``InGeneralRBM_HiddenLayer``, ``InGeneralRBM_PhysHidden``, ``InOrbital``,
``InOrbitalAntiParallel``, ``InOrbitalParallel``, ``InOrbitalGeneral``
をキーワードとして指定することで区別します。なお、ファイルフォーマットは全て共通です。
以下、 ``InJastrow`` ファイルの例を記載します。
``InSpinJastrow`` では2行目にスピンJastrow変分パラメータ数を指定します
(例: ``NSpinJastrowIdx``)。

::

    ======================
    NJastrowIdx  28
    ======================
    == i_j_JastrowIdx  ===
    ======================
    0 -8.909963465082626488e-02  0.000000000000000000e+00
    1  5.521681211878626955e-02  0.000000000000000000e+00
    (continue...)
    27 -9.017586139930480749e-02  0.000000000000000000e+00

ファイル形式
^^^^^^^^^^^^

以下のように行数に応じ異なる形式をとります( :math:`N_v` は変分パラメータの種類の総数)。

-  1行: ヘッダ(何が書かれても問題ありません)。

-  2行: [string01] [int01]

-  3-5行: ヘッダ(何が書かれても問題ありません)。

-  6 - (5+ :math:`N_v`)行: [int02] [double01] [double02]

パラメータ
^^^^^^^^^^

-  [ string01 ]

   **形式 :** string型 (空白不可)

   **説明 :** 変分パラメータ総数のキーワード名を指定します(任意)。

-  [ int01 ]

   **形式 :** int型 (空白不可)

   **説明 :** 変分パラメータ総数を指定します。

-  [ int02 ]

   **形式 :** int型 (空白不可)

   **説明 :** 変分パラメータの種類を指定する整数。0以上
   [ int01 ] で指定します。

-  [ double01 ]

   **形式 :** double型 (空白不可)

   **説明 :** 変分パラメータの初期値の実部を与えます。

-  [ double02 ] **形式 :** double型 (空白不可)

   **説明 :** 変分パラメータの初期値の虚部を与えます。

使用ルール
^^^^^^^^^^

本ファイルを使用するにあたってのルールは以下の通りです。

- RBM相関因子のインプット機能は、現バージョンではベータ版です。今後のバージョンアップにより、ファイル形式や実装が変更される可能性があります。本機能を使用する際は、十分な検証を行ってください。

-  行数固定で読み込みを行う為、ヘッダの省略はできません。

-  [ int01 ] と定義されている変分パラメータの総数が異なる場合はエラー終了します。

OneBodyG指定ファイル(greenone.def)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

計算・出力する一体グリーン関数 :math:`\langle c_{i\sigma_1}^{\dagger}c_{j\sigma_2}\rangle` を指定します。以下にファイル例を記載します。

::

    ===============================
    NCisAjs         24
    ===============================
    ======== Green functions ======
    ===============================
        0     0     0     0
        0     1     0     1
        1     0     1     0
        1     1     1     1
        2     0     2     0
        2     1     2     1
        3     0     3     0
        3     1     3     1
        4     0     4     0
        4     1     4     1
        5     0     5     0
        5     1     5     1
        6     0     6     0
        6     1     6     1
        7     0     7     0
        7     1     7     1
        8     0     8     0
        8     1     8     1
        9     0     9     0
        9     1     9     1
       10     0    10     0
       10     1    10     1
       11     0    11     0
       11     1    11     1

ファイル形式
^^^^^^^^^^^^

以下のように行数に応じ異なる形式をとります。

-  1行: ヘッダ(何が書かれても問題ありません)。

-  2行: [string01] [int01]

-  3-5行: ヘッダ(何が書かれても問題ありません)。

-  6行以降: [int02]  [int03]  [int04]  [int05]

パラメータ
^^^^^^^^^^

-  [ string01 ]

   **形式 :** string型 (空白不可)

   **説明 :** 一体グリーン関数成分総数のキーワード名を指定します(任意)。

-  [ int01 ]

   **形式 :** int型 (空白不可)

   **説明 :** 一体グリーン関数成分の総数を指定します。

-  [ int02 ], [ int04 ]

   **形式 :** int型 (空白不可)

   **説明 :**
   サイト番号を指定する整数。0以上 ``Nsite`` 未満で指定します。

-  [ int03 ], [ int05 ]

   **形式 :** int型 (空白不可)

   | **説明 :** スピンを指定する整数。
   | 0: アップスピン
   | 1: ダウンスピン
   | を選択することが出来ます。

使用ルール
^^^^^^^^^^

本ファイルを使用するにあたってのルールは以下の通りです。

-  行数固定で読み込みを行う為、ヘッダの省略はできません。

-  [ int01 ] と定義されている一体グリーン関数成分の総数が異なる場合はエラー終了します。

-  [ int02 ]-[ int05 ] を指定する際、範囲外の整数を指定した場合はエラー終了します。

TwoBodyG指定ファイル(greentwo.def)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

計算・出力する二体グリーン関数 :math:`\langle c_{i\sigma_1}^{\dagger}c_{j\sigma_2}c_{k\sigma_3}^{\dagger}c_{l\sigma_4}\rangle` を指定します。
以下にファイル例を記載します。

::

    =============================================
    NCisAjsCktAltDC        576
    =============================================
    ======== Green functions for Sq AND Nq ======
    =============================================
        0     0     0     0     0     0     0     0
        0     0     0     0     0     1     0     1
        0     0     0     0     1     0     1     0
        0     0     0     0     1     1     1     1
        0     0     0     0     2     0     2     0
        0     0     0     0     2     1     2     1
        0     0     0     0     3     0     3     0
        0     0     0     0     3     1     3     1
        0     0     0     0     4     0     4     0
        0     0     0     0     4     1     4     1
        0     0     0     0     5     0     5     0
        0     0     0     0     5     1     5     1
        0     0     0     0     6     0     6     0
        0     0     0     0     6     1     6     1
        0     0     0     0     7     0     7     0
        0     0     0     0     7     1     7     1
        0     0     0     0     8     0     8     0
        0     0     0     0     8     1     8     1
        0     0     0     0     9     0     9     0
        0     0     0     0     9     1     9     1
        0     0     0     0    10     0    10     0
        0     0     0     0    10     1    10     1
        0     0     0     0    11     0    11     0
        0     0     0     0    11     1    11     1
        0     1     0     1     0     0     0     0
        …

ファイル形式
^^^^^^^^^^^^

以下のように行数に応じ異なる形式をとります。

-  1行: ヘッダ(何が書かれても問題ありません)。

-  2行: [string01] [int01]

-  3-5行: ヘッダ(何が書かれても問題ありません)。

-  6行以降:
   [int02]  [int03]  [int04]  [int05]  [int06]  [int07]  [int08]  [int09]

パラメータ
^^^^^^^^^^

-  [ string01 ]

   **形式 :** string型 (空白不可)

   **説明 :** 二体グリーン関数成分総数のキーワード名を指定します(任意)。

-  [ int01 ]

   **形式 :** int型 (空白不可)

   **説明 :** 二体グリーン関数成分の総数を指定します。

-  [ int02 ], [ int04 ],
   [ int06 ], [ int08 ]

   **形式 :** int型 (空白不可)

   **説明 :**
   サイト番号を指定する整数。0以上 ``Nsite`` 未満で指定します。

-  [ int03 ], [ int05 ],
   [ int07 ], [ int09 ]

   **形式 :** int型 (空白不可)

   | **説明 :** スピンを指定する整数。
   | 0: アップスピン
   | 1: ダウンスピン
   | を選択することが出来ます。

使用ルール
^^^^^^^^^^

本ファイルを使用するにあたってのルールは以下の通りです。

-  行数固定で読み込みを行う為、ヘッダの省略はできません。

-  [ int01 ] と定義されている二体グリーン関数成分の総数が異なる場合はエラー終了します。

-  [ int02 ]-[ int09 ] を指定する際、範囲外の整数を指定した場合はエラー終了します。

.. [3]
   使用メモリ量が、 :math:`O(N_\text{p}^2)` から
   :math:`O(N_\text{p}^2) + O(N_\text{p}N_\text{MCS})` になります。

.. [4]
   使用メモリ量は、 :math:`O(N_\text{p}) + O(N_\text{p}N_\text{MCS})`
   です。

NBodyG指定ファイル(nbodyg.def)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

計算・出力する :math:`N` 体相関関数
:math:`\langle \prod_{a=1}^{N} c_{i_a\sigma_a}^{\dagger} c_{j_a\tau_a} \rangle`
を指定します。以下にファイル例を記載します。

::

    =============================================
    NNBodyG        3
    =============================================
    ======== NBodyG correlation functions =======
    =============================================
        1     0     0     0     0
        2     0     0     0     0     1     1     1     1
        3     0     0     0     0     1     0     1     0     2     1     2     1

ファイル形式
^^^^^^^^^^^^

以下のように行数に応じ異なる形式をとります。

-  1行: ヘッダ(何が書かれても問題ありません)。

-  2行: [string01] [int01]

-  3-5行: ヘッダ(何が書かれても問題ありません)。

-  6行以降:
   [N] [i_1] [sigma_1] [j_1] [tau_1] ... [i_N] [sigma_N] [j_N] [tau_N]

パラメータ
^^^^^^^^^^

-  [ string01 ]

   **形式 :** string型 (空白不可)

   **説明 :** :math:`N` 体相関関数成分総数のキーワード名を指定します(任意)。

-  [ int01 ]

   **形式 :** int型 (空白不可)

   **説明 :** :math:`N` 体相関関数成分の総数を指定します。

-  [ N ]

   **形式 :** int型 (空白不可)

   **説明 :** 相関関数の次数を指定します。正の整数である必要があります。

-  [ i_a ], [ j_a ]

   **形式 :** int型 (空白不可)

   **説明 :** サイト番号を指定する整数。0以上 ``Nsite`` 未満で指定します。

-  [ sigma_a ], [ tau_a ]

   **形式 :** int型 (空白不可)

   | **説明 :** スピンを指定する整数。
   | 0: アップスピン
   | 1: ダウンスピン。

使用ルール
^^^^^^^^^^

本ファイルを使用するにあたってのルールは以下の通りです。

-  行数固定で読み込みを行う為、ヘッダの省略はできません。

-  [ int01 ] と定義されている :math:`N` 体相関関数成分の総数が異なる場合はエラー終了します。

-  各成分行は :math:`1 + 4N` 個の整数を持つ必要があります。

-  サイト番号またはスピン番号に範囲外の整数を指定した場合はエラー終了します。

-  スピンを変える因子 :math:`\sigma_a \neq \tau_a` は orbital-general モードでのみ指定できます。
   それ以外のモードでは各因子が :math:`\sigma_a = \tau_a` を満たす必要があります。

-  本出力はBackFlowの有無によらず物理量計算で使用できます。BackFlow
   :math:`N` 体測定では複素変分パラメータ、``NLanczosMode=0``、
   ``reweight=0`` が必要です。real BackFlow :math:`N` 体kernelと
   Lanczos補正付き ``ls_NBodyG`` 出力は未実装です。

-  通常の ``Orbital`` / ``OrbitalAntiParallel`` BackFlowでは全因子が
   スピンを保存する必要があります。``OrbitalGeneral`` / FSZで
   スピンを変える因子を使用できるのは ``2Sz=-1`` の場合だけです。

-  入力次数には任意の正の整数を指定できます。縮約後、non-FSZ BackFlowは
   実効次数1をone-body kernelへdispatchし、実効次数2以上ではcandidate状態を
   完全に再構築します。BF-FSZは検証済みの次数1/2 dispatchを使用し、
   それより高い次数を再構築します。真正な実効次数3以上では、残った各成分ごとに
   BackFlow Slater/Pfaffianを完全構築する計算量が必要です。

-  出力ファイル名は従来どおり ``xxx_NBodyG_%03d.dat`` です。

Twist指定ファイル(twist.def)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

計算・出力するTwist演算子 :math:`P^{(\alpha)} = \langle \exp ( i 2\pi \sum_{i\sigma} \sum_{\mu=x,y,z} c^{(\alpha)\mu }_{i\sigma } \mu_{i} n_{i\sigma} ) \rangle` を指定します。位置演算子 :math:`\mu_i` はLattice指定ファイル(lattice.def)で定義します。
以下にファイル例を記載します。

::

    --------------------
    NTwist	2
    --------------------
    idx_site_s_x_y_z
    --------------------
    0	0	0	0.0	0.0	0.0
    0	1	0	0.0	0.0	0.0
    0	2	0	0.0	0.0	0.0
    0	3	0	0.0	0.0	0.0
    0	4	0	0.0	0.0	0.0
    0	5	0	0.0	0.0	0.0
    0	0	1	0.3333333333	0.0	0.0
    0	1	1	0.3333333333	0.0	0.0
    0	2	1	0.3333333333	0.0	0.0
    0	3	1	0.3333333333	0.0	0.0
    0	4	1	0.3333333333	0.0	0.0
    0	5	1	0.3333333333	0.0	0.0
    1	0	0	0.3333333333	0.0	0.0
    1	1	0	0.3333333333	0.0	0.0
        …

ファイル形式
^^^^^^^^^^^^

以下のように行数に応じ異なる形式をとります。

-  1行: ヘッダ(何が書かれても問題ありません)。

-  2行: [string01] [int01]

-  3-5行: ヘッダ(何が書かれても問題ありません)。

-  6行以降:
   [int02]  [int03]  [int04]  [double01]  [double02]  [double03]

パラメータ
^^^^^^^^^^

-  [ string01 ]

   **形式 :** string型 (空白不可)

   **説明 :** Twist演算子の総数のキーワード名を指定します(任意)。

-  [ int01 ]

   **形式 :** int型 (空白不可)

   **説明 :** Twist演算子の総数を指定します。

-  [ int02 ]

   **形式 :** int型 (空白不可)

   **説明 :**
   Twist演算子の番号 :math:`\alpha` を指定する整数。0以上 [ int01 ] 未満で指定します。

-  [ int03 ]

   **形式 :** int型 (空白不可)

   **説明 :**
   サイト番号を指定する整数。0以上 ``Nsite`` 未満で指定します。

-  [ int04 ]

   **形式 :** int型 (空白不可)

   | **説明 :** スピンを指定する整数。
   | 0: アップスピン
   | 1: ダウンスピン
   | を選択することが出来ます。

-  [ double01 ], [ double02 ],
   [ double03 ]

   **形式 :** double型 (空白不可)

   **説明 :** 各座標方向 :math:`\mu = x, y, z` に対応するTwist成分
   :math:`c^{(\alpha)\mu}_{i\sigma}` を指定します。 :math:`x` 方向成分を [ double01 ]、 :math:`y` 方向成分を [ double02 ]、 :math:`z` 方向成分を [ double03 ] に指定します。

使用ルール
^^^^^^^^^^

本ファイルを使用するにあたってのルールは以下の通りです。

-  行数固定で読み込みを行う為、ヘッダの省略はできません。

-  各Twist演算子を指定するにあたり、全てのサイト・スピンの組み合わせを指定する必要があります。

Lattice指定ファイル(lattice.def)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

サイト番号 :math:`i` に対応する位置演算子 :math:`\mu_i` とユニットセル内の軌道番号を定義します。
以下にファイル例を記載します。

::

    --------------------
    NLattice  4 4 4 2
    --------------------
    i_x_y_z_orb
    --------------------
    0 0 0 0 0
    1 0 0 0 1
    2 1 0 0 0
    3 1 0 0 1
    4 2 0 0 0
    5 2 0 0 1
        …

ファイル形式
^^^^^^^^^^^^

以下のように行数に応じ異なる形式をとります。

-  1行: ヘッダ(何が書かれても問題ありません)。

-  2行: [string01] [int01] [int02]  [int03]  [int04]

-  3-5行: ヘッダ(何が書かれても問題ありません)。

-  6行以降:
   [int05]  [int06]  [int07]  [int08]  [int09]

パラメータ
^^^^^^^^^^

-  [ string01 ]

   **形式 :** string型 (空白不可)

   **説明 :** キーワード名を指定します(任意)。

-  [ int01 ]

   **形式 :** int型 (空白不可)

   **説明 :** 位置演算子 :math:`\mu` の :math:`x` 方向の最大値を指定します。

-  [ int02 ]

   **形式 :** int型 (空白不可)

   **説明 :** 位置演算子 :math:`\mu` の :math:`y` 方向の最大値を指定します。

-  [ int03 ]

   **形式 :** int型 (空白不可)

   **説明 :** 位置演算子 :math:`\mu` の :math:`z` 方向の最大値を指定します。

-  [ int04 ]

   **形式 :** int型 (空白不可)

   | **説明 :** ユニットセル内の軌道数を指定します。

-  [ int05 ]

   **形式 :** int型 (空白不可)

   **説明 :**
   サイト番号を指定する整数。0以上 ``Nsite`` 未満で指定します。

-  [ int06 ], [ int07 ],
   [ int08 ]

   **形式 :** int型 (空白不可)

   **説明 :** サイト番号 [ int05 ] に対応するユニットセル位置演算子 :math:`\mu_i` の各座標方向 :math:`\mu = x, y, z` の値を指定します。 :math:`x, y, z` 方向成分をそれぞれ 0以上かつ [ int01 ], [ int02 ], [ int03 ] 未満で指定します。

-  [ int09 ]

   **形式 :** int型 (空白不可)

   **説明 :** サイト番号 [ int05 ] に対応するユニットセル内の軌道番号を指定します。0以上 [ int04 ] 未満で指定します。この値はスピン自由度ではありません。

使用ルール
^^^^^^^^^^

本ファイルを使用するにあたってのルールは以下の通りです。

-  行数固定で読み込みを行う為、ヘッダの省略はできません。

-  全てのサイト番号を指定する必要があります。

-  各行は 1 つのサイトの格子座標を指定します。アップスピン・ダウンスピンに
   対して lattice.def の行を重複して指定しないでください。スピンに依存する
   twist 成分は twist.def で指定します。

-  通常の lattice.def は各サイトに対して 1 つの座標を指定する形式なので、
   ``[int01] * [int02] * [int03] * [int04]`` (すなわち
   ``Nx * Ny * Nz * Norb``) は通常 ``Nsite`` と一致します。一致しない場合も
   互換性のため計算は継続しますが、 warning が出力されます。その場合は
   lattice.def の座標と twist.def の位相が意図した幾何を表すように設定してください。

-  各サイト番号は一度だけ指定する必要があります。
   重複や欠番は入力エラーとして拒否されます。
