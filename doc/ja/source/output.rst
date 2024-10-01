.. include:: ../../bib/ref.txt

.. _outputfile:
             
出力ファイル
------------

出力ファイルの一覧は下記の通りです。\*\*\*には ``ModPara`` ファイルの ``CParaFileHead`` で指定されるヘッダが、xxxには ``CDataFileHead`` で指定されるヘッダが、yyyには ``ModPara`` ファイルの ``NDataIdxStart``,
``NDataQtySmp`` に従い ``NDataIdxStart`` :math:`\cdots` ``NDataIdxStart`` + ``NDataQtySmp`` の順に記載されます。また、zzzには ``ModPara`` ファイルの ``NDataIdxStart`` が記載されます。

+------------------------------------------+----------------------------------------------------------------------+
| ファイル名                               | 対応するファイルの中身                                               |
+==========================================+======================================================================+
| \*\*\*\_opt.dat                          | 最適化された全パラメータ.                                            |
+------------------------------------------+----------------------------------------------------------------------+
| \*\*\*\_gutzwiller\_opt.dat              | 最適化されたGutzwiller因子.                                          |
+------------------------------------------+----------------------------------------------------------------------+
| \*\*\*\_jastrow\_opt.dat                 | 最適化されたJastrow因子.                                             |
+------------------------------------------+----------------------------------------------------------------------+
| \*\*\*\_doublonHolon2site\_opt.dat       | 最適化された2サイトダブロン-ホロン相関因子.                          |
+------------------------------------------+----------------------------------------------------------------------+
| \*\*\*\_doublonHolon4site\_opt.dat       | 最適化された4サイトダブロン-ホロン相関因子.                          |
+------------------------------------------+----------------------------------------------------------------------+
| \*\*\*\_generalRBM\_physlayer\_opt.dat   | 最適化されたRBM相関因子のうち、物理層のみに関わる変分パラメータ.     |
+------------------------------------------+----------------------------------------------------------------------+
| \*\*\*\_generalRBM\_hiddenlayer\_opt.dat | 最適化されたRBM相関因子のうち、隠れ層のみに関わる変分パラメータ.     |
+------------------------------------------+----------------------------------------------------------------------+
| \*\*\*\_generalRBM\_physhidden\_opt.dat  | 最適化されたRBM相関因子のうち、物理層と隠れ層を繋ぐ変分パラメータ.   |
+------------------------------------------+----------------------------------------------------------------------+
| \*\*\*\_orbital\_opt.dat                 | 最適化されたペア軌道因子.                                            |
+------------------------------------------+----------------------------------------------------------------------+
| xxx\_out\_yyy.dat                        | エネルギーとその分散.                                                |
+------------------------------------------+----------------------------------------------------------------------+
| xxx\_SRinfo.dat                          | SR法に関連した情報.                                                  |
+------------------------------------------+----------------------------------------------------------------------+
| xxx\_var\_yyy.dat                        | パラメータ最適化過程の情報.                                          |
+------------------------------------------+----------------------------------------------------------------------+
| xxx\_CalcTimer.dat                       | 各プロセスに対する計算時間に関する情報.                              |
+------------------------------------------+----------------------------------------------------------------------+
| xxx\_time\_zzz.dat                       | モンテカルロサンプリングの過程に関する情報.                          |
+------------------------------------------+----------------------------------------------------------------------+
| xxx\_cisajs\_yyy.dat                     | 一体グリーン関数.                                                    |
+------------------------------------------+----------------------------------------------------------------------+
| xxx\_cisajscktalt\_yyy.dat               | 二体グリーン関数.                                                    |
+------------------------------------------+----------------------------------------------------------------------+

変分パラメータ出力ファイル(\*\*\*\_opt.dat)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

SR法で最適化された変分パラメータとエネルギーが一斉出力されます。
変分パラメータが一斉に読み込めるため、変分パラメータ最適化後の物理量の計算を行う場合に使用します。
出力されるデータは

.. math:: \langle H \rangle, \langle H^2 \rangle, g_i, v_{ij}, \alpha_{2nt}^{d(h)}, \alpha_{4nt}^{d(h)}, a_{i\sigma}, b_{k}, W_{i\sigma k}, f_{ij} \nonumber


で、それぞれの平均値と標準偏差が出力されます(平均値は実数、虚数の順に、標準偏差は実数のみ出力)。
なお、全データが1行で出力され、数値の間は半角空白で区切られます。
\*\*\*には ``ModPara`` ファイルの ``CParaFileHead`` で指定されるヘッダが記載されます。

ステップ別変分パラメータ出力ファイル(xxx\_var\_yyy.dat)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

SR 法の各ステップにおける変分パラメータとエネルギーがzqp\_opt.dat
と同形式で“追記しながら” 出力されます(標準偏差はゼロ)が出力されます。
xxxには ``CDataFileHead`` で指定されるヘッダが、yyyには ``ModPara`` ファイルの ``NDataIdxStart``,
``NDataQtySmp`` に従い ``NDataIdxStart`` :math:`\cdots` ``NDataIdxStart`` + ``NDataQtySmp`` の順に記載されます。

Gutzwiller因子出力ファイル(\*\*\*\_gutzwiller\_opt.dat)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

SR法で最適化されたGutzwiller因子が出力されます。
出力形式は :ref:`InputParam` のInGutzwiller指定ファイルと同じです。

Jastrow因子出力ファイル(\*\*\*\_jastrow\_opt.dat)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

SR法で最適化されたJastrow因子が出力されます。
出力形式は :ref:`InputParam` のInJastrow指定ファイルと同じです。

2サイトダブロンホロン相関因子出力ファイル(\*\*\*\_doublonHolon2site\_opt.dat)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

SR法で最適化された2サイトのdoublon-holon相関因子が出力されます。
出力形式は :ref:`InputParam` のInDH2指定ファイルと同じです。

4サイトダブロンホロン相関因子出力ファイル(\*\*\*\_doublonHolon4site\_opt.dat)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

SR法で最適化された4サイトのdoublon-holon相関因子が出力されます。
出力形式は :ref:`InputParam` のInDH4指定ファイルと同じです。

RBM相関因子出力ファイル(\*\*\*\_generalRBM\_physlayer\_opt.dat, \*\*\*\_generalRBM\_hiddenlayer\_opt.dat, \*\*\*\_generalRBM\_physhidden\_opt.dat)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

SR法で最適化されたRBM相関因子が出力されます。
それぞれの出力形式は :ref:`InputParam` のInGeneralRBM_PhysLayer, InGeneralRBM_HiddenLayer, InGeneralRBM_PhysHidden指定ファイルと同じです。

ペア軌道出力ファイル(\*\*\*\_orbital\_opt.dat)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

SR法で最適化されたペア軌道の変分パラメータが出力されます。
出力形式は :ref:`InputParam` のInOrbital指定ファイルと同じです。

xxx\_out\_yyy.dat
~~~~~~~~~~~~~~~~~

ビン毎の計算情報として、

.. math::

   \langle H \rangle, \langle H^2 \rangle, \frac{\langle H^2 \rangle- \langle H \rangle^2 }{\langle H \rangle^2},
   \langle S^z \rangle, \langle (S^z)^2 \rangle
   \nonumber.


が順に出力されます。 :math:`\langle H \rangle` については実部と虚部がそれぞれ出力され、それ以外は実部のみ出力されます。
xxxには ``CDataFileHead`` で指定されるヘッダが、yyyには ``ModPara`` ファイルの ``NDataIdxStart``,
``NDataQtySmp`` に従い ``NDataIdxStart`` :math:`\cdots` ``NDataIdxStart`` + ``NDataQtySmp`` の順に記載されます。以下に出力例を記載します。

::

    1.151983765704212992e+01  8.124622418360909482e-01  \
    1.619082955438887268e+02  2.019905203939084959e-01 
    1.288482613817423150e+01  5.006903733262847433e-01   
    1.972000325276957824e+02  1.824505193695792893e-01
    1.308897206011880421e+01  5.701244886956570168e-01  \
    2.072610167083121837e+02  2.029162857569105916e-01
    …

xxx\_SRinfo.dat.dat
~~~~~~~~~~~~~~~~~~~~~

各最適化ステップごとに、
変分パラメータ数 ``Npara`` 、S行列の次元数 ``Msize`` 、 ``Orbital`` などの入力ファイルで最適化しないと指定した実パラメータ数 ``optCut`` 、 ``ModPara`` ファイルの ``DSROptRedCut`` の値により削減された実パラメータ数 ``diagCut`` 、S行列の対角成分の最大値 ``sDiagMax`` と最小値 ``sDiagMin`` 、変分パラメータの変化量の最大値 ``absRmax`` とそのパラメータインデックス ``imax`` といったSR法による最適化情報が順に出力されます。 
xxxには ``ModPara`` ファイルの ``CDataFileHead`` で指定されるヘッダが記載されます。以下に出力例を記載します。なお、変分パラメータが実変数( ``ComplexType=0`` )として扱われた場合、虚部のパラメータ数は ``optCut`` にカウントされます。

::

    #Npara Msize optCut diagCut sDiagMax  sDiagMin    absRmax      imax
    4     4     4     0  4.17626e-02  0.00000e+00 -1.60883e-01     4
    4     4     4     0  3.53941e-02  0.00000e+00  1.63056e-01     0
    4     4     4     0  3.28032e-02  0.00000e+00  1.69939e-01     0
    4     4     4     0  3.31451e-02  0.00000e+00  1.92363e-01     0
    …


xxx\_CalcTimer.dat 
~~~~~~~~~~~~~~~~~~~

計算終了後に、処理毎の計算処理時間が処理名、処理に割り当てられた識別番号、実行秒数の順に出力されます。出力例は以下の通りです。

::

    All                         [0]     15.90724
    Initialization              [1]      0.04357
      read options             [10]      0.00012
      ReadDefFile              [11]      0.00082
      SetMemory                [12]      0.00002
      InitParameter            [13]      0.03026
    VMCParaOpt                  [2]     15.86367
      VMCMakeSample             [3]     12.85650
        makeInitialSample      [30]      0.20219
        make candidate         [31]      0.02553
        hopping update         [32]     12.51967
          UpdateProjCnt        [60]      7.41864
          CalculateNewPfM2     [61]      3.67098
          CalculateLogIP       [62]      0.07599
          UpdateMAll           [63]      1.27466
        exchange update        [33]      0.00000
          UpdateProjCnt        [65]      0.00000
          CalculateNewPfMTwo2  [66]      0.00000
          CalculateLogIP       [67]      0.00000
          UpdateMAllTwo        [68]      0.00000
        recal PfM and InvM     [34]      0.08294
        save electron config   [35]      0.00232
      VMCMainCal                [4]      2.45481
        CalculateMAll          [40]      0.47556
        LocEnergyCal           [41]      0.79754
          CalHamiltonian0      [70]      0.00259
          CalHamiltonian1      [71]      0.18765
          CalHamiltonian2      [72]      0.00107
        ReturnSlaterElmDiff    [42]      0.40035
        calculate OO and HO    [43]      0.68045
      StochasticOpt             [5]      0.30489
        preprocess             [50]      0.02587
        stcOptMain             [51]      0.25471
          initBLACS            [55]      0.06564
          calculate S and g    [56]      0.05603
          DPOSV                [57]      0.09833
          gatherParaChange     [58]      0.02774
        postprocess            [52]      0.02372
      UpdateSlaterElm          [20]      0.02556
      WeightAverage            [21]      0.06676
      outputData               [22]      0.10554
      SyncModifiedParameter    [23]      0.02151

xxx\_time\_zzz.dat 
~~~~~~~~~~~~~~~~~~~

| 計算情報としてビン毎にサンプリング数、hoppingおよびexchangeのアップデートに対するacceptance
  ratio(acc\_hopp, acc\_ex)、それぞれのアップデートの試行回数(n\_hopp,
  n\_ex)および実行した際の時間を順に出力します。xxxには ``CDataFileHead`` で指定されるヘッダが、zzzには ``ModPara`` ファイルの ``NDataIdxStart`` が記載されます。
  出力例は以下の通りです。

::

    00000  acc_hop acc_ex  n_hop    n_ex     : Mon Jul 25 14:03:29 2016
    00001  0.59688 0.00000 320      0        : Mon Jul 25 14:03:30 2016
    00002  0.47727 0.00000 176      0        : Mon Jul 25 14:03:30 2016
    00003  0.50000 0.00000 176      0        : Mon Jul 25 14:03:30 2016
    00004  0.49432 0.00000 176      0        : Mon Jul 25 14:03:30 2016
    00005  0.57386 0.00000 176      0        : Mon Jul 25 14:03:30 2016
    00006  0.55114 0.00000 176      0        : Mon Jul 25 14:03:30 2016    
    …

xxx\_cisajs\_yyy.dat 
~~~~~~~~~~~~~~~~~~~~~

OneBodyG指定ファイルで指定された一体グリーン関数 :math:`\langle c_{i\sigma_1}^{\dagger}c_{j\sigma_2}\rangle` の計算結果を出力します。
xxxには ``CDataFileHead`` で指定されるヘッダが、yyyには ``ModPara`` ファイルの ``NDataIdxStart``,
``NDataQtySmp`` に従い ``NDataIdxStart`` :math:`\cdots` ``NDataIdxStart`` + ``NDataQtySmp`` の順に記載されます。以下にファイル例を記載します。

::

        0    0    0    0 0.4452776740 0.0000000000
        0    1    0    1 0.4452776740 0.0000000000
        1    0    1    0 0.5000000000 0.0000000000
        1    1    1    1 0.5000000000 0.0000000000
        2    0    2    0 0.4452776740 0.0000000000
        2    1    2    1 0.4452776740 0.0000000000
        3    0    3    0 0.5000000000 0.0000000000
        3    1    3    1 0.5000000000 0.0000000000
        …

ファイル形式
^^^^^^^^^^^^

-  :math:`[` int01 :math:`]`  :math:`[` int02 :math:`]`  :math:`[` int03 :math:`]`  :math:`[` int04 :math:`]`  :math:`[` double01 :math:`]`  :math:`[` double02 :math:`]`

パラメータ
^^^^^^^^^^

-  :math:`[` int01 :math:`]`, :math:`[` int03 :math:`]`

   **形式 :** int型

   **説明 :**
   サイト番号を指定する整数。 :math:`[` int01 :math:`]` が :math:`i` サイト、 :math:`[` int03 :math:`]` が :math:`j` サイトを表します。

-  :math:`[` int02 :math:`]`, :math:`[` int04 :math:`]`

   **形式 :** int型

   | **説明 :**
     スピンを指定する整数。 :math:`[` int02 :math:`]` が :math:`\sigma_1` 、 :math:`[` int03 :math:`]` が :math:`\sigma_2` に対応します。
   | 0: アップスピン
   | 1: ダウンスピン
   | を表します。

-  :math:`[` double01 :math:`]`, :math:`[` double02 :math:`]`

   **形式 :** double型

   | **説明 :**
     :math:`\langle c_{i\sigma_1}^{\dagger}c_{j\sigma_2}\rangle` の値を表します。
   | :math:`[` double01 :math:`]` が実部、 :math:`[` double02 :math:`]` が虚部を表します。

xxx\_cisajscktalt\_yyy.dat 
~~~~~~~~~~~~~~~~~~~~~~~~~~~

TwoBodyG指定ファイルで指定された二体グリーン関数 :math:`\langle c_{i\sigma_1}^{\dagger}c_{j\sigma_2}c_{k\sigma_3}^{\dagger}c_{l\sigma_4}\rangle` の計算結果を出力します。xxxには ``CDataFileHead`` で指定されるヘッダが、yyyには ``ModPara`` ファイルの ``NDataIdxStart``,
``NDataQtySmp`` に従い ``NDataIdxStart`` :math:`\cdots` ``NDataIdxStart`` + ``NDataQtySmp`` の順に記載されます。以下にファイル例を記載します。

15cm

::

        0    0    0    0    0    0    0    0 0.4452776740 0.0000000000
        0    0    0    0    0    1    0    1 0.1843355815 0.0000000000
        0    0    0    0    1    0    1    0 0.1812412105 0.0000000000
        0    0    0    0    1    1    1    1 0.2640364635 0.0000000000
        0    0    0    0    2    0    2    0 0.0279690007 0.0000000000
        0    0    0    0    2    1    2    1 0.2009271524 0.0000000000
        0    0    0    0    3    0    3    0 0.2512810778 0.0000000000
        0    0    0    0    3    1    3    1 0.1939965962 0.0000000000
        …

ファイル形式
^^^^^^^^^^^^

-  :math:`[` int01 :math:`]`  :math:`[` int02 :math:`]`  :math:`[` int03 :math:`]`  :math:`[` int04 :math:`]`  :math:`[` int05 :math:`]`  :math:`[` int06 :math:`]`  :math:`[` int07 :math:`]`  :math:`[` int08 :math:`]`  :math:`[` double01 :math:`]`  :math:`[` double02 :math:`]`

パラメータ
^^^^^^^^^^

-  :math:`[` int01 :math:`]`,
   :math:`[` int03 :math:`]`, :math:`[` int05 :math:`]`,
   :math:`[` int07 :math:`]`

   **形式 :** int型

   **説明 :** サイト番号を指定する整数。
   :math:`[` int01 :math:`]` が :math:`i` サイト、 :math:`[` int03 :math:`]` が :math:`j` サイト、 :math:`[` int05 :math:`]` が :math:`k` サイト、 :math:`[` int07 :math:`]` が :math:`l` サイトを表します。

-  :math:`[` int02 :math:`]`,
   :math:`[` int04 :math:`]`, :math:`[` int06 :math:`]`,
   :math:`[` int08 :math:`]`

   **形式 :** int型

   | **説明 :** スピンを指定する整数。
     :math:`[` int02 :math:`]` が :math:`\sigma_1` 、 :math:`[` int04 :math:`]` が :math:`\sigma_2` 、 :math:`[` int06 :math:`]` が :math:`\sigma_3` 、 :math:`[` int08 :math:`]` が :math:`\sigma_4` に対応します。
   | 0: アップスピン
   | 1: ダウンスピン
   | を表します。

-  :math:`[` double01 :math:`]`, :math:`[` double02 :math:`]`

   **形式 :** double型

   | **説明 :**
     :math:`\langle c_{i\sigma_1}^{\dagger}c_{j\sigma_2}c_{k\sigma_3}^{\dagger}c_{l\sigma_4}\rangle` の値を表します。
   | :math:`[` double01 :math:`]` が実部、 :math:`[` double02 :math:`]` が虚部を表します。

xxx\_ls\_out\_yyy.dat 
~~~~~~~~~~~~~~~~~~~~~~

Power Lanczos法により求めた :math:`\langle H \rangle`,
:math:`(\langle H^2\rangle - \langle H \rangle^2)/\langle H \rangle^2` および最適化パラメータ :math:`\alpha` の順に出力されます。
``ModPara`` 指定ファイルで ``NVMCCalMode`` =1, ``NLanczosmode`` =1,
2に設定することで計算されます。
xxxには ``CDataFileHead`` で指定されるヘッダが、yyyには ``ModPara`` ファイルの ``NDataIdxStart``,
``NDataQtySmp`` に従い ``NDataIdxStart`` :math:`\cdots` ``NDataIdxStart`` + ``NDataQtySmp`` の順に記載されます。

xxx\_ls\_cisajs\_yyy.dat 
~~~~~~~~~~~~~~~~~~~~~~~~~

Power Lanczos法により求めた,
OneBodyG指定ファイルで指定された一体グリーン関数 :math:`\langle c_{i\sigma_1}^{\dagger}c_{j\sigma_2}\rangle` の計算結果を出力します。
ファイル形式はxxx\_cisajs\_yyy.datファイルと同じです。 ``ModPara`` 指定ファイルで ``NVMCCalMode`` =1,
``NLanczosmode`` =2に設定することで計算されます。

xxx\_ls\_cisajscktalt\_yyy.dat 
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Power Lanczos法により求めた,
TwoBodyG指定ファイルで指定された二体グリーン関数 :math:`\langle c_{i\sigma_1}^{\dagger}c_{j\sigma_2}c_{k\sigma_3}^{\dagger}c_{l\sigma_4}\rangle` の計算結果を出力します。ファイル形式はxxx\_cisajscktalt\_yyy.datファイルと同じです。 ``ModPara`` 指定ファイルで ``NVMCCalMode``
= 1, ``NLanczosmode`` = 2 に設定することで計算されます。
