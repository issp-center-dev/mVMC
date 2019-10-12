チュートリアル
==============

サンプルファイル一覧
--------------------

mVMCでは ``sample/Standard/`` 以下に次のサンプルを用意しています。

-  2次元正方格子Hubbardモデル

   (``sample/Standard/Hubbard/square/``)

-  2次元三角格子Hubbardモデル

   (``sample/Standard/Hubbard/triangular/``)

-  1次元近藤格子モデル

   (``sample/Standard/Kondo/chain/``)

-  1次元反強磁性的Heisenbergモデル

   (``sample/Standard/Spin/HeisenbergChain/``)

-  2次元正方格子反強磁性的Heisenbergモデル

   (``sample/Standard/Spin/HeisenbergSquare/``)

-  2次元カゴメ格子反強磁性的Heisenbergモデル

   (``sample/Standard/Spin/Kagome/``)

これらのチュートリアルの実行方法は全て同じ手順で実行することが可能です。
以下ではHeisenberg模型について説明します。

Heisenberg模型
--------------

以下のチュートリアルはディレクトリ

::

    sample/Standard/Spin/HeisenbergChain/

内で行います。 このディレクトには以下のファイルがあります.

Heisenberg模型におけるサンプル入力ファイル: ``StdFace.def``

参照用出力ディレクトリ: ``reference/``

この例では1次元のHeisenberg鎖(最近接サイト間の反強磁性的スピン結合のみを持つ)を考察します。

.. math::

   \begin{aligned}
     {\hat H} = J \sum_{i=1}^{L} {\hat {\boldsymbol S}}_i \cdot {\hat {\boldsymbol S}}_{i+1}\end{aligned}

| インプットファイルの中身は次のとおりです。

10cm

::

    L = 16
    Lsub=4
    model = "Spin"
    lattice = "chain lattice"
    J = 1.0
    2Sz = 0
    NMPtrans=1

| 
| この例ではスピン結合 :math:`J=1` (任意単位)とし、サイト数は16としました。

計算実行
^^^^^^^^

実行コマンドは次のとおりです。

``$ mpiexec -np `` `` `` ``/vmc.out -s StdFace.def``

使っているシステムによっては ``mpiexec`` コマンドではなく ``mpirun`` や ``mpijob`` 、
``poe`` となる場合もあります。

この実行による標準出力は次のとおりです。

::

    ######  Standard Intarface Mode STARTS  ######

      Open Standard-Mode Inputfile StdFace.def 

      KEYWORD : l                    | VALUE : 16 
      KEYWORD : lsub                 | VALUE : 4 
      KEYWORD : model                | VALUE : spin 
      KEYWORD : lattice              | VALUE : chain 
      KEYWORD : j                    | VALUE : 1.0 
      KEYWORD : nmptrans             | VALUE : 1 

    #######  Parameter Summary  #######

      @ Lattice Size & Shape

                    L = 16 
                 Lsub = 4         
                    L = 16        
                    W = 1         
               phase0 = 0.00000     ######  DEFAULT VALUE IS USED  ######

      @ Hamiltonian 

                   2S = 1           ######  DEFAULT VALUE IS USED  ######
                    h = 0.00000     ######  DEFAULT VALUE IS USED  ######
                Gamma = 0.00000     ######  DEFAULT VALUE IS USED  ######
                    D = 0.00000     ######  DEFAULT VALUE IS USED  ######
                  J0x = 1.00000   
                  J0y = 1.00000   
                  J0z = 1.00000   

      @ Numerical conditions

                 Lsub = 4         
                 Wsub = 1         
          ioutputmode = 1           ######  DEFAULT VALUE IS USED  ######

    ######  Print Expert input files  ######

        qptransidx.def is written.
             filehead = zvo         ######  DEFAULT VALUE IS USED  ######
             filehead = zqp         ######  DEFAULT VALUE IS USED  ######
          NVMCCalMode = 0           ######  DEFAULT VALUE IS USED  ######
         NLanczosMode = 0           ######  DEFAULT VALUE IS USED  ######
        NDataIdxStart = 1           ######  DEFAULT VALUE IS USED  ######
          NDataQtySmp = 1           ######  DEFAULT VALUE IS USED  ######
          NSPGaussLeg = 8           ######  DEFAULT VALUE IS USED  ######
             NMPTrans = 1         
        NSROptItrStep = 1000        ######  DEFAULT VALUE IS USED  ######
         NSROptItrSmp = 100         ######  DEFAULT VALUE IS USED  ######
           NVMCWarmUp = 10          ######  DEFAULT VALUE IS USED  ######
         NVMCInterval = 1           ######  DEFAULT VALUE IS USED  ######
           NVMCSample = 1000        ######  DEFAULT VALUE IS USED  ######
        NExUpdatePath = 2         
              RndSeed = 123456789   ######  DEFAULT VALUE IS USED  ######
           NSplitSize = 1           ######  DEFAULT VALUE IS USED  ######
               NStore = 0           ######  DEFAULT VALUE IS USED  ######
         DSROptRedCut = 0.00100     ######  DEFAULT VALUE IS USED  ######
         DSROptStaDel = 0.02000     ######  DEFAULT VALUE IS USED  ######
         DSROptStepDt = 0.02000     ######  DEFAULT VALUE IS USED  ######
              NSPStot = 0           ######  DEFAULT VALUE IS USED  ######
          ComplexType = 0           ######  DEFAULT VALUE IS USED  ######
        locspn.def is written.
        trans.def is written.
        interall.def is written.
        jastrowidx.def is written.
        coulombintra.def is written.
        coulombinter.def is written.
        hund.def is written.
        exchange.def is written.
        orbitalidx.def is written.
        gutzwilleridx.def is written.
        namelist.def is written.
        modpara.def is written.
        greenone.def is written.
        greentwo.def is written.

    ######  Input files are generated.  ######

    -----------
    Start: Read *def files.
      Read File namelist.def .
      Read File 'modpara.def' for ModPara.
      Read File 'locspn.def' for LocSpin.
      Read File 'trans.def' for Trans.
      Read File 'coulombintra.def' for CoulombIntra.
      Read File 'coulombinter.def' for CoulombInter.
      Read File 'hund.def' for Hund.
      Read File 'exchange.def' for Exchange.
      Read File 'gutzwilleridx.def' for Gutzwiller.
      Read File 'jastrowidx.def' for Jastrow.
      Read File 'orbitalidx.def' for Orbital.
      Read File 'qptransidx.def' for TransSym.
      Read File 'greenone.def' for OneBodyG.
      Read File 'greentwo.def' for TwoBodyG.
    End  : Read *def files.
    Start: Read parameters from *def files.
    End  : Read parameters from *def files.
    Start: Set memories.
    End  : Set memories.
    Start: Initialize parameters.
    End  : Initialize parameters.
    Start: Initialize variables for quantum projection.
    End  : Initialize variables for quantum projection.
    Start: Optimize VMC parameters.
    End  : Optimize VMC parameters.
    -----------

この実行でははじめにエキスパートモード用の入力ファイルとして、
ハミルトニアンの詳細を記述するファイル

-  ``locspin.def``

-  ``trans.def``

-  ``coulombinter.def``

-  ``coulombintra.def``

-  ``exchange.def``

-  ``hund.def``

-  ``namelist.def``

-  ``modpara.def``

と、変分パラメータを設定するファイル

-  ``gutzwilleridx.def``

-  ``jastrowidx.def``

-  ``orbitalidx.def``

-  ``qptransidx.def``

結果として出力する相関関数の要素を指定するファイル

-  ``greenone.def``

-  ``greentwo.def``

が生成されます。 各ファイルの詳細については :ref:`HowToExpert` をご覧ください。

その後実際に計算が行われ、以下のファイルが情報として ``output/`` ディレクトリに出力されます。

12cm

::

    zvo_SRinfo.dat
    zvo_out_001.dat
    zvo_time_001.dat
    zvo_var_001.dat
    zvo_CalcTimer.dat

なお、 ``zvo_out_001.dat`` には、ビン毎の計算情報として、

.. math:: \langle H \rangle, \langle H^2 \rangle, \frac{\langle H^2 \rangle- \langle H \rangle^2 }{\langle H \rangle^2} \nonumber

が順に出力されますので、収束性の目安として利用することが可能です。
gnuplotを用いる場合には、次のようにして表示することが出来ます( :math:`\langle H \rangle` の場合)。

::

    plot "zvo_out_001.dat" u 1

| 各ファイルの詳細については :ref:`outputfile` をご覧ください。

計算結果出力
^^^^^^^^^^^^

| 計算が正常終了すると、エネルギー、エネルギーの分散、
  変分パラメータおよび計算実行時間を記載したファイルが ``output/``
  ディレクトリに出力されます。
  以下に、このサンプルでの出力ファイルを記載します。

12cm

::

    gutzwiller_opt.dat
    jastrow_opt.dat
    orbital_opt.dat
    zqp_opt.dat
    ClacTimer.dat

各ファイルの詳細については :ref:`outputfile` をご覧ください。

Green関数の計算
^^^^^^^^^^^^^^^

``modpara.def`` ファイル中の ``NVMCCalMode`` を0から1に変更の上、以下のコマンドを実行します。
下記のように 実行時のコマンドライン引数として
``"namelist.dat"`` の後ろに ``"zqp_opt.dat"`` を付け加えることで、
一つ前の計算で最適化された変分パラメータを使用した計算が行われます。

``$ `` ``/vmc.out -e namelist.def output/zqp_opt.dat``

| 計算が終了すると以下のファイルが ``output/``
  ディレクトリに出力されます。

12cm

::

    zvo_cisajs_001.dat
    zvo_cisajscktalt_001.dat

| 
| 各ファイルの詳細については :ref:`outputfile` をご覧ください。

エキスパートユーザー向け
------------------------

mVMCでは、以下の6つに分類される入力ファイルを読み込み、計算実行を行います。

(1) List:
    詳細入力ファイルの種類と名前を指定するファイル

(2) Basic parameters:
    基本的なパラメータを指定するファイル

(3) Set Hamiltonian:
    ハミルトニアンを指定するファイル

(4) Set condition of variational parameters :
    最適化する変分パラメータを指定するファイル

(5) Initial variational parameters:
    変分パラメータの初期値を指定するファイル

(6) Output:
    出力する一体・二体グリーン関数の成分を指定するファイル

上記で分類されるファイルを直接作成・指定することで、より複雑な計算を行うことが可能です。
ファイルの詳細については :ref:`HowToExpert` をご覧ください。

相関関数のフーリエ変換
----------------------

このパッケージには、上で求めた相関関数をフーリエ変換し、プロットするユーティリティーが付属しています。
詳しくは :ref:`fourier` を参照してください。
