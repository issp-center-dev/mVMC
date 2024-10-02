How to use mVMC?
================

要件
----

mVMCのコンパイル :math:`\cdot` 使用には次のものが必要です。

-  Cコンパイラ (インテル、富士通、GNUなど)

-  MPIライブラリ

-  LAPACKライブラリ (インテルMKL, 富士通, ATLASなど)

-  オプション：ScaLAPACKライブラリ

.. note::

   **intelコンパイラーでの設定**

   intelコンパイラを使用する場合には、コンパイラに付属の設定用スクリプトを使用するのが簡単です。

   64ビットOSでbashを使っている場合には

   .. code-block:: bash

      $ source /opt/intel/bin/compilervars.sh intel64

   または

   .. code-block:: bash

      $ source /opt/intel/bin/iccvars.sh intel64
      $ source /opt/intel/mkl/bin/mklvars.sh

   等を ``~/.bashrc`` に記載してください。
   詳しくはお手持ちのコンパイラ、ライブラリのマニュアルをお読みください。

インストール方法
----------------

mVMC は次の場所からダウンロードできます。
https://github.com/issp-center-dev/mVMC/releases

ダウンロードしたファイルを次のように展開してください。

.. code-block:: bash

   $ tar xzvf mVMC-xxx.tar.gz

mVMCは次の2通りの方法でインストールできます。

``mVMCconfig.sh`` を使う方法
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

展開したディレクトリのなかにある ``mVMCconfig.sh`` スクリプトを次のように実行してください。
(物性研システムB”sekirei”の場合)

.. code-block:: bash

   $ bash mVMCconfig.sh sekirei

これによりコンパイル環境設定ファイル ``make.sys`` が ``src/`` ディレクトリに作られます。
``mVMCconfig.sh`` の引数は次のものに対応しています。

-  ``sekirei`` : 物性研究所システムB "sekirei"

-  ``kei`` : 京コンピューターおよび物性研究所システムC "maki"(FX10)

-  ``intel-openmpi`` : Intel コンパイラ + OpenMPI

-  ``intel-mpich`` : Intelコンパイラ + MPICH2

-  ``intel-intelmpi`` : Intelコンパイラ + IntelMPI

-  ``gcc-mpich-mkl`` : GCC + MPICH + MKL

-  ``gcc-openmpi`` : GCC + OpenMPI

``make.sys`` の中身は次のようになっています(物性研システムB
"sekirei"の場合)。

.. code-block:: makefile

   CC = mpicc
   F90 = mpif90
   CFLAGS = -O3 -no-prec-div -xHost -qopenmp -Wno-unknown-pragmas
   FFLAGS = -O3 -implicitnone -xHost
   LIBS = -L $(MKLROOT)/lib/intel64 -lmkl_scalapack_lp64 -lmkl_intel_lp64 \
           -lmkl_intel_thread -lmkl_core -lmkl_blacs_sgimpt_lp64 -lpthread -lm
   SFMTFLAGS = -no-ansi-alias -DHAVE_SSE2

となります。それぞれのマクロ(変数)の説明は次のとおりです。

-  ``CC`` : C コンパイラー( ``mpicc``, ``mpifccpx`` など)

-  ``F90`` : fortran コンパイラー( ``ifort``, ``frtpx`` など)

-  ``LIBS`` : リンカーオプション。

-  ``CFLAGS`` : C コンパイルオプション。

-  ``FFLAGS`` : fortran コンパイルオプション。

これでコンパイルのための準備が整います。その後

.. code-block:: bash

    $ make mvmc

とすることで実行可能ファイル ``vmc.out`` 、 ``vmcdry.out`` が ``src/内に`` 生成されるので、
このディレクトリにパスを通すか、
パスの通っている場所にシンボリックリンクを作ってください。

実行ファイルにパスを通す時には、次のようにします。
``$ export PATH=${PATH}:`` ``/src/``
この設定を常に残すには、例えばログインシェルが ``bash`` の場合には
``~/.bashrc`` ファイルに上記のコマンドを記載します。

cmakeを使う方法
~~~~~~~~~~~~~~~

mVMCを展開したディレクトリのパスを$PathTomvmc、ビルドディレクトリを$HOME/build/mvmc
(任意の場所を指定可能)とした場合に、

.. code-block:: bash

   cd $HOME/build/mvmc
   cmake -DCONFIG=gcc $PathTomvmc
   make

でコンパイルすることができます。コンパイル後、 ``$HOME/build/mvmc``
直下に ``src`` フォルダが作成され、実行ファイルである ``vmc.out`` がそのフォルダ内に作成されます。

なお、上の例ではgccコンパイラを前提としたコンパイルになっていますが、

-  ``sekirei`` : 物性研究所システムB "sekirei"

-  ``fujitsu`` : 富士通コンパイラ

-  ``intel`` : intelコンパイラ + Linux PC

-  ``gcc`` : GCC + Linux PC

のオプションが用意されています。以下、mVMCを展開したディレクトリでビルドする例を示します(intelコンパイラの場合)。

.. code-block:: bash

   mkdir ./build
   cd ./build
   cmake -DCONFIG=intel ../
   make

実行後、 ``build/`` フォルダ直下に ``src/`` フォルダが作成され、 ``vmc.out`` が ``src/`` フォルダ内に作成されます。
また、LAPACKに代わりScaLAPACKを計算に使用することが可能です。その場合には、

::

    -DUSE_SCALAPACK=ON -DSCALAPACK_LIBRARIES="xxx"

をcmakeをする際に付け加えてください(xxxにはScaLAPACKを利用するためのライブラリ一式を指定します)。
なお、コンパイラを変更しコンパイルし直したい場合には、都度buildフォルダごと削除を行った上で、新規に上記作業を行うことをお薦めします。

.. note::

   sekirei で cmake を利用するには

   .. code-block:: bash

      $ source /home/issp/materiapps/tool/env.sh

   をあらかじめ実行する必要があります。

   またScaLAPACK を利用するには

   .. code-block:: bash

      cmake -DCONFIG=sekirei ../ -DUSE_SCALAPACK=ON

   を行うと、デフォルトで

   ::

      -DSCALAPACK_LIBRARIES="\${MKLROOT}/lib/intel64 -lmkl_scalapack_lp64 
      -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core
      -lmkl_blacs_sgimpt_lp64"

が指定されます。ライブラリへのパスが異なる場合には、 ``-DSCALAPACK_LIBRARIES`` を適宜変更してください。

ディレクトリ構成
----------------

mVMC-xxx.gzを解凍後に構成されるディレクトリ構成を以下に示します。::

   .
   |-- CMakeLists.txt
   |-- COPYING
   |-- README.md
   |-- config/
   |   |-- fujitsu.cmake
   |   |-- gcc.cmake
   |   |-- intel.cmake
   |   `-- sekirei.cmake
   |-- dist.sh*
   |-- doc/
   |   |-- CMakeLists.txt
   |   |-- bib/
   |   |-- en/
   |   |-- figs/
   |   |-- ja/
   |   |-- package/
   |   `-- userguide.html
   |-- mVMCconfig.sh
   |-- samples/
   |   |-- Standard/
   |   |   |-- Hubbard/
   |   |   |-- Kondo/
   |   |   `-- Spin/
   |   `-- Wannier/
   |       |-- Sr2CuO3/
   |       `-- Sr2VO4/
   |-- src/
   |   |-- ComplexUHF/
   |   |-- StdFace/
   |   |-- common/
   |   |-- mVMC/
   |   |-- pfapack/
   |   `-- sfmt/
   |-- test/
   |   |-- CMakeLists.txt
   |   |-- python/
   |   `-- tool/
   `-- tool/

基本的な使い方
--------------

mVMCは次の二つのいずれかのモードで動作します。

-  エキスパートモード

   mVMCでは一般的な格子フェルミオン/スピン系に対応しており、
   各サイト毎にホッピング等を別々に指定することが出来ます。
   これにより計算対象の柔軟な指定が可能となりますが、用意する入力ファイルは多く、
   計算のセットアップは比較的煩雑になります。

-  スタンダードモード

   典型的なモデル(正方格子上のHeisenbergモデルなど)では、
   計算するセルのサイズや共通の相互作用項の大きさなど少数のパラメーターのみを入力して
   エキスパートモード用の入力ファイルを自動生成し、計算をすることが出来ます。
   計算対象はエキスパートモードに比べて限られますが、比較的容易に計算をセットアップすることが出来ます。
   また、エキスパートモード用の入力ファイルを自動生成した後、計算をする前にそれらを手動で編集して
   より広範なモデルに対応させることも可能です。

これらのモードを用いて次のように計算を行います。

#. 計算用ディレクトリの作成

   計算シナリオ名を記載したディレクトリを作成します。

#. スタンダードモードの入力ファイルの作成

   あらかじめ用意されたいくつかのモデル(HeisenbergモデルやHubbardモデル)や格子(正方格子など)を指定し、
   それらに対するいくつかのパラメーター(最近接 :math:`\cdot` 次近接スピン結合やオンサイトクーロン積分など)を設定します。
   各ファイルは :ref:`HowToStandard` に従い記載してください。

#. 実行

   作成した入力ファイル名を引数として ``vmc.out`` を実行します。
   このとき入力ファイル名の前にオプション ``-s`` を付けます。

   .. code-block:: bash
                
      $ mpiexec -np プロセス数 パス/vmc.out -s 入力ファイル

   ワークステーションやスパコン等でキューイングシステムを利用している場合は
   プロセス数をジョブ投入コマンドの引数として与える場合があります。
   詳しくはお使いのシステムのマニュアルをご参照ください。

#. 途中経過

   計算実行の経過についてカレントディレクトリ直下の ``output/``
   ディレクトリ(無ければ作られる)にログファイルが出力されます。
   出力されるファイルの詳細に関しては :ref:`outputfile` を参考にしてください。

#. 最終結果

   計算が正常終了した場合、 計算モードに従い ``output/``
   ディレクトリに計算結果ファイルが出力されます。
   出力されるファイルの詳細に関しては :ref:`outputfile` を参考にしてください。

#. エキスパートモードの入力ファイルの作成と実行

   上の例ではエキスパートモードのファイルを自動生成した後そのまま計算を開始していますが、
   エキスパートモードのファイルの生成のみを行う場合には ``vmcdry.out`` を実行します。
   MPIは使用しません。

   .. code-block:: bash
                
      $ パス/vmcdry.out 入力ファイル

   このとき生成されたファイルを必要に応じて手動で編集したのち、 ``-e``
   というオプションの後に
   ``namelist/def`` というファイルを引数として ``vmcd.out``
   を実行します。

   .. code-block:: bash
                
      $ mpiexec -np プロセス数 パス/vmc.out -e namelist.def

   以降はスタンダードモードと同様です。

**OpenMPスレッド数の指定**

実行時のOpenMPのスレッド数を指定する場合は、
``vmc.out`` を実行する前に以下の様にしてください(16スレッドの場合)。

.. code-block:: bash

   export OMP_NUM_THREADS=16

バージョン番号の確認
~~~~~~~~~~~~~~~~~~~~

次のように ``-v`` オプションをつけて ``vmc.out``,
``vmcdry.out`` を実行すると, バージョン番号を標準出力した後終了します。

.. code-block:: bash

    $ パス/vmcdry.out -v
    $ パス/vmc.out -v
