.. _tutorial:

チュートリアル
==============

このチュートリアルは ``sample/Standard/Spin/HeisenbergSquare/``
にあるインプットファイルを用いて行う.

HPhi/vmc.out の実行
-------------------

- :math:`{\mathcal H}\Phi` の場合

  基底状態および相関関数の計算を行う.
  
  .. code-block:: bash

     $ ../../../../src/HPhi -s StdFace.def

- mVMC の場合

  変分波動関数の最適化を行う.
  
  .. code-block:: bash

     $ ../../../../src/vmc.out -s StdFace.def

  相関関数を計算するために, ``StdFace.def`` に以下の行を付け加える.

  ::

     NVMCCalMode = 1

  相関関数を計算する.
  
  .. code-block:: bash

     $ ../../../../src/vmc.out -s StdFace.def zqp_opt.dat
         
これにより, カレントディレクトリの ``output/`` 以下に
1体および2体の相関関数が出力される.

関連するファイル

- StdFace.def (mVMC/:math:`{\mathcal H}\Phi` のマニュアル参照)
- zqp_opt.dat (mVMCのマニュアル参照)
- greenone.def (:ref:`greenindex`)
- greentwo.def (:ref:`greenindex`)

相関関数のフーリエ変換
----------------------

ユーティリティプログラム ``fourier`` を使って,
相関関数をフーリエ変関する.

.. code-block:: bash

   $ ../../../../tool/fourier namelist.def geometry.dat
     
これにより, カレントディレクトリの ``output/`` 以下に
フーリエ変換された相関関数が出力される.

関連するファイル

- output/zvo_cisajs_001.dat (:ref:`zvocisajs`)
- output/zvo_cisajs.dat (:ref:`zvocisajs`)
- output/zvo_cisajscktalt_001.dat (:ref:`zvocisajs`)
- output/zvo_cisajscktalt.dat (:ref:`zvocisajs`)
- geometry.dat (:ref:`geometry`)
- output/zvo_corr.dat (:ref:`zvocorr`)

相関関数のプロット
------------------

ユーティリティプログラム ``corplot`` を使って,
相関関数を :math:`k` 空間でプロットする.

.. code-block:: bash

   $ ../../../../tool/corplot output/zvo_corr.dat

関連するファイル

- kpoint.dat (:ref:`kpoint`)
- correlation.gp (:ref:`gnuplot`)
- correlation.dat (:ref:`correlation`)
