.. _tutorialwannier:

チュートリアル
==============

このチュートリアルでは Sr\ :sub:`2`\ VO\ :sub:`4`
を2次元3軌道Hubbardモデルにダウンフォールドして,
それをHPhi/mVMCで計算する.
DFT計算はQuantumESPRESSOで行う.

電荷密度のSCF計算
-----------------

:download:`scf.in <../../../samples/Wannier/Sr2VO4/scf.in>`

.. literalinclude:: ../../../samples/Wannier/Sr2VO4/scf.in

(Optional) バンド計算と描画
---------------------------

:download:`band.in <../../../samples/Wannier/Sr2VO4/band.in>`

.. literalinclude:: ../../../samples/Wannier/Sr2VO4/band.in
                    
:download:`bands.in <../../../samples/Wannier/Sr2VO4/bands.in>`

.. literalinclude:: ../../../samples/Wannier/Sr2VO4/bands.in
                    
Kohn-Sham軌道の計算
-------------------

:download:`nscf.in <../../../samples/Wannier/Sr2VO4/nscf.in>`

.. literalinclude:: ../../../samples/Wannier/Sr2VO4/nscf.in

Wannier関数, 誘電関数, 有効相互作用の計算
-----------------------------------------

:download:`respack.in <../../../samples/Wannier/Sr2VO4/respack.in>`

.. literalinclude:: ../../../samples/Wannier/Sr2VO4/respack.in

HPhi/mVMCによるモデル計算
-------------------------

:download:`respack.in <../../../samples/Wannier/Sr2VO4/stan.in>`

.. literalinclude:: ../../../samples/Wannier/Sr2VO4/stan.in
