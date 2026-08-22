リリースノート
==============

未リリース
----------

Power-Lanczos estimatorの選択
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Power-Lanczos inputはcorrected estimatorをデフォルトで選択するようになりました。
``NLanczosMode`` を1または2に設定し、新しい ``NLanczosEstimatorMode`` を省略した
既存inputも対象です。段階導入中はcorrected production経路をまだ有効化していない
ため、該当inputはavailability diagnosticを出し、exit status 20でfail-fastします。
``NLanczosMode=2`` ではP6 observable censusより前に終了するため、census errorや
resource-limit errorがcorrected経路の未提供を覆い隠すことはありません。

従来のbiased base-support estimatorを一時的に再現するには、ModPara fileへ次の行を
追加してください::

   NLanczosEstimatorMode 1

このlegacy経路はwarningを出します。出力はdiagnostic専用であり、corrected release
resultとして扱わないでください。他の新しいP6 estimator controlはすべて0にする
必要があります。キーワードの完全な契約は :ref:`HowToExpert` を参照してください。
