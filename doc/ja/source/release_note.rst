リリースノート
==============

未リリース
----------

Power-Lanczos estimatorの選択
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Power-Lanczos inputはcorrected estimatorをデフォルトで選択するようになりました。
``NLanczosMode`` を設定し、新しい ``NLanczosEstimatorMode`` を省略した既存inputも
対象です。``NLanczosMode=1`` はfull-support、scale-normalized corrected
energy/variance経路を独立なcoefficient/final chainで実行し、compactな安定化JSONを
1ファイル出力します。これは数値安定化のevidenceであり、release統計認証では
ありません。corrected ``NLanczosMode=2`` は引き続きscope外で、メモリ確保前に
拒否します。追加observableのproduction enableは変更しません。

corrected launcherは :ref:`HowToExpert` に記載したsource/input/binary/environment/
seedの正確なidentity変数を指定する必要があります。固定default policyは
4096+4096のscale pilot、各scored stageのwarmup 4096・saved 16384、32 blocksを
使用し、sample-level traceを永続保存しません。

従来のbiased base-support estimatorを一時的に再現するには、ModPara fileへ次の行を
追加してください::

   NLanczosEstimatorMode 1

このlegacy経路はwarningを出します。出力はdiagnostic専用であり、corrected release
resultとして扱わないでください。他の新しいP6 estimator controlはすべて0にする
必要があります。キーワードの完全な契約は :ref:`HowToExpert` を参照してください。
