src以下にコードが入っています。
コンパイルはsrc/com.shで出来るようにしてあります。

Test_Hubbard_4by4以下にHubbard model 4by4, U/t=4
でのサンプルが入っています。

UHF自体はVMCと同じように相互作用パラメータ[zinterall.defは未対応]を用意
すれば、

UHF UHFnamelist.def

で動きます。
VMCとの違いは
UHFmodpara.defとzinitial.defです。

1. UHFmodpara.defはVMCの
zmodpara.defと大体同じです。
以下の二つがUHF独特ののパラメータです。

mix-> linear mixingのパラメータmix=1で完全に
新しいGreen関数に置き換える
eps-> 収束判定条件。残渣が10^(-eps)で計算を打ち切る。

2.zinitial.def
Green functionの初期値を与えます。

site_i site_j spin G_{ij}

の形で指定しています。
値を指定しないグリーン関数には０が入るようにしています。

3.出力ファイル。
zvo_result.dat -> energyと粒子数
zvo_check.dat-> 収束過程のエネルギー、残渣
zvo_UHF_cisajs.dat -> 収束したGreen関数
zvo_eigen.dat -> 収束したハミルトニアンの固有値
zvo_gap.dat -> エネルギーギャップ
zvo_fij.dat -> slater determinantから生成したfij

4. X.shが一連の動きをまとめたスクリプトです。
パールスクリプト内で
input.txtという、格子のジオメトリーを
指定したファイルを参照しています。

./IniGreen.pl -> 反強磁性のzinitial.defを生成

./UHF UHFnamelist.def -> UHFを実行

./UHF_RealSpin.pl ->
zvo_UHF_cisajs.datから実空間の電荷・スピンを作成　( Result_UHF_local_0.dat)
+zvo_UHF_cisajs.datからFourier変換したの電荷・スピンを作成
(Result_UHF_Q_0.dat)

./MakeFij.pl -> zvo_fij.dat + VMCの変分パラメータの指定ファイルから
VMCのinitialファイル VMC_zinitial.defを生成。
