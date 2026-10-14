<script type="text/javascript" async src="https://cdnjs.cloudflare.com/ajax/libs/mathjax/2.7.7/MathJax.js?config=TeX-MML-AM_CHTML">
</script>
<script type='text/x-mathjax-config'>
	MathJax.Hub.Config({
		tex2jax: {
			inlineMath: [['$', '$']],
			displayMath: [ ['$$','$$'],["\\[","\\]"] ]
		}
	});
</script>

# ActFEM v1.4 マニュアル
## 有限要素法による ACTIVE 順計算・インバージョン コード

#### 2026年5月30日 現在

南 拓人（Takuto MINAMI）
tminami@port.kobe-u.ac.jp
（神戸大学 理学部）

### 目次

1. [必要な環境](#Required_environment)
2. [座標系](#Coordinate_system)
3. [順計算サンプルの実行](#Run_the_sample)
4. [インバージョンサンプルの実行](#Run_samples_for_inversions)
5. [ソースフォルダ構成](#Source_folder_tree)
6. [支配方程式と離散化](#Governing_equation_and_discretization)
7. [ファイル形式](#File_formats)
8. [メッシュ生成チュートリアル](#Mesh_generation_tutorial)
9. [順計算・インバージョン用比抵抗モデルの作成](#Generation_of_resistivity_model)
10. [バージョン履歴](#Version_notes)
11. [参考文献](#References)

- [付録A：支配方程式の離散化形式の導出](#AppendixA)
- [付録B：ActFEM 計算のための支配方程式のスケーリング](#AppendixB)



## 1 必要な環境 <a id="Required_environment"></a>
ActFEM には、OpenMP 用の Intel コンパイラ（ifx または ifort、mpiifx または mpiifort）と Intel MPI が必要です。順計算コードには OpenMP のみのシリアルモードがあります。一方、インバージョンには OpenMP と MPI のハイブリッド並列版のみが用意されています。必要な環境は以下のとおりです。

### 順計算・インバージョン
- Intel コンパイラ（連立一次方程式は Intel MKL ライブラリの PARDISO で解きます）
- gmsh
  - Mac → `brew install gmsh`
  - Ubuntu → `sudo apt install gmsh`

### 可視化
- GMT6
  - Mac → `brew install gmt`
  - Ubuntu → `sudo apt install gmt`
- Gnuplot
  - Mac → `brew install gnuplot`
  - Ubuntu → `sudo apt install gnuplot`

メッシュ生成（src/src_mesh）と後処理（src/src_post）には gfortran も使用できます。

## 2 座標系 <a id="Coordinate_system"></a>
ActFEM では以下の座標系を用います。

- X : 東向き
- Y : 北向き
- Z : 上向き

## 3 順計算サンプルの実行 <a id="Run_the_sample"></a>
順計算のサンプルフォルダ ActFEM/Fwd_test があります。ここで、ActFEM による ACTIVE と MT の順計算を確認できます。

### 3.1 Step 1: メッシュ生成（地形なし）
```bash
$ cd ActFEM/Fwd_test/mesh_flat
$ chmod +x meshgen.sh
$ ./meshgen.sh
```
meshgen.sh は、まず ActFEM/src/src_mesh のソースを gfortran でコンパイルし、gmsh を使ってプログラムを実行します。gmsh が正しくインストールされていれば、ターミナルに次のような出力が表示されます。

<img src="./images/v1.4_image1.png" width="50%">

次に、gmsh ビューアでメッシュを確認してください。
```bash
$ gmsh nakadake3d.msh
```

<img src="./images/v1.4_image2.png" width="30%"> <img src="./images/v1.4_image3.png" width="27%"> <img src="./images/v1.4_image4.png" width="36%">

図. Tools > Options > Mesh で「3D element faces」にチェックを入れると中央の図のように表示されます。「Tools > Visibility > Volume 2 を選択 > Apply」で陸側の要素のみが表示され、観測点付近の高解像度領域を確認できます。

### 3.2 Step 2: ACTIVE の解析解の準備
順計算の精度を評価するため、順計算の前に解析解を用意しておきます。
```bash
$ cd ActFEM/Fwd_test/active_ana/
$ chmod +x run.sh
$ run.sh
```
デフォルトでは、run.sh は ifx で src/ana のソースをコンパイルします。run.sh が正常に終了すると、次のように表示されます。

<img src="./images/v1.4_image5.png" width="52%">

### 3.3 Step 3: ACTIVE 順計算
メッシュを生成したら、ACTIVE と MT の順計算を実行できます。
```bash
$ cd ActFEM/Fwd_test/active
$ chmod +x run.sh
$ ./run.sh
```
ACTIVE 順計算が正常に終了すると、次のように表示されます。

<img src="./images/v1.4_image6.png" width="43%">

### 3.4 Step 4: ACTIVE 順計算結果と解析解の比較
ACTIVE の計算が終わったら、Fwd_test/active_ana/result の解析解と比較して順計算の精度を確認できます。

Fwd_test/active/ で
```bash
$ chmod +x compana.sh
$ ./compana.sh
```
以下の図が表示されれば、ACTIVE 順計算と解析解の計算はどちらも成功しています。

<img src="./images/v1.4_image7.png" width="72%">

図. シンボルは Fwd_test/active/ で計算した ACTIVE 応答、線は Fwd_test/active_ana/ で計算した解析解です。（上：振幅 [nT/A]、下：位相 [deg]）

### 3.5 Step 5: 一様比抵抗構造での MT 順計算（MPI なし）
MT 計算には、シリアル（OpenMP）モードとハイブリッド（OpenMP + MPI）モードがあります。まずシリアルモードを確認します。
```bash
$ cd ActFEM/Fwd_test/MTfwd_test
$ chmod +x run.sh
$ ./run_homo.sh
```
計算が成功すると、次のような出力が得られます。

<img src="./images/v1.4_image8.png" width="63%">

### 3.6 Step 6: MT 順計算の精度確認
計算が終わったら、MT 計算の結果を確認してください。
```bash
$ chmod +x plot_rho.sh
$ ./plot_rho.sh
```

<img src="./images/v1.4_image9.png" width="56%">

図のとおり、ActFEM の MT 計算は 1〜100 Hz の周波数範囲で、見かけ比抵抗 100 Ohm.m・位相 45 度を数値誤差 1% 未満で再現できています。

### 3.7 Step 7: 水平方向に不均質な構造での MT 順計算（MPI なし）
MT 順計算ではティッパーも計算できます。水平方向に不均質な構造でその動作を確認できます。まず不均質構造を作成します。ActFEM/Fwd_test/MTfwd_test/ で
```bash
$ chmod +x change_model2cond.sh
$ ./change_model2cond.sh
```
次のような出力とともに **cond_lateral.msh** が生成されます。

<img src="./images/v1.4_image10.png" width="64%">

生成された比抵抗構造を確認するには、以下のコマンドを実行します。
```bash
$ chmod +x view_slcie.sh
$ view_slcie.sh
```

<img src="./images/v1.4_image11.png" width="59%">

東西・南北断面に見られるように、生成された比抵抗構造は x=0 を境に、西側の 10 Ohm.m から東側の 100 Ohm.m へと水平方向に不連続に変化しています。

次に、この構造で MT 順計算を実行します。
```bash
$ chmod +x run_lateral.sh
$ ./run_lateral.sh
```
計算が終わったら、ティッパーを確認します。
```bash
$ chmod +x plot_tip.sh
$ ./plot_tip.sh
```

<img src="./images/v1.4_image12.png" width="51%">

図のとおり、西側の比抵抗が低いため、インダクションベクトルの実部は西を向いています。

### 3.8 Step 8: MPI + OpenMP による MT 順計算
現在、MPI で並列化されているのは MT 順計算のみです。Fwd_test/MTfwd_test と同じ計算を、Fwd_test/MTfwd_test_mpi/ で MPI を使って実行できます。
```bash
$ cd Fwd_test/MTfwd_test_mpi/
$ chmod +x run_homo.sh
$ ./run_homo.sh
```
計算が成功すると、次のような出力が得られます。

<img src="./images/v1.4_image13.png" width="55%">

ip=1（ノード1）が freq=3.0Hz を、ip=3（ノード3）が freq=30Hz を計算している、というように、周波数ごとの計算が各 MPI プロセスに割り当てられます。そのため、シリアル版よりも大幅に高速です。MTfwd_test/result_homo/ と MTfwd_test_mpi/result_homo/ の結果が一致することを確認してください。

## 4 インバージョンサンプルの実行 <a id="Run_samples_for_inversions"></a>
ActFEM には、ACTIVE 2件、MT 2件、ACTIVE+MT 同時インバージョン 2件の、計6件のインバージョンサンプルがあります。インバージョンを実行するには、まず ACTIVE と MT のサンプルデータセットを用意します。

### 4.1 Step 1: メッシュの準備
サンプルインバージョンは阿蘇の実地形を用いて行います。まずメッシュを生成しましょう。
```bash
$ cd ActFEM/Joint_test/mesh_joint
$ chmod +x meshgen.sh
$ ./meshgen.sh
```
メッシュを生成したら、gmsh ビューアで確認してください。
```bash
$ gmth nakadake3d.msh
```

<img src="./images/v1.4_image14.png" width="34%"> <img src="./images/v1.4_image15.png" width="36%">

図. 阿蘇火山（九州中部）中岳第一火口周辺の四面体メッシュ。

### 4.2 Step 2: 真の比抵抗構造の準備
インバージョンのチュートリアルでは、目標とする比抵抗構造を仮定する必要があります。ActFEM では ActFEM/Joint_test/structure に次の目標構造を用意しています。
```bash
$ cd ActFEM/Joint_test/structure
$ ./cond_G3.sh
$ ./cond_homo.sh
```
成功すると **cond_G3.msh** と **cond_homo.msh** が生成されます。cond_G3.msh は下図の目標構造、cond_homo.msh は以下のインバージョンで初期モデルおよび参照モデルとして使う一様構造です。cond_\*.msh は元の mesh_joint/nakadake3d.msh と一対一に対応している点に注意してください。mesh_joint/nakadake3d.msh を更新した場合は、cond_G3.msh と cond_homo.msh も作り直してください。

cond_G3.msh の異常体は次のコマンドで確認できます。
```bash
$ ./view_slice.sh
```

<img src="./images/v1.4_image16.png" width="79%">

図のとおり、異常体は火口底の直下に広がる直方体です。

### 4.3 Step 3: インバージョン用 ACTIVE サンプルデータの作成
真のモデルができたので、順計算によってインバージョン用のサンプルデータを作成します。まず ACTIVE データを次のように用意します。
```bash
$ cd ActFEM/Joint_test/fwd_active/
$ chmod +x run.sh
$ ./run.sh
```
順計算が終わったら、ActFEM/Joint_test/fwd_active で誤差付きのデータを作成します。
```bash
$ chmod + generr.sh
$ ./generr.sh
```
その結果、Joint_test/fwd_active/data/ に `A0*_S?_amp.dat` と `A0*_S?_pha.dat` が生成されます。`A0*` は active.ctl に記載された観測点名、`S?` は active.ctl に記載された電気ダイポール送信源の名前です。

### 4.4 Step 4: インバージョン用 MT サンプルデータの作成
次に、ActFEM/Joint_test/fwd_3DMT/ で MT データを用意します。
```bash
$ cd ActFEM/Joint_test/fwd_3DMT/
$ chmod +x run.sh
$ ./run.sh
```
これにより **result/MT1_MT.dat** などが生成されます。誤差付きデータは次のコマンドで作成します。
```bash
$ chmod +x generr.sh
$ ./generr.sh
```
すると **mt_err/MT1_MT_imp.err** などの誤差ファイルが作成されます。

### 4.5 Step 5: サンプルインバージョンの実行
以上で、すべてのインバージョンの準備が整いました。どのサンプルも手順は同じで、**run.sh** を実行し、**plot_rms.sh** と **view_slice.sh** で結果を確認します。たとえば Joint_test/inv_act/ では次のようにします。
```bash
$ cd ActFEM/Joint_test/inv_act/
$ ./run.sh
```
run.sh 内の OpenMP スレッド数と MPI プロセス数は、お使いの環境に合わせて編集してください。

以下では、得られたインバージョン結果のみを示します。

#### Joint_test/inv_act の結果
```bash
$ plot_rms.sh
$ ./view_slice.sh 06
```
ここで 06 は、収束したモデル model06 を表示することを意味します。

<img src="./images/v1.4_image17.png" width="28%"> <img src="./images/v1.4_image18.png" width="66%">

#### Joint_test/inv_mt/ の結果
```bash
$ ./plot_rms.sh
$ ./view_slice.sh 09
```

<img src="./images/v1.4_image19.png" width="31%"> <img src="./images/v1.4_image20.png" width="63%">

#### Joint_test/inv_joint/ の結果
```bash
$ ./plot_rms.sh
$ ./view_slice.sh 08
```

<img src="./images/v1.4_image21.png" width="30%"> <img src="./images/v1.4_image22.png" width="64%">

> **以下は icombine = -1 の結果です。この設定では、インバージョンで内側・外側のモデルパラメータ化を使用できます。詳細は次節で説明します。**

#### Joint_test/inv_act_icomb-1 の結果
```bash
$ plot_rms.sh
$ ./view_slice.sh 06
```

<img src="./images/v1.4_image23.png" width="28%"> <img src="./images/v1.4_image24.png" width="66%">

#### Joint_test/inv_mt_icomb-1 の結果
```bash
$ plot_rms.sh
$ ./view_slice.sh 10
```

<img src="./images/v1.4_image25.png" width="31%"> <img src="./images/v1.4_image26.png" width="63%">

#### Joint_test/inv_joint_icomb-1 の結果

<img src="./images/v1.4_image27.png" width="29%"> <img src="./images/v1.4_image28.png" width="66%">

#### Joint_test/inv_mt_tip の結果
mt impedanceのerrorfloorは0.1、tipperのerrorfloorは0.1に設定し、αはpower -0.3で減少させた。-0.3が大きいと、構造の変動が大きくなり収束しない。

<img src="./images/v1.4_image29.png" width="33%"> <img src="./images/v1.4_image30.png" width="62%">

## 5 ソースフォルダ構成 <a id="Source_folder_tree"></a>

| フォルダ             | 内容                              | テストフォルダ      |
|----------------------|-----------------------------------|---------------------|
| src/solver/          | シリアル版 ACTIVE 順計算ソルバー  | Fwd_test/active     |
| src/solver_mpi/      | MPI 版 ACTIVE 順計算ソルバー      | Fwd_test/active_mpi |
| src/src_3DMT/        | シリアル版 3次元 MT ソルバー      |                     |
| src/src_3DMT_mpi/    | MPI 版 3次元 MT ソルバー          |                     |
| src/src_inv_joint/   | ACTIVE/MT/同時インバージョンコード |                     |

## 6 支配方程式と離散化 <a id="Governing_equation_and_discretization"></a>
ActFEM の順計算・インバージョンの詳細は Minami et al. (2018) を参照してください。ベクトルポテンシャル $\mathbf{A}$ に関する誘導方程式

$$\nabla \times \left( \nabla \times \mathbf{A} \right) + i\omega\mu\widehat{\sigma}\mathbf{A} = \mu\mathbf{i}_{s},$$

を辺要素有限要素法で解きます。ここで $\widehat{\sigma} = \sigma + i\omega\epsilon$ であり、$\epsilon$ は誘電率、$\omega$ は角周波数、$\mu$ は透磁率、$\mathbf{i}_{s}$ は送信源の電流密度です。なお、$\mathbf{E} = i\omega\mathbf{A}$ の関係を用いれば、この式は電場 $\mathbf{E}$ に関する誘導方程式

$$\nabla \times \nabla \times \mathbf{E} + i\omega\mu\widehat{\sigma}\mathbf{E} = i\omega\mu\mathbf{i}_{s}$$

に容易に書き換えられます。

磁場は、$\mathbf{A}$ について方程式を解いた後に $\mathbf{B} = \nabla \times \mathbf{A}$ から求めます。

## 7 ファイル形式 <a id="File_formats"></a>

### 7.1 制御ファイル
順計算は制御ファイル “nakaya.ctl” で制御します（括弧内は説明）。

```
## lines starting with "##" work as lines for comments ! 2020.09.28
## this file is forward control file
##-----10!-------20!
itopofile 0 or 1 !1 (0:地形ファイルを使わずメッシュ生成、1:地形ファイルを使う。)
# of topofile !1 (使う地形ファイルの数。以下、狭い領域のファイルから、広い領域のファイルの順で、入力する。)
topofile !../topo/topo127_134_29_36.xyz (地形ファイル名)
lon lat shift !0.0 0.0 (lon latの実際の値からのshift値)
mesh file !../mesh_aso_A04/nakadake3d.msh (順計算で使う3次元メッシュファイル)
2d triangle z file !../mesh_aso_A04/nakadake2dz.msh (地形情報付きの2次元メッシュ)
local line file !../mesh_aso_A04/lineinfo.dat (辺情報ファイル)
output folder !./result_fwd/
header2d (a50) !nakadake2d
header3d (a50) !nakadake3d
# of frequency !2 (計算する周波数の数)
Frequency [Hz] !1.d0 (1番目の周波数)
Frequency [Hz] !3.d0 (2番目の周波数)
##Frequency [Hz] !7.d0 (“##”で始まる行はコメントアウト)
##Frequency [Hz] !11.d0 (“##”で始まる行はコメントアウト)
##Frequency [Hz] !21.d0 (“##”で始まる行はコメントアウト)
##Frequency [Hz] !41.d0 (“##”で始まる行はコメントアウト)
##Frequency [Hz] !61.d0 (“##”で始まる行はコメントアウト)
##Frequency [Hz] !99.d0 (“##”で始まる行はコメントアウト)
west bound !-1.7 (メッシュ用パラメータ。順計算では未使用)
east bound !1.7 (メッシュ用パラメータ。順計算では未使用)
south bound !-1.5 (メッシュ用パラメータ。順計算では未使用)
north bound !1.5 (メッシュ用パラメータ。順計算では未使用)
lenout [km] !50.0 (メッシュ用パラメータ。順計算では未使用)
upz in [km] (>0) !1.3 (メッシュ用パラメータ。順計算では未使用)
downz in [km](<0) !-1.1 (メッシュ用パラメータ。順計算では未使用)
zmax [km] !50.0 (メッシュ用パラメータ。順計算では未使用)
zmin [km] !-50.0 (メッシュ用パラメータ。順計算では未使用)
sizein [km] !0.15 (メッシュ用パラメータ。順計算では未使用)
sizebo [km] !10.0 (メッシュ用パラメータ。順計算では未使用)
sigma_obs [km] !0.4 (メッシュ用パラメータ。順計算では未使用)
A_obs [km] !0.01 (メッシュ用パラメータ。順計算では未使用)
dlen_source [km] !0.1 (メッシュ用パラメータ。順計算では未使用)
sigma_src [km] !0.3 (メッシュ用パラメータ。順計算では未使用)
A_src [km] !0.005 (メッシュ用パラメータ。順計算では未使用)
# of observatory !4 (磁場・電場を計算する観測点の数)
lonlat(1),xyz (2) !1 (1: 経度緯度で入力、2: xyz で入力)
UTM ZONE |52S (経度緯度入力時の UTM ゾーン)
lonlatorigin !131.084782 32.884882 (メッシュ原点の経度緯度)
1 Name !A02 (観測点名 1)
1 xyz !131.083411 32.886706 -0.001 (経度 緯度 高度。高度は地表面からの値)
2 Name !A04 (観測点 2)
2 xyz,sigma,A[km] !131.081939 32.884808 -0.001 (経度 緯度 高度。高度は地表面からの値)
3 Name !A01 (観測点 3)
3 xyz,sigma,A[km] !131.083367 32.882725 -0.001 (経度 緯度 高度。高度は地表面からの値)
4 Name !A03 (観測点 4)
4 xyz,sigma,A[km] !131.086847 32.881981 -0.001 (経度 緯度 高度。高度は地表面からの値)
ixyflg 0:no,1:surfv!0 (地表面マップのフラグ)
# of sources !1 (電気ダイポール送信源の数)
Source Name !S1 (1番目の送信源の名前)
source start point !131.0784333 32.8908028 -0.001 (ダイポールの始点)
source end point !131.0814639 32.8912333 -0.001 (ダイポールの終点)
Elcetric current[A]! 1.0 (電気ダイポールに流す電流値)
sigma_air [S/m] !1.e-8 (空気の電気伝導度)
condflag 0:home,1: !0 (地下の電気伝導度フラグ。0: 一様大地、1: ファイル入力)
nvolume !1 (フラグが0のときのボリューム数)
cond !0.01 (地表面下の1番目のボリュームの電気伝導度)
```

## 8 メッシュ生成チュートリアル <a id="Mesh_generation_tutorial"></a>
有限要素計算のためのメッシュ生成は、順計算・インバージョンの準備の中で最も手間のかかる作業の一つです。ここでは、ActFEM の順計算・インバージョンでメッシュ生成コードを使うための簡単なチュートリアルを示します。チュートリアルは ActFEMv1.0/mesh_demo_Aso/ フォルダで始めてください。

ActFEM/mesh_demo_Aso/ では、meshgen.sh と meshgen_add.sh を使って2種類の3次元メッシュを作成できます。両者の違いは、使用する制御ファイル（aso.ctl と aso_add.ctl）だけです。ActFEM/mesh_demo_Aso/ で次のいずれかを実行してください。
```bash
$ ./meshgen.sh
```
または
```bash
$ ./meshgen_add.sh
```
meshgen.sh と meshgen_add.sh の出力ファイル名はすべて同じである点に注意してください。meshgen.sh または meshgen_add.sh で生成した出力ファイルをすべて削除したい場合は、ActFEM/mesh_demo_Aso/ で
```bash
$ ./clean.sh
```
を実行してください。2つのスクリプトの違い、すなわち2つの制御ファイルの違いは観測点の数です。aso_add.ctl にのみ観測点 B02 が追加されており、観測点数は aso.ctl の7点に対して aso_add.ctl では8点になります。違いを下図に示します。

<img src="./images/v1.4_image31.png" width="50%"> <img src="./images/v1.4_image32.png" width="46%">

図. メッシュ細分化のために aso.ctl（左）と aso_add.ctl（右）で設定した観測点。aso_add.ctl では B02 が追加されています。

meshgen.sh または meshgen_add.sh の実行に成功し、ActFEM の順計算・インバージョン用の3次元メッシュが生成されると、次のファイルが得られます。

<img src="./images/v1.4_image33.png" width="30%">

各出力ファイルがどの段階で生成されるかは meshgen.sh を参照してください。最終的な3次元四面体メッシュは nakadake3d.msh です。生成された3次元四面体メッシュの詳細は gmsh の GUI で確認できます。「Tools > Visibility > Volume 2 を選択 > Apply」で、生成された3次元メッシュの地表面が表示されます。

![](./images/v1.4_image34.png)

次に、生成された3次元メッシュの地表面を示す下の2枚の図を確認してください。meshgen_add.sh（aso_add.ctl）で生成したメッシュでは、B02 に対応する細分化領域が見られます。このように aso.ctl から aso_add.ctl へ変更するのと同様に、制御ファイルを編集することで生成される3次元メッシュを調整できます。

<img src="./images/v1.4_image35.png" width="40%"> <img src="./images/v1.4_image36.png" width="45%">

図. meshgen.sh で生成した nakadake3d.msh（左）と meshgen_add.sh で生成したもの（右）。

## 9 順計算・インバージョン用比抵抗モデルの作成 <a id="Generation_of_resistivity_model"></a>
ActFEM で順計算・インバージョンを行うには、与えられた3次元四面体メッシュに対応する比抵抗モデルを用意する必要があります。ここでは、前節で生成した3次元メッシュに対して比抵抗モデルを作成する簡単なチュートリアルを示します。このチュートリアルを始める前に、ActFEMv1.0/mesh_demo_Aso/ で meshgen.sh または meshgen_add.sh を実行しておいてください（詳細は[8章](#Mesh_generation_tutorial)を参照）。

まず ActFEMv1.0/mesh_demo_rec2tet/ に移動してください。ここでは ActFEMv1.0/mesh_demo_Aso/nakadake3d.msh をもとに比抵抗モデルを作成します。ActFEMv1.0/mesh_demo_rec2tet/ で次を実行します。
```bash
$ ./condmesh.sh
```
condmesh.sh の実行に成功すると、次のファイルが得られます。

<img src="./images/v1.4_image37.png" width="35%">

cond_init.ctl は出力ファイルではなく、condmesh.sh 実行用の制御ファイルです。cond_init.msh は ActFEM の順計算・インバージョンで使用できるファイルです。生成される cond_init.msh は、Kanda et al. (2018, JpGU) に基づき Minami et al. (2018, EPS) の初期モデルとして使われたものです。cond_init.msh には3次元メッシュの情報が含まれていないため、`gmsh cond_init.msh` だけでは表示できない点に注意してください。init.msh は mesh_demo_rec2tet/cond_init.msh と mesh_demo_Aso/nakadake3d.msh を結合したファイルです。cond_init.msh と init.msh の構造の違いは次のとおりです。

![](./images/v1.4_image38.png)

図. nakadake3d.msh、cond_init.msh、init.msh の構造

生成された比抵抗モデルは
```bash
$ gmsh init.msh
```
で表示し、下図のようにビューアのパラメータを設定してください。

![](./images/v1.4_image39.png)

図. gmsh GUI による init.msh の表示。カラースケールを反転したい場合（高比抵抗を寒色で表示したい場合）は、Option パネルの Map タブをクリックして「i」を押してください。

condmesh.sh では、ActFEMv1.0/src/src_mesh/gridmdl2mesh.exe を使って cond_init.msh を生成しています。condmesh.sh の処理の概略を下図に示します。出力する比抵抗ファイルを変更するには、cond_init.ctl を編集してください。

![](./images/v1.4_image40.png)

図. condmesh.sh の概要

## 10 バージョン履歴 <a id="Version_notes"></a>

#### Version 1.4
2026年5月30日：ACTIVE（CSEM）および MT インピーダンスとの同時インバージョンで、ティッパーインバージョンが利用可能になりました。

#### Version 1.0
最終更新日は2021年9月30日。順計算とインバージョンが利用可能です。

## 11 参考文献 <a id="References"></a>
Minami, T., Utsugi, M., Utada, H., Kagiyama, T., & Inoue, H. (2018). Temporal variation in the resistivity structure of the first Nakadake crater, Aso volcano, Japan, during the magmatic eruptions from November 2014 to May 2015, as inferred by the ACTIVE electromagnetic monitoring system. *Earth, Planets and Space*, *70*(1), 1-10.

## 付録A：支配方程式の離散化形式の導出 <a id="AppendixA"></a>
支配方程式

$$\nabla \times \nabla \times \mathbf{A} + i\omega\mu\sigma\mathbf{A} = \mu\mathbf{i}_{s}$$

に対し、ガラーキン法による有限要素法では次式を数値的に解きます。

$$\int_{\Omega}{\mathbf{w} \cdot \left( \nabla \times \nabla \times \mathbf{A} + i\omega\mu\sigma\mathbf{A} - \mu\mathbf{i}_{s} \right)dV} = 0,\  \cdot \cdot \cdot \cdot (1)$$

ここで $w$ は任意の試験関数であり、ガラーキン法では $\mathbf{w}$ に辺要素基底関数を採用します。$\Omega$ は計算領域です。

関係式 $\nabla \cdot \left( \mathbf{A} \times \mathbf{B} \right) = \mathbf{B} \cdot \left( \nabla \times \mathbf{A} \right) - \mathbf{A} \cdot (\nabla \times \mathbf{B})$ を用いると、式(1)左辺の第1項は次のように変形できます。

$$\int_{\Omega}{\mathbf{w} \cdot (\nabla \times \nabla \times A)dV} = \int_{\Omega}{\left\{ \nabla \cdot \left( (\nabla \times A) \times w \right) + (\nabla \times w) \cdot (\nabla \times A) \right\} dV}$$

$$= \int_{\Omega}{\left\{ \nabla \cdot \left( \left( \nabla \times \mathbf{A} \right) \times w \right) + (\nabla \times w) \cdot (\nabla \times A) \right\} dV}$$

$$= \int_{\partial\Omega}{\left( \left( \nabla \times \mathbf{A} \right) \times \mathbf{w} \right) \cdot \mathbf{n}dS} + \int_{\Omega}{\left\{ \left( \nabla \times \mathbf{w} \right) \cdot \left( \nabla \times \mathbf{A} \right) \right\} dV}$$

$$= \int_{\partial\Omega}{\left( \mathbf{B} \times \mathbf{w} \right) \cdot \mathbf{n}dS} + \int_{\Omega}{\left\{ \left( \nabla \times \mathbf{w} \right) \cdot \left( \nabla \times \mathbf{A} \right) \right\} dV} \cdot \cdot \cdot (2)$$

ここで $\mathbf{B} = \nabla \times \mathbf{A}$ は磁場、$\partial\Omega$ は計算領域の境界面、$\mathbf{n}$ は境界面 $\partial\Omega$ の外向き単位法線ベクトルです。$\partial\Omega$ において

$$\mathbf{n} \cdot (B \times w) = w \cdot \left( \mathbf{n} \times \mathbf{B} \right) = 0$$

とおけることを考慮します。ここで $\nabla \cdot \left( \mathbf{n} \times \mathbf{B} \right) = B \cdot \left( \nabla \times \mathbf{n} \right) - n \cdot (\nabla \times B) = - n \cdot \left( \mu\mathbf{j} \right) = 0$ は、境界面を貫く電流がないことを表します。すると式(2)は次のようになります。

$$\int_{\Omega}{\mathbf{w} \cdot \left( \nabla \times \nabla \times \mathbf{A} \right)dV = \int_{\Omega}{\left\{ \left( \nabla \times \mathbf{w} \right) \cdot \left( \nabla \times \mathbf{A} \right) \right\} dV}} \cdot \cdot \cdot (3)$$

式(3)を式(1)に代入すると、次式が得られます。

$$\int_{\Omega}{\left\{ \left( \nabla \times \mathbf{w} \right) \cdot \left( \nabla \times \mathbf{A} \right) \right\} dV} + \int_{\Omega}{i\omega\mu\sigma\mathbf{w} \cdot \mathbf{A}dV} = \mu\int_{\Omega}{\mathbf{w} \cdot \mathbf{i}_{s}dV}. \cdot \cdot \cdot (4)$$

有限要素法では、計算領域を有限個の四面体に分割し、各四面体要素について体積積分を計算します。以下では e 番目の要素 $\Omega_{e}$ の体積積分のみを考えると、式(4)は e 番目の要素について次のように書き直されます。

$$\int_{\Omega_{e}}{\left\{ \left( \nabla \times \mathbf{w}_{i} \right) \cdot \left( \nabla \times \mathbf{A}_{e} \right) \right\} dV} + \int_{\Omega_{e}}{i\omega\mu\sigma\mathbf{w}_{i} \cdot \mathbf{A}_{e}dV} = \mu\int_{\Omega_{e}}{\mathbf{w}_{i} \cdot \mathbf{i}_{s_{e}}dV} \cdot \cdot \cdot (5)$$

ここで $w_{i}\ (i = 1, \cdot \cdot \cdot ,6)$ は e 番目の四面体の各辺に対する辺要素基底関数です。e 番目の要素内のベクトルポテンシャルを次のように表します。

$$\mathbf{A}_{e}=\sum_{j = 1}^{6}{\mathbf{w}_{j}{\lbrack Al\rbrack}_{j}} \cdot \cdot \cdot \cdot (6)$$

ここで

$$\mathbf{w}_{j} = \mathbf{w}_{kl} = \lambda_{k}\nabla\lambda_{l} - \lambda_{l}\nabla\lambda_{k}$$

であり、$\lambda_{k}$ は節点基底関数です。j と k, l の対応関係を下図と表に示します。

<img src="./images/v1.4_image41.png" width="55%">

表A1.

| i, j（辺番号） | k   | l   | m   | n   |
|:--------------:|:---:|:---:|:---:|:---:|
| 1              | 1   | 2   | 3   | 4   |
| 2              | 2   | 3   | 1   | 4   |
| 3              | 1   | 3   | 4   | 2   |
| 4              | 1   | 4   | 2   | 3   |
| 5              | 2   | 4   | 3   | 1   |
| 6              | 3   | 4   | 1   | 2   |

式(6)を用いると、式(5)は次のようになります。

$$\sum_{j = 1}^{6}{\int_{\Omega_{e}}{\left\{ \left( \nabla \times \mathbf{w}_{i} \right) \cdot \left( \nabla \times \mathbf{w}_{j} \right) \right\} dV}\lbrack Al\rbrack_{j}} + \sum_{j = 1}^{6}{\int_{\Omega_{e}}{i\omega\mu\sigma\mathbf{w}_{i} \cdot \mathbf{w}_{j}dV}\lbrack Al\rbrack_{j}} = \mu\int_{\Omega_{e}}{\mathbf{w}_{i} \cdot \mathbf{i}_{s_{e}}dV},\  \cdot \cdot \cdot (7)$$

これを行列形式で書くと次のようになります。

$$\left\{ \left\lbrack M_{e} \right\rbrack + i\omega\mu\sigma_{e}\left\lbrack N_{e} \right\rbrack \right\}\left\lbrack \mathbf{Al}_{e} \right\rbrack = \mathbf{s}_{e} \cdot \cdot \cdot (8)$$

$$\left\lbrack \mathbf{Al}_{e} \right\rbrack = \left\lbrack {Al}_{1}\ {Al}_{2}\ {Al}_{3}\ {Al}_{4}\ {Al}_{5}\ {Al}_{6} \right\rbrack^{T} \cdot \cdot \cdot (9)$$

$$\left\lbrack M_{e} \right\rbrack_{ij} = \int_{\Omega_{e}}{\left\{ \left( \nabla \times \mathbf{w}_{i} \right) \cdot \left( \nabla \times \mathbf{w}_{j} \right) \right\} dV} = \frac{1}{3v}\mathbf{x}_{mn} \cdot \frac{1}{3v}\mathbf{x}_{m^{'}n^{'}} \times v = \frac{1}{9v}\mathbf{x}_{mn} \cdot \mathbf{x}_{m^{'}n^{'}} \cdot \cdot \cdot (10)$$

$$\left\lbrack N_{e} \right\rbrack_{ij} = \int_{\Omega_{e}}{\mathbf{w}_{i} \cdot \mathbf{w}_{j}dV} = \int_{\Omega_{e}}{\left( \lambda_{k}\nabla\lambda_{l} - \lambda_{l}\nabla\lambda_{k} \right) \cdot \left( \lambda_{k^{'}}\nabla\lambda_{l^{'}} - \lambda_{l^{'}}\nabla\lambda_{k^{'}} \right)dV}$$

$$= \int_{\Omega_{e}}{\left( \lambda_{k}\lambda_{k^{'}}\nabla\lambda_{l} \cdot \nabla\lambda_{l^{'}} - \lambda_{k}\lambda_{l^{'}}\nabla\lambda_{l} \cdot \nabla\lambda_{k^{'}} - \lambda_{l}\lambda_{k^{'}}\nabla\lambda_{k} \cdot \nabla\lambda_{l^{'}} + \lambda_{l}\lambda_{l^{'}}\nabla\lambda_{k} \cdot \nabla\lambda_{k^{'}} \right)dV}$$

$$= \nabla\lambda_{l} \cdot \nabla\lambda_{l^{'}}\int_{\Omega_{e}}{\lambda_{k}\lambda_{k^{'}}dV} - \nabla\lambda_{l} \cdot \nabla\lambda_{k^{'}}\int_{\Omega_{e}}{\lambda_{k}\lambda_{l^{'}}dV}$$

$$- \nabla\lambda_{k} \cdot \nabla\lambda_{l^{'}}\int_{\Omega_{e}}{\lambda_{l}\lambda_{k^{'}}dV} + \nabla\lambda_{k} \cdot \nabla\lambda_{k^{'}}\int_{\Omega_{e}}{\lambda_{l}\lambda_{l^{'}}dV} \cdot \cdot \cdot (11)$$

$$s_{i}^{e} = \mu\int_{\Omega_{e}}{\mathbf{w}_{i} \cdot \mathbf{i}_{s_{e}}dV} = \mu I_{s}\left( \frac{\mathbf{w}_{i}^{e}\left( r_{p1} \right) + \mathbf{w}_{i}^{e}\left( r_{p2} \right)}{2} \right) \cdot \overrightarrow{p_{1}p_{2}}\ \ (12)$$

$$= \frac{\mu I_{s}\left( \lambda_{k}\left( p_{1} \right)\nabla\lambda_{l} - \lambda_{l}\left( p_{1} \right)\nabla\lambda_{k} + \lambda_{k}\left( p_{2} \right)\nabla\lambda_{l} - \lambda_{l}\left( p_{2} \right)\nabla\lambda_{k} \right)}{2} \cdot \overrightarrow{p_{1}p_{2}}$$

$$= \frac{\mu I_{s}\left( (\lambda_{k}\left( p_{1} \right) + \lambda_{k}\left( p_{2} \right))\nabla\lambda_{l} - \left( \lambda_{l}\left( p_{1} \right) + \lambda_{l}\left( p_{2} \right) \right)\nabla\lambda_{k} \right)}{2} \cdot \overrightarrow{p_{1}p_{2}}$$

式(12)の導出では、下図のように、ヘヴィサイド関数で定義される直線電流が四面体を貫くとしています。

<img src="./images/v1.4_image42.png" width="70%">

式(9)の導出では、次の関係を用いました。

$$\nabla \times \mathbf{w}_{j} = 2\nabla\lambda_{k} \times \nabla\lambda_{l} = \frac{1}{3v}\mathbf{x}_{mn}\  \cdot \cdot \cdot (13)$$

(k,l) と (m,n) の対応は表A1に示しています。

式(11)の計算では、次の関係を利用できます。

$$\nabla\lambda_{k} = \frac{1}{6v}\mathbf{x}_{ln} \times \mathbf{x}_{mn} \cdot \cdot \cdot (14)$$

$$\int_{\Omega_{e}}{\lambda_{k}\lambda_{l}dV} = \left\{ \begin{matrix}
6v\frac{1!1!}{(1 + 1 + 3)!} = \frac{v}{20}\ (k \neq l) \\
6v\frac{2!}{(2 + 3)!} = \frac{v}{10}\ \ \ \ \ \ \ \ (k = l)
\end{matrix} \right.\  \cdot \cdot \cdot (15)$$

なお、四面体における節点基底関数の積分には次の公式があります。

$$\int_{\Omega_{e}}{\lambda_{1}^{k}\lambda_{2}^{l}\lambda_{3}^{m}\lambda_{4}^{n}dV} = 6v\frac{k!l!m!n!}{(k + l + m + n + 3)!}\ \  \cdot \cdot \cdot (16)$$

## 付録B：ActFEM 計算のための支配方程式のスケーリング <a id="AppendixB"></a>
式(7)と式(12)を組み合わせた離散化形式

$$\sum_{j = 1}^{6}{\int_{\Omega_{e}}{\left\{ \left( \nabla \times \mathbf{w}_{i} \right) \cdot \left( \nabla \times \mathbf{w}_{j} \right) \right\} dV}\lbrack Al\rbrack_{j}} + \sum_{j = 1}^{6}{\int_{\Omega_{e}}{i\omega\mu\sigma\mathbf{w}_{i} \cdot \mathbf{w}_{j}dV}\lbrack Al\rbrack_{j}} = \mu I_{s}\left( \frac{\mathbf{w}_{i}^{e}\left( r_{p1} \right) + \mathbf{w}_{i}^{e}\left( r_{p2} \right)}{2} \right) \cdot \overrightarrow{p_{1}p_{2}}\  \cdot \cdot \cdot (17)$$

を解く際、ActFEM ではスケーリングを行っています。スケーリング長 $L$ を用い、ActFEM では $L$ = 1000 m としています。次の置き換えを用いて式(17)を書き直すことを考えます。

$$\nabla = \frac{1}{L}\nabla^{'},\quad \mathbf{w} = \frac{1}{L}\mathbf{w}^{'},\quad dV = L^{3}dV^{'},\quad \overrightarrow{p_{1}p_{2}} = L\overrightarrow{p_{1}p_{2}'},\quad Al\left\lbrack s \cdot V/m \cdot m \right\rbrack = \int \mathbf{A} \cdot d\mathbf{l} = \frac{1}{L}Al^{'},$$

ここで $Al^{'}$ の単位は $\left\lbrack s \cdot mV/km \cdot km \right\rbrack$、また $\mathbf{E} = i\omega\mathbf{A}$ です。これらの置き換えは、**$\mathbf{w}$ と $\nabla$ が、km 単位の座標系で $\mathbf{w}'$ および $\nabla'$ として構成されている**ことに対応します。

$$\sum_{j = 1}^{6}{\int_{\Omega_{e}}{\left\{ \left( \frac{1}{L}\nabla' \times \frac{1}{L}\mathbf{w}_{i}^{'} \right) \cdot \left( \frac{1}{L}\nabla^{'} \times \frac{1}{L}\mathbf{w}_{j}^{'} \right) \right\} L^{3}dV^{'}}\frac{1}{L}\ \left\lbrack Al^{'} \right\rbrack_{j}} + \sum_{j = 1}^{6}{\int_{\Omega_{e}}{i\omega\mu\sigma\frac{1}{L}\mathbf{w}_{i}^{'} \cdot \frac{1}{L}{\mathbf{w}^{'}}_{j}L^{3}dV'}\frac{1}{L}\left\lbrack Al^{'} \right\rbrack_{j}} = \mu I_{s}\frac{1}{L}\left( \frac{w'\left( r_{p1} \right) + w'\left( r_{p2} \right)}{2} \right) \cdot L\overrightarrow{p_{1}p_{2}'}$$

$$\Longleftrightarrow\sum_{j = 1}^{6}{\int_{\Omega_{e}}{\left\{ \left( \nabla^{'} \times \mathbf{w}_{i}^{'} \right) \cdot \left( \nabla^{'} \times \mathbf{w}_{j}^{'} \right) \right\} dV^{'}}\left\lbrack Al^{'} \right\rbrack_{j}} + \sum_{j = 1}^{6}{\int_{\Omega_{e}}{L^{2}i\omega\mu\sigma\mathbf{w}_{i}^{'} \cdot {\mathbf{w}^{'}}_{j}dV^{'}}\left\lbrack Al^{'} \right\rbrack_{j}} = \mu I_{s}L^{2}\left( \frac{w'\left( r_{p1} \right) + w'\left( r_{p2} \right)}{2} \right) \cdot \overrightarrow{p_{1}p_{2}'}$$

なお、送信源の電流 $I_{s}$ はスケーリングしておらず、単位は [A] のままです。
