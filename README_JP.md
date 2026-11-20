# SPADExp4D

[English](README_EN.md) | 日本語

現在のバージョン：**1.3.0**

SPADExp4Dは、第一原理計算コードOpenMXの出力を利用して、角度分解光電子分光（ARPES）の光電子角度分布（PAD）を4次元空間 `(kx, ky, kz, E)` で計算・可視化する研究用ソフトウェアです。

本ソフトウェアはHiroaki Tanaka氏らによる[SPADExp](https://github.com/Hiroaki-Tanaka-0606/SPADExp)を基に、第一Brillouin zone内の固有状態キャッシュ、3次元運動量グリッド、エネルギー依存PAD、計算結果管理などを追加した派生・改良版です。

## 主な機能

- 第一Brillouin zone内の固有値・LCAO係数をHDF5へ保存し、再利用
- `(kx, ky, kz, E)` 上のPAD計算と断面表示
- 通常のOpenMX軸順序（`a1/a2` が面内、`a3` が表面法線）に対応
- CIFから計算対象を追加し、構造ごとに設定・キャッシュ・結果を管理
- C++一括PAD計算と、比較用Pythonバックエンド
- Gaussian平滑化、Fermi準位シフト、偏光・Euler角・光子エネルギーの設定
- 計算済み固有状態を使った再計算の省略と中断後の再開
- ローカル計算と大型計算機との手動ファイル受け渡し

## 公開範囲

このリポジトリは**ソースコードのみ**を公開します。次のファイルは含まれません。

- コンパイル済み実行ファイル（`bin/`）
- ビルド生成物（`.build/`）
- SPADExp/OpenMXのHDF5データ
- OpenMXのPAO・VPSデータセット
- 計算キャッシュ、ログ、途中結果、ユーザー設定

clone直後の状態だけでは計算を実行できません。以下の依存ソフトウェアとデータを利用者が用意してください。

## 動作環境と依存関係

SPADExp4Dは **LinuxおよびmacOS** で動作します。実行には、OpenMX、C++17対応のコンパイル・実行環境、およびPython 3.12が動作する環境が必要です。

本リリースの動作確認はmacOSで行いました。同梱のMakefileはOSを自動判定し、macOSではlibompとApple Accelerate framework、LinuxではコンパイラのOpenMPとBLAS/LAPACKを使用します。

### C/C++

- C++17対応コンパイラ（Clangを想定）
- GNU Make
- HDF5（C++ APIを有効にしたもの）
- OpenMP runtime（macOSでは通常libomp）
- BLAS/LAPACK（macOSではApple Accelerate framework、Linuxでは利用環境に対応する実装）

### Python

- Python 3.12（動作確認：3.12.7）
- PyQt6（動作確認：6.11.0、Qt 6.11.0）
- Matplotlib 3.11.1（PyQt6と`QtAgg`バックエンドで動作確認済み）
- NumPy
- h5py
- pytest（テストを実行する場合）

### 外部ソフトウェアとデータ

- OpenMX実行ファイル `openmx`
- OpenMXソース配布物に付属する`cif2omx`（コンパイルして使用）
- OpenMXの `DFT_DATA19` または互換PAO・VPSデータ
- `data/Atomic_potentials.hdf5`
- `data/PAO_and_AO_after_opt.hdf5`

HDF5ファイルは元のSPADExp配布物など、利用条件を確認できる正規の入手元から取得してください。

## フォルダ構成

```text
SPADExp4D/
├── cpp/             C++計算エンジンとOpenMX変換ツール
├── python/          GUI、CLI、4D計算制御
├── structures/      CIFと構造ごとの初期設定
├── tests/           自動テスト
├── exchange/        大型計算機との手動受け渡し領域
├── data/            実行時HDF5・OpenMXデータ（Git管理外）
├── bin/             ローカルで生成される実行ファイル（Git管理外）
├── .build/          ビルド中間生成物（Git管理外）
├── Makefile
└── LICENSE
```

## セットアップ（外部ソフトウェアと本体を分離）

### A. 利用者が別途インストールするもの

以下は本リポジトリに含まれません。

- OpenMX本体（ソース配布物に含まれる`cif2omx`もコンパイルします）
- C++17対応コンパイラ、GNU Make、HDF5（C++ API）、OpenMP、BLAS/LAPACK
- Python 3.12、PyQt6、Matplotlib 3.11.1、NumPy、h5py（テストにはpytest）

#### OpenMXと付属cif2omx

OpenMXをインストールした後、OpenMXソースの`source`ディレクトリで付属の`cif2omx`をコンパイルします。

```bash
cd /path/to/openmx/source
gcc -std=c99 -O2 -Wall -Wextra -pedantic \
  -DCIF2OPENMX_MAIN cif2omx.c -lm -o cif2omx
```

生成した`openmx`と`cif2omx`を`PATH`に置くと自動検出します。明示的に指定する場合は次の環境変数を使います。

```bash
export SPADEXP4D_OPENMX=/path/to/openmx
export SPADEXP4D_CIF2OMX=/path/to/cif2omx
```

環境変数と `PATH` のどちらからも見つからない場合は、必要になった時点で場所の指定を求めるエラーになります。

#### Pythonパッケージ

```bash
python3 -m venv .venv
source .venv/bin/activate
python3 -m pip install --upgrade pip
python3 -m pip install PyQt6 'matplotlib==3.11.1' numpy h5py pytest
```

#### HDF5、OpenMP、BLAS/LAPACK

OSまたはパッケージ管理ツールでインストールしてください。HDF5を`pkg-config`で検出できない場合は、ビルド時に指定します。

```bash
make HDF5_PREFIX=/path/to/hdf5 LIBOMP_PREFIX=/path/to/libomp
```

### B. 本パッケージへ配置する実行時データ

```text
data/
├── Atomic_potentials.hdf5
└── PAO_and_AO_after_opt.hdf5
```

OpenMXデータがOpenMX本体の隣の `DFT_DATA19/` にない場合は、`data/openmx/` に配置できます。

### C. SPADExp4D本体のビルドと実行

リポジトリのルートで実行します。

```bash
make
make check
```

実行ファイルは `bin/`、中間生成物は `.build/` に作成され、Gitには登録されません。

## 実行方法

### GUI

```bash
python3 python/_main.py
```

初期設定ファイルがない場合は組み込みの既定値を使用します。GUIで変更した設定は構造ごとの `parameter.json` などへ保存されます。


## 結晶構造の追加

GUIの「Import CIF…」から追加します。
下記のフォルダに計算の状況が保存されます。

```text
structures/<structure_name>/
├── <structure>.cif
├── parameter.json
├── cache/
└── results/
```

## GUIの操作について

### 1. 計算対象の選択と構造の追加

起動後、`計算対象`から構造を選択します。新しい構造を追加する場合は、メニューの`Import CIF…`からCIFファイルを選択し、構造名を入力します。構造ごとの入力ファイルと計算結果は、`structures/<structure_name>/`以下に保存されます。

### 2. 第一Brillouin zoneの準備

`N1, N2, N3`で第一Brillouin zoneの分割数を設定します。分割数を大きくすると精度が向上する一方、OpenMX計算時間とHDF5ファイル容量が増加します。`Calc`を押すと、必要な固有状態がOpenMXで計算されます。既に同じ条件のキャッシュがある場合は再利用され、不足分だけが計算されます。

`OpenMX Parameters`タブでは、計算対象ごとにOpenMXのキーワードを変更できます。

- `scf.XcType`、`scf.SpinPolarization`（off / on / NC）、`scf.SpinOrbit.Coupling`（NCのときのみ）、`scf.energycutoff`、`scf.Kgrid`、`scf.ElectronicTemperature`、`scf.criterion`、`scf.maxIter`
- 初期値はCIFから生成した`.openmx.base.dat`の値です。変更した値だけが太字になり、`parameter.json`の`"openmx"`に保存されます。変更すると、次の`Calc`で新しいSCFとk-space HDF5が計算されます（SCF restartもキーワードの組ごとに別になります）。
- `Reset to recommended`は単位胞から推奨値を設定します（`scf.Kgrid`はk点間隔0.2 Å⁻¹以下、`scf.maxIter`は200以上など。テンプレートより弱い値にはしません）。
- `scf.SpinPolarization on`では上下スピンのPADを合計し、`NC`ではスピノルの両成分を合計します。常磁性でスピン軌道相互作用を入れる場合は`NC`と`scf.SpinOrbit.Coupling`を使います。
- `Stored bands E−EF`でk-space HDF5に保存するバンドのエネルギー範囲を指定できます（`All`は全バンド）。NCではバンド数と係数が増えるため、PADを計算するEB範囲を含む窓を指定するとHDF5とメモリを数分の一にできます。窓の外の状態はPADに含まれません。

OpenMXで計算する前に「Calculate the eigenstates of <計算対象> with OpenMX now?」という確認が表示されます。k点の分割数と総数、中断した計算から保存済みのk点数（あれば）と、これから計算するk点数が示されます。計算を開始するには`Yes`、分割数を見直すなど今は計算しない場合は`No`を選択してください。OpenMXの実行中はログと進捗が画面に表示され、必要に応じて`Cancel`で停止できます。

`Results`タブの`K-space HDF5`には、PADに使ったk-space HDF5ファイル、計算日時、OpenMXとHDF5変換の時間（合計と実時間）、ファイルサイズが表示されます。フルパスはツールチップに表示されます。計算時間はk-space HDF5ファイル自体（属性`CalculationTiming`）に保存されるため、キャッシュを再利用したとき、再起動後、保存したPAD結果を開いたときにも表示されます。中断して再開した計算では、中断前の時間も合計されます。このバージョンより前に作ったキャッシュには計算時間が記録されていないため、計算日時とサイズだけが表示されます。

### 3. PAD条件（`PAD`タブ）

OpenMXのタブより下の設定は、`PAD`タブと`Band Diagram`タブに分かれています。`PAD`タブには次の順に並んでいます。

- `PAD calculation`：`C++ batch`または比較用の`Python per EB (comparison)`
- `hν`：光子エネルギー
- `Polarization`：`lx`、`ly`、`c+`、`c-`、`unpolarized`
- `Incident Angles (zyz Euler)`：入射光のEuler角a, b, c（z-y-z回転）。a・cは−360〜360°、bは0〜180°です。`Incident Angles (zyz Euler, °)…`ボタンを押すとビームラインのプリセット（SPring-8 BL25SU：0, 90, −90、UVSOR PMM：0, 68, −90）を選べ、OKで値が入ります
- `EF shift`：表示時のFermi準位シフト
- `V0`：物質内部ポテンシャル。`V0 <= 0`では自動設定を使用します
- `wfnc`：光電子波動関数のパラメーター
- `Slice direction`：PADとバンドのスライス方向（`kx`または`ky`）
- `ky fixed`（`ky`スライスでは`kx fixed`）：スライスの位置を決める固定座標（4.を参照）
- `KMin, KMax, Kstep`：表示する運動量範囲と刻み
- `EB Begin, End, Step`：結合エネルギーの範囲と刻み
- `Gaussian dK`、`Gaussian dE`：保存済みPADに運動量方向・エネルギー方向のGaussian平滑化を適用します。変更してもOpenMXやPADの再計算は不要です
- `Maximum angle θmax`：この放出角より外側のPADは表示しません
- `Boundary`：第一Brillouin zoneの境界を重ねて表示します

PAD計算条件を変更しても、第一BZの固有状態キャッシュは再利用されます。

### 4. スライスの座標と図（`Band Diagram`タブ）

`Band Diagram`タブはバンド図の設定用で、まだ実装していないため、このバージョンでは選択できません（無効）。このタブにあるのは2つ目の固定座標だけで、`parameter.json`（または開いた結果ファイル）に保存された値のままになります。1つ目の固定座標は`PAD`タブの`Slice direction`の下で入力します。`Slice direction`で選んだスライスの固定座標は次のとおりです。`kx`スライスでは`ky fixed`（`PAD`タブ）でkyを指定します。kzは固定ではなく、各EBで光電子の球面 kz = √((hν + V0 − EB)/(ħ²/2m) − kx² − ky²) に従うため、`kz fixed`は使いません（`ky`スライスも同様で、`kx fixed`を指定します）。`Slice direction`で選べるのは`kx`と`ky`だけです。以前のバージョンで保存した`kz`スライスの設定は`kx`として開き、`kz`スライスの結果ファイルはそのまま表示できます。

図の並びは、`Calc`の上の`Display type`で選びます（計算し直す必要はありません）。

| Display type | 図 |
|---|---|
| `Normal`（既定） | 2×2：左上`Band Diagram EB′`（選択したEBでのバンドの断面）、右上`PAD: EB′`、左下`Band Diagram (PAD cut)`、右下`PAD Slice kx-EB` |
| `With Photoelectron kz` | 3×2：`Band Slice EB′`と`PAD: EB′`（選択したEBでのkx–kyマップ）、`Band Slice kx-EB`と`PAD Slice kx-EB`（スライスのEB依存性）、`Band Diagram (PAD cut)`、kx–kz光電子軌跡 |

`Band Diagram EB′`は、選択したEBでPADのマップが通る光電子の球面上（kz = √((hν + V0 − EB)/(ħ²/2m) − kx² − ky²)）で、EB = −E_n(kx, ky, kz) となるバンドの断面を線で描きます。線の色は、その位置のPAD（`PAD: EB′`と同じGaussian dK/dE適用後の値）を表し、範囲はそのEBで表示されているPADの最小値から最大値までです（レインボー、カラーバー`PAD`）。PADがない部分（光電子の脱出限界やθmaxの外側）の線は灰色です。断面線は計算のたびに全EBについて求めて結果に保存します。断面線を含まない以前の結果は、k-space HDF5があれば表示するときに求めます。

どちらの表示でも、`Band Diagram (PAD cut)`はkzを固定したバンド図ではありません。PADのスライスが実際に通るk点、つまり EB = −E_n(kx, ky, kz(EB)) となる点のバンドを描くので、`PAD Slice kx-EB`とそのまま比較できます。各バンドは、スライス上のその位置のPAD（Gaussian dK/dE適用後）でレインボーカラーに色分けされます。色の範囲は、表示範囲内のスライスのPADの最小値から最大値までです（カラーバー`PAD`）。スライスにPADがない部分（光電子の脱出限界やθmaxの外側）のバンドは灰色です。すべての図は上段のkx–kyマップと同じ幅・同じ横位置で並ぶので、各列のk軸がそろいます。

### 5. 結果の再利用と保存

計算済みの構造・メッシュ・PADは構造ごとの`cache/`と`results/`に保存されます。同じ条件で再度`Calc`を押すと、利用可能なキャッシュが自動的に再利用されます。メニューの`計算結果を再読み込み`で結果一覧を更新できます。

設定は構造ごとの`parameter.json`などへ保存されます。計算結果を別の計算機で処理する場合は、`受け渡し`メニューからアップロード用ジョブを作成し、返却された結果をダウンロード領域から取り込みます。

#### 前回の状態の復帰

終了時に、計算対象、左側のすべてのパラメーター、選択中のタブ、`Display type`、EBの位置、ウィンドウの大きさを`python/spadexp_4d.ini.json`に保存します。また、計算が完了するたびに、その結果を`<計算対象>/cache/last_pad_result.h5`へ保存します。次回の起動時には同じ計算対象が選ばれ、前回のPAD結果が再計算なしで表示されます。計算後に変更したGaussian dK/dE、EF shift、θmax、Boundaryもそのまま戻ります。左側のパラメーターは、終了時に表示していた値がそのまま戻ります。組み合わせが条件を満たさない場合（例：EB End − EB BeginがEB Stepで割り切れない）も入力した値のまま戻り、`Calc`の上に赤字で理由が表示されます（直すまで`Calc`はできません。`parameter.json`には最後の有効な値が残ります）。計算対象を切り替えたときは、その対象に保存したパラメーターが使われます。`Calculation result`で計算結果を選んだときだけ、その計算の条件が読み込まれます。終了時にPAD結果ファイルを開いていた場合は、そのファイルを表示します。

#### PAD結果の保存と読み込み（HDF5）

`File > Save PAD result as HDF5…`で、表示中の結果を1つのHDF5ファイルに保存します。`File > Open PAD result (HDF5)…`で読み込むと、その計算対象が選ばれ、パラメーターも保存時の値に戻ります。計算対象がデータフォルダーにない場合は、ファイルに保存された結晶構造から追加できます。ファイルの内容は次のとおりです。

| グループ | 内容 |
|---|---|
| `/`（属性） | `format = "SPADExp4D PAD result"`、`format_version`、`software`、`created`、`structure_name` |
| `Settings` | GUIとOpenMXのすべての設定（属性`json`）と、各スカラー設定の属性 |
| `PAD` | `binding_energy`（eV）、`k`（Å⁻¹、kxとky共通）、平滑化前の`pad`と`band_weight`（EB, ky, kx、gzip圧縮）、`nearest_eigenstate_distance`。属性`v0_eV`、`request_json` |
| `BandCut` | PADスライスに沿ったバンド：線は`k`と`energy`（E − EF、eV。スライス外はNaN）、同じkで複数回交わる点と`kz`スライスは`points`（k、E − EF） |
| `Structure` | `cif`、`openmx_template`（`.openmx.base.dat`）、`lattice`（Å）、`reciprocal`（Å⁻¹） |
| `BandContours` | 各EBのバンドの断面線（`Band Diagram EB′`）：`binding_energy`、`offsets`（EB iの線分は`offsets[i]:offsets[i+1]`）、`segments`（線分の両端のkx, ky、Å⁻¹） |
| `Calculation` | k-space HDF5のパスとサイズ、OpenMX・変換・PADの計算時間、k-space HDF5の計算記録（`mesh_timing`、JSON）、固有状態までの最大距離 |

Gaussian dK/dEは設定として保存し、表示するときに適用します。そのため、保存されるPADは平滑化前の値です。k-space HDF5（固有状態）は含まれません。条件を変えて再び`Calc`を押す場合は、k-space HDF5が必要です（ない場合は再計算します）。


## 既知の制約

- LinuxとmacOSに対応していますが、本リリースの動作確認はmacOSで行いました。
- Linuxのディストリビューションや導入方法によっては、HDF5、OpenMP、BLAS/LAPACKのライブラリ名や検索先の指定が必要です。
- バイナリと実行時HDF5は配布されません。
- OpenMXおよびそのデータファイルは別途必要です。
- 大きな3次元k点グリッドは計算時間とHDF5容量が急増します。
- 自動CIF変換は部分占有・無秩序サイト・未定義元素に対応しません。使えないCIFの場合、`Import CIF`は問題点を一つずつ英語で説明し、構造を追加しません。対象は、セルの値の不足や誤り、原子の一覧や分率座標がないこと、元素ではない記号（構造のひな形の`T`や`X1`のような仮の記号など）、部分占有のサイト、同じ位置に重なった原子、OpenMXに基底関数がない元素です。cif2omxのメッセージは「Show Details」で見られます。
- `scf.SpinPolarization on`の初期スピン配置はCIFから生成した値（上下同数）のままです。磁気秩序を扱うには`.openmx.base.dat`の初期スピンや磁気単位胞を用意してください。

## 引用

SPADExp4Dを研究で使用した場合は、元のSPADExp論文を引用してください。

H. Tanaka, K. Kuroda, and T. Matsushita, “SPADExp: A photoemission angular distribution simulator directly linked to first-principles calculations,” *Journal of Electron Spectroscopy and Related Phenomena* **264**, 147297 (2023). https://doi.org/10.1016/j.elspec.2023.147297

## 著作権、帰属、ライセンス

SPADExp4Dの2026年の変更部分：

```text
Copyright (C) 2026 Matsushita Tomohiro
Nara Institute of Science and Technology (NAIST)
```

本ソフトウェアはHiroaki Tanaka氏らのSPADExpを基にした改変版です。元の著作権表示、ライセンス表示、各ソース内の第三者表示は保持してください。`sphere_lebedev_rule.cpp` にはGNU LGPLで配布されたJohn Burkardt氏のコードと表示が含まれます。

SPADExp4Dは元のSPADExpと同じ **GNU General Public License version 3** で配布します。無保証です。詳しい条件は[LICENSE](LICENSE)を参照してください。

## ChatGPTの利用

2026年のSPADExp4D開発では、コードレビュー、リファクタリング、テスト、公開準備および文書作成の支援ツールとしてOpenAIのChatGPTを使用しました。ChatGPTは著作権者または保守責任者ではありません。公開内容と研究上の判断に対する責任は著作権者にあります。

## サポートと貢献

問題を報告する際は、OS、Python・コンパイラ・OpenMXの版、使用した構造、再現手順、エラーログを添えて、このリポジトリのIssue trackerを使用してください。大容量HDF5、計算結果、個人情報、認証情報は添付しないでください。

変更を提案する場合は、生成物を含めず、関連テストと説明を添えてください。投稿された変更は本プロジェクトのGNU GPL version 3の条件で配布できるものとします。
