# SPADExp4D

[English](README_EN.md) | 日本語

現在のバージョン：**1.1.0**

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

計算前に確認ダイアログが表示されます。計算を開始するには`Yes`を選択してください。OpenMXの実行中はログと進捗が画面に表示され、必要に応じて`Cancel`で停止できます。

### 3. PAD表示条件

- `hν`：光子エネルギー
- `V0`：物質内部ポテンシャル。`V0 <= 0`では自動設定を使用します
- `EF shift`：表示時のFermi準位シフト
- `wfnc`：光電子波動関数のパラメーター
- `EB Begin, End, Step`：結合エネルギーの範囲と刻み
- `偏光`：`lx`、`ly`、`c+`、`c-`、`unpolarized`
- `Euler a, b, c`：偏光・結晶座標の回転角

PAD計算条件を変更しても、第一BZの固有状態キャッシュは再利用されます。`PAD計算`で`C++ 一括計算`または比較用の`Python EB毎`を選択できます。

### 4. 断面表示と平滑化

`断面方向`で`kx`、`ky`、`kz`のいずれかを選択し、表示する断面の`a`、`b`座標を指定します。`KMin, KMax, Kstep`は表示する運動量範囲と刻みです。`Brillouin zone表示`を有効にすると、第一Brillouin zoneの境界が重ねて表示されます。

`Gaussian dK`と`Gaussian dE`で、保存済みPADに運動量方向・エネルギー方向のGaussian平滑化を適用できます。平滑化条件の変更は、OpenMXやPADの再計算を必要としません。

### 5. 結果の再利用と保存

計算済みの構造・メッシュ・PADは構造ごとの`cache/`と`results/`に保存されます。同じ条件で再度`Calc`を押すと、利用可能なキャッシュが自動的に再利用されます。メニューの`計算結果を再読み込み`で結果一覧を更新できます。

設定は構造ごとの`parameter.json`などへ保存されます。計算結果を別の計算機で処理する場合は、`受け渡し`メニューからアップロード用ジョブを作成し、返却された結果をダウンロード領域から取り込みます。


## 既知の制約

- LinuxとmacOSに対応していますが、本リリースの動作確認はmacOSで行いました。
- Linuxのディストリビューションや導入方法によっては、HDF5、OpenMP、BLAS/LAPACKのライブラリ名や検索先の指定が必要です。
- バイナリと実行時HDF5は配布されません。
- OpenMXおよびそのデータファイルは別途必要です。
- 大きな3次元k点グリッドは計算時間とHDF5容量が急増します。
- 自動CIF変換は部分占有・無秩序サイト・未定義元素に対応しません。
- スピン偏極計算の全バンド4Dキャッシュには現在対応していません。

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
