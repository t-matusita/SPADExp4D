# 計算対象の結晶構造

計算対象ごとにフォルダを作り、その直下へCIFファイルを置きます。ファイル名は
`structure.cif`でなくてもよく、`TaS2.cif`など任意の`.cif`名を使用できます。
例: `structures/graphite_ab/graphite.cif`。
`parameter.json`がない場合は、`_main.py`が初回選択時に初期設定から自動作成します。
PAD・表示パラメーターはStructureを切り替える直前とアプリ終了時に、この
`parameter.json`へ自動保存されます。切替後は選択先の値で画面を更新します。
`_main.py` 起動時にこの単位フォルダを読み込み、パラメーター欄最上段の「計算対象」
ComboBoxへ `graphite_ab` を登録します。

GUIの「ファイル」→「Import CIF…」を使うと、CIF選択後にStructureのフォルダ名を
入力できます。CIF、`parameter.json`、OpenMX入力、`results/openmx/`を自動作成し、
初期パラメーターには最後に保存されたStructureの値を使用します。

フォルダ名がユーザー定義の `structure_name` です。CIFの内容はSHA-256で検証し、
変更された場合は別の計算として扱います。部分占有・無秩序サイト・未定義元素は
自動計算の対象外です。固有状態HDF5は各構造フォルダ内の
`results/openmx/mesh_<hash>/eigenstates.hdf5`、結果登録情報は
`results/<calculation_id>/`、SCF再開データは`results/<structure_name>_scf_rst/`
に保存します。計算完了後、OpenMXの一時作業ファイルは削除します。
大型計算機との受け渡し時だけ、CIFはジョブ内で
`structure.cif`という標準名にコピーされます。
