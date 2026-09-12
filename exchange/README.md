# 4D Plus 計算結果の受け渡し

- `upload/`: この中のジョブフォルダをZIP化して大型計算機へアップロード
- `download/`: 大型計算機から戻ったZIPを展開して置く
- `../structures/<structure_name>/results/<calculation_id>/`: 検証後に展開した正式な結果

ZIPの作成・転送は手動で行います。`upload/<structure_name>_<calculation_id>/` を
丸ごとZIP化し、計算後は同じ名前のフォルダを `download/` に置きます。
OpenMX結果だけを先に戻してこのMacでPADを実行することもできます。受け取ったZIPはmanifestとSHA-256を確認して
から `structures/<structure_name>/results/<calculation_id>/` に配置し、その構造の
`results/active.json` で表示対象を選びます。
