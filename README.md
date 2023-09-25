# lascut.py

.lasポイントデータから、指定した範囲を切り出すスクリプトです。

東京、静岡、長崎等で公開されている、レーザー測量による点群データのlasファイルを使うことを想定しています。

 - 範囲指定は中心座標(緯度経度)と範囲(m)で指定
 - lasファイルの名前を自動計算
 - 9つのファイルにまたがる範囲まで対応
 - 出力はオレオレcsv仕様なので、カスタマイズしてね