# 2025/302 位相残差の独立検証

yi-corr で未適用の対流圏遅延が、提示された三次関数的な位相差の主因であることを強く支持する結果を得た。他の相関器の実装は参照せず、比較にはユーザー提示の実行結果だけを使用。本番コードは変更していない。

## 条件と方法

元の XML はこの環境で見つからず、会話中の観測条件を再構成した。正確な入力値は report.json に保存。搬送波6.6 GHz、サンプリング1.024 GHz、開始MJD60977.34375、3600秒、DUT1・極運動ゼロ、TT−UTC 69.184秒。

probe.rs は yi-corr の geom.rs と model_diag.rs を使用。Python側では C ERFA の CIO 座標変換と光行差を独立計算。IERS式11.9は運動項のみで、重力遅延・重力ポテンシャルを含む完全モデルではない。

大気は ESA の公開された公称天頂経路長とマッピング関数を組み合わせた感度試験：

- Z = 2.3 exp(-0.116e-3 H) + 0.1 [m]
- M(E) = 1.001 / sqrt(0.002001 + sin(E)^2)
- delta_tau = (Z2 M(E2) - Z1 M(E1)) / c
- 位相 = 360 × 6.6e9 × delta_tau [deg]

実測気象を使わず、楕円体高を標高の近似に用いた。仰角は probe 計算値（report.json参照）を使用し、提示表の仰角と完全一致する再現ではない。物理係数は観測に合わせて調整していない。troposphere_free_scale は別の診断値で、以下の予測には使用しない。

各出力が独立に rate/acel 補正されているため、定数・一次・二次の差を除去する。13時刻で二次多項式だけを推定し、別の12時刻では多項式も物理係数も固定。転記値は paired_output_rows.txt と paired_delay_rows.txt、選択時刻と処理は各スクリプトに保存。

## 結果

| 指標 | 大気補正前 | 公称大気補正後 |
|---|---:|---:|
| 13時刻の位相差RMS（二次成分除去） | 57.600° | 1.596° |
| 同peak-to-peak | 190.857° | 5.255° |
| 別の12時刻の固定予測誤差RMS | — | 1.332° |
| 残留遅延差RMS（二次成分除去、13時刻） | 0.025014 sample | 0.001009 sample |

大気基線差は約2.926 nsから6.443 nsへ変化。仰角低下に伴う非線形な大気遅延から二次成分を除くと、三次以上が残る。

360時刻で二次成分を除いた座標変換経路差はpeak-to-peak 0.000091°、光行差の一次近似との差は0.00650°。IERS運動項の差は9.30°で、主要な観測差を説明しない。補間1080試験点での最大誤差は0.191°。

提示出力の25時刻と再構成条件による検証であり、XML読み込みからの再現や元の相関データ全体の再処理ではない。残る約1–2°の原因までは確定していない。

## 再現

リポジトリルートで `python3 docs/validation/vlbi-2025302/run.py` を実行。Rust、既存 target/debug/deps/liberfa-*.rlib、Python numpy・pyerfa・matplotlib が必要。report.json と phase_comparison.png を更新し数値範囲を検査。本環境で成功を確認。

## 公開資料

- [ESA Tropospheric Delay](https://gssc.esa.int/navipedia/index.php/Tropospheric_Delay)：式4・9。
- [ERFA c2t06a](https://raw.githubusercontent.com/liberfa/erfa/master/src/c2t06a.c)
- [ERFA ab](https://raw.githubusercontent.com/liberfa/erfa/master/src/ab.c)
- [IERS Conventions 2010](https://iers-conventions.obspm.fr/conventions/content/tn36.pdf)：式11.9・11.10。
