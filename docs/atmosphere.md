# 電離層・中性大気の補正（3.11.0）

`yi-corr` と `yi-phasedarray` に局別の大気入力を渡せる。電離層は一次の冷たいプラズマ、中性大気（対流圏・成層圏の積分経路）は非分散の遅延として扱う。

```bash
target/release/yi-corr --sc observation.xml --raw /observations/raw --cor cor \
  --vlbi --atmosphere-file atmosphere.csv
```

`--vlbi` は従来どおり PNM06A+GAST、VLBI minus、標準大気の対流圏を選ぶ。入力された天頂遅延・気象値があれば標準大気を置き換え、TEC があれば電離層を有効にする。入力・明示指定がない場合は従来と同じ計算になる。

## 入力

XML の局 `key` と CSV の `station` を一致させる。局名・天体名で推測しない。基線の両局の値が必要であり、ゼロも `0` と明示する。STEC は**その天体方向の視線 TEC**、VTEC は薄層に水平一様と仮定した鉛直 TEC である。両者を混在して指定しない（局 A は STEC、局 B は VTEC という指定は可能）。複数の天体方向に異なる STEC を与える場合は、方向ごとにファイル・実行を分ける。

局別定数の例（数値は説明用で、実測値ではない）:

```csv
station,mjd_utc,stec_tecu,zhd_m,zwd_m
A,,10.0,2.30,0.10
B,,10.2,2.28,0.12
```

`mjd_utc` を空欄にした1行は、局ごとの定数を意味する。日時を入れた場合は各局2行以上が必要。UTC MJD の順に線形補間し、観測範囲を覆わない入力では停止する。外挿はしない。任意順の行を受け付けるが、同一局の重複時刻・途中での列の切り替えはエラーになる。

```csv
station,mjd_utc,vtec_tecu,pressure_hpa,temperature_k,vapour_pressure_hpa
A,60337.000000000000,10.0,1010.0,288.15,10.0
B,60337.000000000000,10.2,1000.0,289.15,11.0
A,60337.003472222222,10.1,1009.8,288.20,10.1
B,60337.003472222222,10.3,999.8,289.20,11.1
```

使える列は上記の `station,mjd_utc,stec_tecu,vtec_tecu,zhd_m,zwd_m,pressure_hpa,temperature_k,vapour_pressure_hpa`。未使用の列を省略できる。`#` 行はコメント、引用符付き CSV は対象外。単位は TECU、m、hPa、K で固定。NaN/Inf、負の TEC、不完全な気象値などは拒否する。

- `stec_tecu` または `vtec_tecu`: 1 TECU = 10¹⁶ electrons/m²。
- `zhd_m` と `zwd_m`: 天頂の静水圧・湿潤の追加光路長。ゼロの湿潤成分は `zwd_m=0`。
- 気象値を使う場合は、気圧・絶対温度・水蒸気分圧の3つを渡す。Saastamoinen の天頂遅延を計算し、Niell で視線方向へ写像する。相対湿度を水蒸気分圧と取り違えない。天頂遅延との同時指定は不可。

CSV のほかに XML の局内定数も指定できる:

```xml
<station key="A">
  <name>ANT1</name>
  <!-- 既存の局座標・terminal 等 -->
  <atmosphere stec-tecu="10.0" zhd-m="2.30" zwd-m="0.10"/>
</station>
<station key="B">
  <name>ANT2</name>
  <atmosphere stec-tecu="10.2" zhd-m="2.28" zwd-m="0.12"/>
</station>
```

XML では `_` を `-` に替えた属性名を使う。CSV をスケジュールから指定する場合は `<schedule>` 直下に次を置く。

```xml
<atmosphere file="atmosphere.csv" ionosphere="input" troposphere="input"
            shell-height-km="450" min-elevation-deg="5"/>
```

XML のファイル参照は XML 所在ディレクトリ基準、CLI のファイル参照は実行ディレクトリ基準。ファイル入力は局内定数を置き換える。CLI は対応する XML 設定より優先する。CLI のみで局 key が得られない入力では `ant1` / `ant2` を使う。

| オプション | 値・既定の動作 |
|---|---|
| `--atmosphere-file CSV` | ファイルと局内定数を置き換える |
| `--ionosphere off\|input` | 未指定なら TEC 入力時に有効 |
| `--troposphere off\|nominal\|input` | 未指定なら気象・天頂遅延を優先し、なければ `--vlbi` の標準大気 |
| `--ionosphere-shell-height-km KM` | VTEC 用の薄層高度、既定450 km |
| `--atmosphere-min-elevation-deg DEG` | VTEC/Niell の仰角下限、既定5°、指定可能3–90° |

## 群遅延と位相

基線を `X1 conj(X2)`、局差を `station2 − station1`、光速を `c`、周波数を Hz とする。

```
D = 40.3 * (STEC2 − STEC1) * 10^16 / c        [s Hz²]
ionosphere group delay = +D / f²              [s]
ionosphere phase delay = −D / f²              [s]
required XCF phase correction = +2π D / f      [rad]
```

1 TECU の差は100 MHzで40.3 m（約134.4 ns）の群経路差になる。TEC を既存の非分散遅延へ足すだけでは位相の符号が逆になる。

この実装は共通受信帯域の中心 `f_ref` における群遅延 `G=D/f_ref²` を、読み出し・整数/小数サンプル遅延へ加える。この経路が作る非分散位相 `−2π f G` を打ち消す分も含めて、各 FFT ビンに

```
extra_phase = +2π D * (1/f + f/f_ref²)
```

を掛ける。両者の合計が必要な `+2πD/f` になり、群遅延と位相を二重に数えない。`f=XML stream/frequency + output_bin*fs/N` の物理 RF を使い、電離層補正によるビン移動・ピークへのフィッティングは行わない。位相合成では、既存の遅延補正を受ける片側の電圧に同じ基線補正を割り当てる。

対流圏は `τ=(ZHD₂ m_h₂ + ZWD₂ m_w₂ − ZHD₁ m_h₁ − ZWD₁ m_w₁)/c`。非分散なので群遅延と位相遅延は同符号。Niell の季節項・南北半球・緯度補間・乾燥成分の高度補正を使う。局高さは WGS84 楕円体高を標高の近似に使い、ジオイド補正・水平勾配は含めない。

入力値は UTC で補間する。幾何投影の1秒節点に入力時刻と端点を追加して、フレーム中央時刻へ補間する。`model-time-offset` も適用し、範囲不足はエラーにする。ZIP `.info` に入力、モデル、基準周波数、補正テーブルを記録する。標準出力には有効モデルと入力ファイルを短く表示する。

## 適用範囲と検証

これは一次電離層と非分散中性大気の補正である。Faraday 回転、高次電離層、振幅シンチレーション、屈折による方向変化は含めない。100 MHz未満などで一次項だけで十分かは観測条件に依存する。FFT 内は既存と同じフレーム中央時刻の近似で、フレーム内の急速な位相変動・チャネル内分散を完全に戻す処理ではない。

IONEX/GIM の取得・読込、IRI、GPT/VMF、気象解析のレイトレーシングは未実装。入力層と遅延・位相計算を分離してあり、これらの provider を追加できる。VTEC の単一薄層は水平一様の近似で、IPP 位置の空間勾配を含めない。実用上は天体方向の STEC を優先する。

電離層は標準の物理 RF キャリア（`YI_FRINGE_STOP_FREQ_MODE=low`）と通常の帯域グリッドで検証している。実験用 `rotation2`/Hermitian グリッドは拒否する。既存の局名依存の経験的グリッドオフセットを、この大気補正に流用していない。

検証は独立に分散を注入した8-bit実電圧で、正負・時間変動する TEC、LSB/USB、異なる局帯域、通常相関・debug 経路・パルサー折り畳み・位相合成を照合する。元のビンでの残留位相を確認する。XML 定数・CSV 定数・UTC 時系列の出力一致、範囲外入力の停止、標準大気の従来挙動も試験する。Niell は [独立な RTKLIB C の72条件](../reference/atmosphere/validation/rtklib-niell.json) と数値照合する。

実観測の低周波 VLBI における改善量は、実際の TEC・気象データを入れたうえで評価する。説明用の値を実測値として使わない。

出版論文・取得版・DOI・ダウンロード状況は [reference/atmosphere/](../reference/atmosphere/README.md) に保存する。
