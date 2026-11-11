# 帯域幅・サンプリング周波数が異なる局の相関（3.12.0）

`yi-corr` は各局の実数サンプリング周波数を XML の局別 terminal から読み、共通 RF 帯域を相関する。帯域幅は各局の `speed/2`。一方が他方を包含する場合、共通部分は狭帯域側の全帯域になる。

例: 広帯域局 6600–7112 MHz（512 MHz）、狭帯域局 6600–6616 MHz（16 MHz）。

| 項目 | 広帯域局 | 狭帯域局 |
|---|---:|---:|
| サンプリング周波数 | 1024 Msamples/s | 32 Msamples/s |
| 入力帯域幅 | 512 MHz | 16 MHz |
| FFT 点数 | 65536 | 2048 |
| FFT が扱う時間 | 64 µs | 64 µs |
| 周波数刻み | 15.625 kHz | 15.625 kHz |

両局の同じ時刻・同じ RF 周波数のビンを取り出し、6600–6616 MHz の16 MHz分を積算する。広帯域側も元の記録速度で FFT する。相関積分に使うのは共通 RF のビンである。

## XML と実行

局ごとに terminal を分ける。`stream/fft` と `--fft` は**狭帯域側の FFT 点数**を指定する。広帯域側は記録速度の比で自動的に増える。

```xml
<terminal name="wide">
  <speed>1024000000</speed><channel>1</channel><bit>2</bit>
  <level>-1.5,-0.5,0.5,1.5</level>
</terminal>
<terminal name="narrow">
  <speed>32000000</speed><channel>1</channel><bit>2</bit>
  <level>-1.5,-0.5,0.5,1.5</level>
</terminal>
```

既存の `<station key="A">` の `<terminal>` を `wide`、局 B を `narrow` として参照する。level、shuffle、sideband は実際の記録方式に合わせる。

```xml
<stream>
  <frequency>6600000000</frequency>
  <fft>2048</fft><output>1</output>
  <special key="A"><rotation>0</rotation><sideband>USB</sideband></special>
  <special key="B"><rotation>0</rotation><sideband>USB</sideband></special>
</stream>
```

```bash
target/release/yi-corr --sc observation.xml --raw /observations/raw --cor cor
```

狭帯域が広帯域の中央にある場合も同じ方法で扱う。例えば狭帯域が6608–6624 MHzなら、上の例で局 B の `rotation` を `8000000` Hzにする。出力の RF 原点は自動的に6608 MHzになる。局の順序を交換しても、出力は狭帯域側の RF グリッドになる。

## 遅延・位相・振幅

各局で `N_i/f_s,i` を等しくし、`Δf=f_s,i/N_i` を等しくする。読み出すサンプル数・整数遅延・小数遅延は局ごとの `f_s,i` を使う。共通の時間差は秒で保持する。

周波数原点が異なると、ネイティブ FFT の位相は各 FFT 窓の開始時刻に依存する。既知の LO 差と実際の読み出し位置から位相原点を共通化する。整数サンプル時刻で剰余を計算し、長い skip による大きな位相の丸め誤差を抑える。

ACF は各局の `N_i²`、XCF は `N_1 N_2` で FFT 振幅を正規化する。片側だけ FFT 点数が増えることによる振幅の32倍・1024倍などの誤差を防ぐ。再量子化は行わない。

`--delay` 等のサンプル単位は狭帯域側の出力サンプリング基準。XML clock の delay/rate は従来どおり秒・秒/秒。Hz で指定する user rate の基準は XML の `stream/frequency`。

## 出力と対応範囲

- この例の `.cor` は sampling=32 MHz、FFT=2048、RF=6600 MHz の16 MHz帯域を持つ。実数 FFT の Nyquist 端点を重複して格納せず、1024個の複素ビンを保存する。
- 標準出力と ZIP `.info` に両局の RF 範囲・帯域幅・記録速度・FFT 点数、共通帯域・時間幅を表示する。
- 部分重複では狭帯域側の出力周波数軸を保ち、XCF の非重複部分はゼロ。共通 RF がなければエラーになる。
- サンプリング周波数の比が整数で、LO 差が選択した `Δf` の整数倍になる組み合わせに対応する。整数倍でなければ FFT 点数を変える等の設定が必要。任意の有理比のリサンプリングや fractional-frequency mixer は未実装であり、未対応の入力はエラーにする。
- 異速度モードは XML の RF 設定でビンを選び、従来の局名依存の経験的オフセットを適用しない。`YI_ANT2_GRID_OFFSET` は0を必要とする。
- 通常相関、debug のマッピング経路、パルサー折り畳み、USB 並列読み出し、入力バッファに対応。一次電離層補正とも組み合わせられる。
- 異なる記録速度の入力は `yi-corr` で対応。`yi-phasedarray` の raw 合成は現在、同じサンプリング周波数の入力を必要とする。実験用 `rotation2` と carrier 切替も異速度モードでは対象外。

検証には合成電圧のみを使い、帯域の下端・中央・上端での包含、局順の交換、LSB/USB、2-bit/8-bit、非ゼロ clock、skip、チャンク境界、電離層位相、帯域外の強い信号を試験する。位相と振幅を独立な注入値と照合する。第三者の観測データは使用・公開しない。
