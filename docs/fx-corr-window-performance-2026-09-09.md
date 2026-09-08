# 整数遅延窓のコピー削減（2026-09-09）

RustFFT 6.4.0 / RealFFT 3.5.0 を継続使用する。FFT ライブラリの追加や切り替えは行っていない。

## 採用した変更

src/main.rs の decode_shifted_frame_from_chunk で、整数遅延によって入力窓が物理32-bit wordの途中から始まる場合も、中央の完全なword群をFFT入力へ直接復号する。先頭と末尾の部分wordだけscratchを経由する。これまでの窓全体のscratch復号とコピーを除去した。遅延・位相補正とFFTサイズは従来どおり。

## 実測

入力は data/test/yi-corr.xml と data/test/ の raw、5秒分、FFT 8192、20計算スレッド。YI_TIMING_SAMPLE_STRIDE=256。Core i7-12700H、release + target-cpu=native。
変更前後をAB/BA交互に各4回、初回を除く3回の中央値で比較した。

| 項目 | 変更前 | 変更後 |
|---|---:|---:|
| 実時間 [s] | 3.8529 | 3.6962 |
| 計算部分 [s] | 3.439 | 3.289 |
| user+system CPU時間 [s] | 59.9204 | 58.5999 |

今回の中央値では実時間4.07%減、CPU時間2.20%減。変更後が遅い組もあり、温度や他の負荷によるばらつきを含む小幅な改善である。過去の3.579秒とは測定時期が異なるため直接比較しない。入力はページキャッシュが効く条件。

別途試したRayon foldによるバッファ再利用は3.7358 → 3.8020秒で改善を示さず、元に戻した。

## 検証

- 4回×3種類、12個の .cor 出力すべてが変更前後でファイル全体のバイト一致。
- cargo test --offline: 35単体テスト×2バイナリ、5結合テストが通過。復号窓テストは1/2/4/8-bit、USB/LSB、開始偶奇、シャッフル、整数遅延、端数raw、境界ゼロ埋めを比較する。
- cargo build --offline --release --bin yi-corr --bin yi-phasedarray 完了。
- 長時間・長基線の三次状位相残差を解決したことを示す検証ではない。

ログ・出力・集計: /tmp/fx-corr-validation/window-20260909/
不採用の並列化案のログ: /tmp/fx-corr-validation/fold-20260909/

実行ファイルSHA-256:

- 変更前: ab240a24bb9360796d1a8b8f3d0b7bacf3310f7fb04154adb380dd8171173e56
- 変更後: df43b23c25b4cb5f17ec4ee3cbbd921ca70b62a48d0b16ffa90eca1620fd1d64
