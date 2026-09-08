# 復号処理の追加高速化（2026-09-09）

20計算スレッドで5秒分の相関処理を **4.602 → 3.579秒** に短縮した。同時期の変更前後比較で **1.286倍**、実行時間 **22.2%減**。以前報告した4.65秒とは別の反復測定であり、過去の最良値と比較した数字ではない。

## 変更

- src/utils.rs: 2-bit入力の共通レコーダー配列（恒等、VSREC、ペア入替、同じバイト内ビット順を持つバイト並べ替え）を256項の小さな変換表で復号する。USBとLSBの偶数・奇数開始用の計3表、約12 KiB。恒等配列で使っていた計4 MiBの16-bit表を置き換え、VSRECの汎用32-bitシャッフルもこの表へ統合した。任意のバイト間ビット並べ替えは従来の汎用経路を使う。
- src/main.rs: 内部の遅延窓では入力rawを直接参照し、物理word境界に整列した完全なFFT窓はFFT入力へ直接復号する。その他の窓では既存scratchを使い、出力のうち欠損境界だけをゼロ埋めする。

遅延モデル、FFTサイズ、フリンジ回転の位置、整数遅延の符号、LSBの絶対サンプル偶奇は変更していない。GICO3のcomplex方式との差や長時間の三次状位相残差を解決する変更ではない。

## 実測

環境：Core i7-12700H、20論理CPU、Rust 1.94.0、release + target-cpu=native。
入力：data/test/yi-corr.xml、data/test/ のYAMAGU32/USUDA64 raw。
1.024 Gsample/s、FFT=8192、重なり3200ビン、処理長5秒。
--cpu 20、YI_TIMING_SAMPLE_STRIDE=256、chunk/pipelineは既存の自動設定。

変更前後を交互に各4回実行し、各版の初回を除く3回の中央値を使用した。入力がページキャッシュに載った条件であり、温度や他プロセスの負荷により時間は変動する。

| 項目 | 変更前 | 変更後 |
|---|---:|---:|
| プロセス実時間 [s] | 4.602 | 3.579 |
| 計算部分 [s] | 4.265 | 3.237 |
| user+system CPU時間 [s] | 79.118 | 58.938 |
| 全20論理CPUに対する平均使用率 | 86.0% | 81.9% |

初回のサンプル計測では復号のCPU時間が0.103000 → 0.020241秒、計算中の割合が34.7% → 9.0%になった。FFTは変更後のサンプル計算時間の61.9%、積算は29.1%。これは抽出フレームの累積CPU時間であり、実時間ではない。

## 検証

- 4回×3種類の出力、計12個の.corファイルが変更前後でファイル全体のバイト一致。
- USB/LSB、開始偶奇、6種類のシャッフル、端数サンプル、不完全な入力wordを独立したビット抽出基準と比較。
- 1/2/4/8-bit入力で、整数遅延、前後のゼロ埋め、未整列窓、再利用バッファを全体復号から切り出す基準と比較。
- cargo test --offline: 35単体テスト×2バイナリと5結合テストが通過。
- cargo build --offline --release --bin yi-corr --bin yi-phasedarray 完了。

新しい実行ファイルは target/release/yi-corr と target/release/yi-phasedarray。

    YI_TIMING_SAMPLE_STRIDE=256 target/release/yi-corr --sc data/test/yi-corr.xml --raw data/test --cor /tmp/fx-corr-fast --cpu 20

測定ログ、各回の出力、timings.json、comparison.json、binaries.json は /tmp/fx-corr-validation/decode-20260909/ に保存。

SHA-256：

- 変更前: 40431b32bfff0953b27d6c5c5c31e185d5b60f9e4c7068e5d7e6a53828209f2f
- 変更後: ab240a24bb9360796d1a8b8f3d0b7bacf3310f7fb04154adb380dd8171173e56
