# 入力バッファ再利用と I/O 用 CPU の分離（2026-09-09）

ユーザー提供の Ryzen 9 7950X ログは60秒観測を32.507秒で処理し、consumer入力待ち18.057秒、計算8.979秒、遅延配列生成3.389秒だった。RAM使用率の変動が指摘された。従来は chunk=32768 と表示しながら通常相関では1秒積算全体を読み込み、2局合計約512 MBのraw配列と100万フレーム分の遅延配列をまとめて確保していた。先読み数もCPU数に依存していた。

## 動作

- --cpu N は I/O 用1論理CPUを含む総数。計算は N−1。--cpu 30 なら計算29＋I/O 1。--cpu 1 は共用。利用可能なCPU数を超える指定は上限に丸める。
- 既存のaffinity設定とOSの許可CPU集合を尊重する。Linuxでは読み込みスレッドと出力処理・最終flushをI/O CPUに固定し、Rayon計算スレッドは別の論理CPUへ固定する。SMTの物理コア共有までは排除しない。出力終了後に元のCPUマスクを復元し、子プロセスに1CPU制約を引き継がせない。
- 内部I/OをXML積算時間から分離。自動設定は2局合計16 MiBを目安（最大32768フレーム）、ready queue 2、再利用バッファ4組。CPU機種・個数に依存しない。--chunk-frames と --pipeline-depth で上書きできる。
- rawと遅延配列を固定数のスロットで循環再利用する。各フレームの整数遅延から必要な前後サンプル範囲を求め、物理32-bit word境界で読む。元の積算sector端のゼロ埋めと固定/adaptive read-alignを維持する。
- 1秒積算はチャンク結果を合算して従来と同じ時刻・積算秒数で出力する。RustFFT / RealFFT を継続使用する。
- ファイル全体への POSIX_FADV_WILLNEED ヒントを廃止し、SEQUENTIAL ヒントを維持する。OSのページキャッシュは相関器自身のバッファとは別であり、システム全体のRAM使用量を固定する機能ではない。
- 遅延計算はreader内で行い計算と重なるため、consumerのSynth timing summaryのdelayは0。実作業時間は Input delay preparation として別記する。I/O reader summaryのavgは処理全体時間に対する入力量であり、ディスク単独の転送速度ではない。

## 実データ比較

Core i7-12700H、Rust release + target-cpu=native。data/test/yi-corr.xml の5秒分のYAMAGU32–USUDA64実データを使用し、FFTだけ1024/8192に変更。各条件の前後を交互に4回実行、初回を除く中央値。両版とも --cpu 20 を指定したため、旧版は計算20、新版は計算19＋I/O 1。入力はページキャッシュが効く条件であり、Ryzen実機での60秒データの速度を保証する測定ではない。

| FFT | 実時間 旧→新 [s] | 最大RSS 旧→新 [MiB] |
|---|---:|---:|
| 1024 | 2.979 → 2.738 | 1691.0 → 88.6 |
| 8192 | 3.428 → 3.439 | 1721.1 → 94.8 |

FFT1024は約8%短縮、FFT8192はほぼ同等。いずれもプロセスの最大RSSを約95%削減した。最初の試作ではFFT8192で遅くなったため、小チャンクでの計算ジョブ数を増やして待ちを減らした。

## 検証

- 前後各4回×2つのFFT×3種類、24個の相関出力がファイル全体でバイト一致。
- 1/2/4/8-bit、USB/LSB、開始偶奇、整数遅延、チャンク前後のゼロ埋めを全sector復号基準で比較。
- 3フレーム分割とsector一括処理を、正負のclock rate・adaptive read-align・USB並列読み込み・queue深さ1で比較。1秒積算の相関出力と合成rawが一致。
- producer待機中のconsumer終了とファイルオープン失敗がハングせず終了するテストを追加。
- --cpu 4 実行中の /proc の Cpus_allowed_list を確認し、計算workerが0/1/2、readerが3に固定されていることを確認した（affinity.json）。
- cargo test --offline: 単体38件×2バイナリ、結合6件通過。
- cargo build --offline --release --bin yi-corr --bin yi-phasedarray 完了。

最終測定ログ・各回のRSSと時間・出力・SHA-256・比較結果: /tmp/fx-corr-validation/streaming-final-20260909/
最初の試作の測定: /tmp/fx-corr-validation/streaming-20260909/

Ryzen側では再ビルド後、既存のコマンドの --cpu 30 をそのまま使用できる。ログの CPU allocation で I/O CPU と計算29スレッドを確認する。
