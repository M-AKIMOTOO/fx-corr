# 出力 CPU 切り替えの修正とバッファ容量比較

Ryzen 9 7950Xのユーザーログでは、同じ --cpu 32 の計算時間は9.512〜9.713秒に対し、consumerの入力待ちは10.389〜25.641秒だった。--cpu 5 の計算34.183秒、--cpu 10 の計算16.468秒とも比較すると、計算自体は並列化で短縮している。全体時間が似ていても、CPU数が無効とは言えない。

## 修正

src/main.rs の通常相関では、実際には出力しない中間チャンクでも IoAffinityGuard を使ってI/O CPUへ移動していた。読み込みスレッドが実行中の同じCPUに集計処理まで移し、スケジューリング待ちを生む構造だった。移動を積算完了時の実際の出力に限定した。

ユーザーの条件（100万frame/sector、32768frame/chunk、60sector）では、通常相関の出力用CPU移動を1860回から60回へ減らす。raw出力と最後のflushのI/O affinityは維持する。中間チャンクの集計は引き続きoutput計時に含まれるが、そこではCPUを移動しない。

また、Input read calls のログを追加した。これはreaderの読み込み区間の経過時間で、ファイルシステム・ページキャッシュ・seek・読み込み待ちを含む。遅延配列の準備とqueue待ちは除外する。ディスク装置単独の性能測定ではない。

## バッファ容量の考え方

既定の約88 MiBという値はプロセスのRSSである。rawバッファは16 MiB×4組、そのうちready queueは2組。バッファ容量とMB/s・Gbpsの転送速度は別の量であり、容量を増やして改善するのは主に先読み量や転送速度の揺れを吸収できる範囲である。

同じHDD/RAID上の2ファイルを交互に読む場合、chunkを大きくするとファイル切り替え頻度も下がる。queue深さだけを増やした場合とは効果が異なる。実際の配置は未確認なので、ユーザー環境での改善量は断定しない。

既存オプションで比較可能（FFT=1024、2-bit、2局）：

| 設定 | 2局合計raw/chunk | rawスロット容量（遅延余白除く） |
|---|---:|---:|
| --chunk-frames 32768 --pipeline-depth 2 | 16 MiB | 64 MiB |
| --chunk-frames 65536 --pipeline-depth 4 | 32 MiB | 192 MiB |
| --chunk-frames 131072 --pipeline-depth 4 | 64 MiB | 384 MiB |

プロセスRSSには遅延・FFT・積算配列等が追加される。OSのファイルキャッシュは別。CPU数による自動バッファ増加は導入していない。

## 手元での反復測定

Core i7-12700H、--cpu 20（計算19＋I/O 1）、FFT1024、data/test/yi-corr.xmlの5秒実データ。初回を除く各3回の中央値。4条件の実行順を回ごとに循環させた。入力はページキャッシュが効くため、RyzenのHDD/RAIDでの入力待ちは再現していない。

| 条件 | 実時間 [s] | 最大RSS [MiB] | output [s] |
|---|---:|---:|---:|
| 修正前・標準容量 | 2.681 | 87.8 | 0.023 |
| 修正後・標準容量 | 2.677 | 88.1 | 0.007 |
| 修正後・raw192 MiB | 2.633 | 238.1 | 0.020 |
| 修正後・raw384 MiB | 2.686 | 462.6 | 0.033 |

この測定では全体時間の大幅改善は確認できないため、既定容量は変更していない。大容量条件のoutputにはより大きい再利用配列の管理・中間集計時間等も含まれる。

- 修正前に対し、3条件×4回×3種類の相関出力36個がファイル全体でバイト一致。
- cargo test --offline: 単体38件×2バイナリ、結合6件通過。
- releaseのyi-corr/yi-phasedarrayを更新済み。
- 測定ログ、RSS、時間、出力、比較結果: /tmp/fx-corr-validation/io-buffer-comparison/

Ryzen側で修正版を再ビルドした後の比較コマンド例：

    time yi-corr --sc test.xml --raw raw --cor test-small --cpu 32 --chunk-frames 32768 --pipeline-depth 2
    time yi-corr --sc test.xml --raw raw --cor test-large --cpu 32 --chunk-frames 131072 --pipeline-depth 4

同じキャッシュ条件にできない場合は実行順を交互にし、全体時間だけでなくInput read calls、Input delay preparation、recv-wait、outputを比較する。
