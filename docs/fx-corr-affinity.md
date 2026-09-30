# CPU affinity の比較（2026-09-30）

`--cpu` の worker 数制限と CPU 固定を切り分けるため、`--no-affinity` を追加した。

```bash
time target/release/yi-corr --sc test.xml --raw raw --cor cor-pinned --cpu 10
time target/release/yi-corr --sc test.xml --raw raw --cor cor-unpinned --cpu 10 --no-affinity
```

通常は計算 worker を1論理 CPUずつ、reader と出力処理を I/O 用 CPU に固定する。
固定する CPU は番号順に選ぶため、物理コアを均等に選ぶ方式ではない。
SMT では別の論理 CPU 番号でも同じ物理コアを共有する。

`--no-affinity` の動作:

- worker 数は通常と同じ `min(N, OS許可CPU数)−1`、最小1。省略時の N は物理コア数。
- reader、Rayon worker、出力処理の CPU 固定を解除する。
- affinity 設定ファイルと `YI_READER_CORE` は読み込まない。
- `$HOME/.yi-corr` のプロセス間 CPU 予約を使わない。同時実行の CPU 重複は防がない。
- 外部から継承した `taskset` / cpuset 等の CPU マスクは尊重する。
  このモードの `--cpu` は worker 数を制限し、使用可能な CPU 番号を N 個に絞らない。

## 手元での測定

Core i7-12700H（14物理 / 20論理 CPU）、release + `target-cpu=native`。
`data/test/yi-corr.xml` の5秒分、FFT8192、`--cpu 10`（9 worker）、
`--chunk-frames 32768 --pipeline-depth 2`。
affinity 設定ファイルと reader 指定なし、他の yi-corr との CPU 予約競合なし。
初回を各モード1回除き、固定→解除、解除→固定、固定→解除の順で各3回。
両モードとも `/proc/PID/task/TID/status` を20 ms間隔で読み取ってマスクを確認した。

| 指標（各3回の中央値） | 固定あり | 固定解除 |
|---|---:|---:|
| プロセス全体の実行時間 [s] | 4.777 | 3.640 |
| 計算区間 [s] | 4.611 | 3.476 |
| consumer 入力待ち [s] | 0.103 | 0.119 |
| input read 呼び出し区間 [s] | 1.479 | 1.469 |

全体で約24%短縮。入力 read 区間はほぼ同じで、差は主に計算区間に現れた。
入力 read 区間は計算と重なり、全体時間には単純加算できない。
これはファイルキャッシュが効く条件であり、HDD / USB / RAID の I/O 待ちは再現していない。
ユーザー環境で同じ短縮率になることは示していない。

このホストでは固定時の worker は CPU 0〜8、reader は CPU 9。
0/1、2/3、4/5、6/7、8/9 がそれぞれ SMT の同一物理コアである。
そのため9 worker は物理5コアに集中し、reader と worker も物理コアを共有していた。
この配置は時間差の一因と考えられるが、SMT の影響だけを独立に測ったものではない。
解除時はメイン、9 worker、reader の全スレッドが `0-19` のマスクを継承した。

別途、OS マスクを `0-3` に制限した `--cpu 10 --no-affinity` は3 workerへ丸められ、
全スレッドのマスクが `0-3` のままであることを確認した。

## 正しさの確認

- 実データの固定 / 解除計8回で、全 `.cor` の SHA-256 が一致。
- `cargo test --offline --bin yi-corr`: 53件成功。
- 追加の実行テスト: 両バイナリで `--cpu 1/3`、小チャンク、USB 並列読み込みを比較し、
  `.cor` / phased `.raw` がバイト一致。無効な affinity 設定と reader 指定を無視し、
  CPU 予約ファイルを変更しないことも確認。
- `cargo build --offline --release --bin yi-corr --bin yi-phasedarray`: 成功。

測定ログと CPU マスクは `/tmp/fx-corr-affinity-check-x5dd7t5a/` に保存。
