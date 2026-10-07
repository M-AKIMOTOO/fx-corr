# CPU affinity と同時実行（3.10.0）

`yi-corr --cpu 6` は、この yi-corr の計算5スレッドと I/O を合わせて
6論理 CPU に割り当てる。Linux ではメインスレッドも同じ CPU 群に制限し、
後から生成する補助スレッドはそのマスクを継承する。計算 worker と reader
は、その中の各 CPU に固定する。他プロセスの CPU マスクや OS 全体の設定を
変更する処理はない。

3.10.0 では次のように選ぶ:

- まず物理コアごとに1論理 CPU を選び、必要なら同じ割り当て内の SMT を使う。
- 選んだ物理コアの SMT sibling 全体を `$HOME/.yi-corr` に予約する。
  使用 CPU は N 個だが、予約ファイルの CPU ID は N より多くなることがある。
- 既に起動したプロセスの割り当てを変更せず、後のプロセスに空いている
  物理コアを割り当てる。古い予約が一方の sibling だけでも、その物理コアを避ける。
- `yi-phasedarray` と `--cpu` 省略時にも同じ予約を使う。
- OS から見える論理 CPU を最低2個、予約外に残す。確保できなければ
  先行プロセスの CPU を奪わず、説明付きのエラーで終了する。

`--cpu 1` では計算と I/O が1 CPUを共有する。
`taskset` / cpuset や affinity ファイルの許可範囲を尊重する。
SMT の判定は Linux の `thread_siblings_list` を用いる。
Unix の予約は同じ HOME を使う本プログラム間の協調であり、他のアプリケーション
を排除する OS 専用 CPU ではない。Linux 以外で topology が得られなければ
論理 CPU 単位の予約となる。Unix 以外ではこの予約を行わない。

例えばこの Core i7-12700H（14物理 / 20論理 CPU）では、設定ファイルなしの
2本の `--cpu 6` は次の割り当てになる:

| 起動順 | 計算 CPU | I/O CPU | 予約する CPU ID |
|---|---|---:|---|
| 1本目 | 0,2,4,6,8 | 10 | 0–11（SMT sibling を含む） |
| 2本目 | 12,13,14,15,16 | 17 | 12–17 |

18,19 は予約外に残る。CPU番号や P/E コア構成はマシンごとに異なる。
この例では2本目が E コアなので、両方の処理速度が同じになる保証はない。
メモリ帯域、キャッシュ、CPU の電力/温度上限、同じ RAID の読み出しも共有する。
affinity 修正だけで同時実行の速度低下すべてを防げるわけではない。

CPU 固定とスレッド数を切り分けるには `--no-affinity` を使える。
worker 数は同じだが、CPU 予約と固定を解除するため同時実行のコア重複は防がない。

## 時間差起動の検証（2026-10-07）

合成2-bit RAW、268435456 sample/s、16秒、FFT8192、出力1 Hzを使用。
各局は1 GiB、同じ内容のファイルを hard link で参照し、キャッシュが効く条件で
実行した。`--cpu 6 --chunk-frames 32768 --pipeline-depth 2` で、2本目を
約0.6秒後に開始。3.9.0 と修正版をそれぞれ単独・時間差起動で測った。
これは各条件1回の確認であり、RAID の連続読み出し性能は測っていない。

| 全体実行時間 [s] | 3.9.0 | 修正版 |
|---|---:|---:|
| 単独 | 5.568 | 4.002 |
| 同時実行の1本目 | 7.114 | 4.441 |
| 同時実行の2本目 | 7.235 | 8.344 |

修正版の2本目は上の表の E コア配置なので遅い。1本目も共有資源の影響で
単独より約11%長くかかり、速度低下が完全になくなるわけではない。
旧版の1本目の増加はこの確認では約28%だった。

40 ms 間隔の `/proc/PID/task/TID/status` 確認で、修正版は処理中の全スレッドが
割り当てた6 CPU以内で動作した。2本目開始後も1本目の割り当ては変わらず、
両プロセスは物理コアを共有しなかった。監視側プロセスのマスクは変化しなかった。
開始前と処理終了後のガード復元は測定対象から除いた。
全6実行の ACF/XCF `.cor` は SHA-256 が一致し、終了時に予約が解放された。
測定ログ・マスク・集計は `/tmp/fx-corr-cpu-isolation-l4gp8ecp/` に保存。

なお検証前に PATH の `/home/akimoto/.cargo/bin/yi-corr` は 3.7.0 だった。
ビルドしただけでは普段の `yi-corr` は更新されないため、インストール先の
`yi-corr --version` も確認する必要がある。

## 以前の固定/解除比較（2026-09-30、3.8.0）

`--cpu` の worker 数制限と CPU 固定を切り分けるため、`--no-affinity` を追加した。

```bash
time target/release/yi-corr --sc test.xml --raw raw --cor cor-pinned --cpu 10
time target/release/yi-corr --sc test.xml --raw raw --cor cor-unpinned --cpu 10 --no-affinity
```

この時点では計算 worker を1論理 CPUずつ、reader と出力処理を I/O 用 CPU に固定していた。
固定する CPU は番号順に選んでおり、物理コアを均等に選ぶ方式ではなかった。
SMT では別の論理 CPU 番号でも同じ物理コアを共有する。

`--no-affinity` の動作:

- worker 数は通常と同じ `min(N, OS許可CPU数)−1`、最小1。省略時の N は物理コア数。
- reader、Rayon worker、出力処理の CPU 固定を解除する。
- affinity 設定ファイルと `YI_READER_CORE` は読み込まない。
- `$HOME/.yi-corr` のプロセス間 CPU 予約を使わない。同時実行の CPU 重複は防がない。
- 外部から継承した `taskset` / cpuset 等の CPU マスクは尊重する。
  このモードの `--cpu` は worker 数を制限し、使用可能な CPU 番号を N 個に絞らない。

### 手元での測定

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

### 正しさの確認

- 実データの固定 / 解除計8回で、全 `.cor` の SHA-256 が一致。
- `cargo test --offline --bin yi-corr`: 53件成功。
- 追加の実行テスト: 両バイナリで `--cpu 1/3`、小チャンク、USB 並列読み込みを比較し、
  `.cor` / phased `.raw` がバイト一致。無効な affinity 設定と reader 指定を無視し、
  CPU 予約ファイルを変更しないことも確認。
- `cargo build --offline --release --bin yi-corr --bin yi-phasedarray`: 成功。

測定ログと CPU マスクは `/tmp/fx-corr-affinity-check-x5dd7t5a/` に保存。
