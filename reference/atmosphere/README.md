# 大気補正の出版文献と検証資料

使用方法・符号規約は [docs/atmosphere.md](../../docs/atmosphere.md)。観測 raw はこのディレクトリへコピーしない。

`literature/` に PDF を保存し、`manifest.json` と各 PDF の隣の JSON に出版社・掲載誌・DOI・実際のダウンロード URL・取得版・SHA-256 を記録する。本文は PDF 冒頭の著者・題名と照合した。arXiv の著者版と出版社版は区別する。

| 論文 | 出版社・用途 | 取得状況 |
|---|---|---|
| [Landskron & Böhm (2018), VMF3/GPT3](https://doi.org/10.1007/s00190-017-1066-2) | Springer, *Journal of Geodesy* 92, 349–360。改訂静水圧天頂遅延式、地球中性大気・マッピングの整理 | [出版社版 PDF](literature/landskron-boehm-2018-vmf3-gpt3.pdf) |
| [Yan et al. (2021), Revisiting the light time correction…](https://doi.org/10.1007/s00190-021-01498-5) | Springer, *Journal of Geodesy* 95, 48。一次電離層の位相前進と群遅延の符号 | [出版社版 PDF](literature/yan-2021-ionosphere-propagation.pdf) |
| [Hofmeister & Böhm (2017), Application of ray-traced tropospheric slant delays…](https://doi.org/10.1007/s00190-017-1000-7) | Springer, *Journal of Geodesy* 91, 945–964。VLBI の中性大気の視線遅延、天頂遅延・マッピングの限界 | [出版社版 PDF](literature/hofmeister-boehm-2017-neutral-atmosphere.pdf) |
| [Soja et al. (2015), Tropospheric delay determination by Kalman filtering VLBI data](https://doi.org/10.1186/s40623-015-0293-0) | SpringerOpen, *Earth, Planets and Space* 67, 144。湿潤遅延の時間変動と観測入力 | [出版社版 PDF](literature/soja-2015-vlbi-wet-delay.pdf) |
| [de Gasperin et al. (2018), The effect of the ionosphere on ultra-low-frequency radio-interferometric observations](https://doi.org/10.1051/0004-6361/201833012) | EDP Sciences, *Astronomy & Astrophysics* 615, A179。低周波の一次・高次電離層、Faraday 回転、シンチレーション | [著者公開版 PDF](literature/degasperin-2018-ionosphere.pdf)。出版社の PDF は HTTP 403。保存物は arXiv v1、出版社版ではない |
| [Niell (1996), Global mapping functions for the atmosphere delay at radio wavelengths](https://doi.org/10.1029/95JB03048) | AGU/Wiley, *JGR: Solid Earth* 101(B2), 3227–3246。実装の Niell 係数・季節・高度項 | **PDF 未取得**。出版社・ADS・NRAO は HTTP 403、旧著者公開 URL は404。書誌情報と失敗記録を保存 |
| [Saastamoinen (1972), Atmospheric Correction for the Troposphere and Stratosphere in Radio Ranging Satellites](https://doi.org/10.1029/GM015p0247) | AGU/Wiley, *Geophysical Monograph Series* 15, 247–251。気象値からの大気遅延の原典 | **PDF 未取得**。出版社・ADS は HTTP 403。書誌情報と失敗記録を保存 |

VMF3/GPT3、レイトレーシング、Kalman 推定、高次電離層は拡張のための参考文献であり、今回の実装には含めていない。中性大気は静水圧＋湿潤の積分遅延をモデル化し、成層圏を別の任意遅延として加えていない。

Niell 原論文の式(4)/(5)には既知の訂正があるため、実装は訂正済みの正規化連分数を使う。`validation/rtklib-niell.{csv,json}` は T. Takasu の [RTKLIB 原著者 C 実装](https://github.com/tomojitakasu/RTKLIB/blob/master/src/rtkcmn.c) と独立に照合した72条件（南北半球、季節、緯度補間、標高、仰角）。参照ソースの SHA-256 も保存し、Rust の単体試験で数値照合する。C ソース自体は `/tmp` に置く。

再取得・既存 PDF の検証:

```bash
python3 tools/download_atmosphere_references.py
cd reference/atmosphere
sha256sum -c SHA256SUMS
```

ダウンロードできなかった文献を成功扱いにしない。文献を追加で入手した場合も PDF の題名・DOI・取得版を確認してから記録を更新する。
