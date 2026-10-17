# KATMER: katmanlı levha çözücü

Tarayıcıda çalışan, kurulumsuz bir katmanlı kompozit levha programı. Katmanları, malzemeyi ve yükü girersiniz; program klasik laminasyon teorisiyle gerilmeleri, Tsai-Hill ölçütüyle ilk katman hasarını hesaplar ve levhayı 3B olarak renklendirir. Baretli Kafa YouTube kanalının 12. bölümünde yapay zekaya yazdırıldı.

*A zero-install, single-file browser app for composite laminates: classical lamination theory stresses and Tsai-Hill first-ply failure.*

Canlı sürüm: https://baretlikafa.pages.dev/katmer/

## Kullanım

`index.html` dosyasını bir tarayıcıda açın. Tek dosya; sunucu, paket ya da internet bağlantısı gerekmez.

## Sınama

`sinama/` klasöründeki betikler, `index.html` içindeki hesap çekirdeğini (`<script id="cekirdek">`) Node ile çalıştırıp programdan bağımsız yazılmış Python hesabıyla karşılaştırır. Gereksinim: Node.js, Python 3, `numpy`, `sympy`.

```
cd sinama
python sina.py          # örnek levha, izotrop ve denge sınamaları, rastgele 400 levha; sonuç sina_sonuc.txt
python karsilastir.py   # dört durumu referans.json ile karşılaştırır
python bagimsiz.py      # referans.json dosyasını yeniden üretir
```

## Sınırlar

Eğitim ve deneme programıdır. Yalnız ilk katman hasarına bakar, tek malzemeyle çalışır; sıcaklık gerilmesini, delikleri ve kenar etkilerini hesaba katmaz. Gerçek bir tasarımı yalnız bu programa bakarak yapmayın.

## Lisans

MIT ([LICENSE](LICENSE)).
