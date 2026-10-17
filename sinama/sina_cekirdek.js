// index.html içindeki hesap çekirdeğini (<script id="cekirdek">) dosyadan çeker ve çalıştırır.
// Girdi (stdin, JSON): { hesap: [g, g, ...], sayi: ["53,78", ...] }
// Çıktı (stdout, JSON): { hesap: [sonuç...], sayi: [sonuç...] }
// Böylece sınanan kod, tarayıcıda çalışan kodun AYNISIDIR; kopyası değil.
const fs = require('fs');
const path = require('path');
// Hangi dosya sınanacak: LEVHA_DOSYA ortam değişkeni (verilmezse ../index.html).
const html = fs.readFileSync(path.join(__dirname, process.env.LEVHA_DOSYA || '../index.html'), 'utf8');
const m = html.match(/<script id="cekirdek">([\s\S]*?)<\/script>/);
if (!m) { console.error('cekirdek betiği bulunamadı'); process.exit(2); }
const C = new Function(m[1] + '\nreturn LevhaCekirdek;')();

let veri = '';
process.stdin.setEncoding('utf8');
process.stdin.on('data', d => veri += d);
process.stdin.on('end', () => {
  const g = JSON.parse(veri);
  // Infinity JSON'da yok: null'a döner; Python tarafı bunu bilir.
  const cikti = {
    hesap: (g.hesap || []).map(x => C.hesapla(x)),
    sayi: (g.sayi || []).map(x => C.sayiOku(x))
  };
  process.stdout.write(JSON.stringify(cikti));
});
