// Verilen HTML dosyasindaki <script id="cekirdek"> betigini calistirir (sinama klasorunun kendi koprusu).
// Kullanim: node kopru.js <html-yolu>   (stdin: {hesap:[...]}  stdout: {hesap:[...]})
const fs = require('fs');
const html = fs.readFileSync(process.argv[2], 'utf8');
const m = html.match(/<script id="cekirdek">([\s\S]*?)<\/script>/);
if (!m) { console.error('cekirdek betigi bulunamadi'); process.exit(2); }
const C = new Function(m[1] + '\nreturn LevhaCekirdek;')();
let veri = '';
process.stdin.setEncoding('utf8');
process.stdin.on('data', d => veri += d);
process.stdin.on('end', () => process.stdout.write(JSON.stringify({ hesap: JSON.parse(veri).hesap.map(x => C.hesapla(x)) })));
