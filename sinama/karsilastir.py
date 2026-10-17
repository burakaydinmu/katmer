# Programin cekirdegini (index.html icindeki betik, Node ile) bagimsiz referansla (referans.json) karsilastirir.
import json, subprocess, sys, numpy as np
HTML = sys.argv[1] if len(sys.argv) > 1 else "../index.html"
print("sinanan dosya:", HTML)
ref = json.load(open("referans.json", encoding="utf-8"))
cam = dict(E1=53.78e9, E2=17.93e9, nu12=0.25, G12=8.62e9, Xt=1035e6, Xc=1035e6, Yt=27.6e6, Yc=138e6, S=41.4e6)
t = 0.127e-3
def kat(dizi): return [dict(aci=a, t=t, n=n) for a, n in dizi]
girdi = {
 "1-odev [0/90x10/0], Nx=1": (kat([(0,1),(90,10),(0,1)]), [1,0,0,0,0,0]),
 "2-yari-izotrop [0/45/-45/90]s, Nx=1000": (kat([(a,1) for a in (0,45,-45,90,90,-45,45,0)]), [1000,0,0,0,0,0]),
 "3-simetrik-olmayan [0/90], Nx=1000, Mx=1": (kat([(0,1),(90,1)]), [1000,0,0,1,0,0]),
 "4-aci-kat [30/-30]s, Nx=5000, Nxy=2000, My=0.5": (kat([(30,1),(-30,1),(-30,1),(30,1)]), [5000,0,2000,0,0.5,0]),
}
adlar = list(girdi)
g = dict(hesap=[dict(cam, katmanlar=girdi[a][0], yuk=girdi[a][1]) for a in adlar])
out = json.loads(subprocess.run(["node", "kopru.js", HTML], input=json.dumps(g), capture_output=True, text=True, encoding="utf-8").stdout)
def bagil(a, b):
    a, b = np.array(a, float), np.array(b, float); olcek = np.abs(b).max()
    return 0.0 if olcek == 0 else float(np.abs(a - b).max() / olcek)
enb = 0
for ad, p in zip(adlar, out["hesap"]):
    r = ref[ad]; assert not p["hatalar"], p["hatalar"]
    f = {k: bagil(p[k], r[k]) for k in ("A", "D")}
    f["B"] = float(np.abs(np.array(p["B"]) - np.array(r["B"])).max() / (np.abs(np.array(r["A"])).max() * r["h_mm"] * 1e-3))
    f["eps0+kappa"] = max(bagil(p["eps0"], r["e0"]) if max(map(abs, r["e0"])) > 0 else 0, bagil(np.array(p["kappa"]) * p["h"], np.array(r["kap"]) * p["h"] + 0) if max(map(abs, r["kap"])) > 1e-12 else 0)
    f["yuk carpani"] = abs(p["carpan"] - r["yuk_carpani"]) / r["yuk_carpani"]
    print(ad); print("   program: carpan %.6f | referans: %.6f | kritik katman (program) %d %s, mod: %s" % (p["carpan"], r["yuk_carpani"], p["kritikKatman"] + 1, p["kritikNokta"]["yuz"], p["mod"]))
    print("   bagil farklar:", {k: "%.1e" % v for k, v in f.items()}); enb = max(enb, *f.values())
p = out["hesap"][0]
print("odev: sigma_x 0-kat %.4f Pa (Ansys 2022: 1490,5) | 90-kat %.4f Pa (Ansys 2022: 489,29) | ilk hasar %.2f N/m (ChatGPT: 56 415,09)" % (p["katmanlar"][0]["ust"]["sxy"][0], p["katmanlar"][1]["ust"]["sxy"][0], p["carpan"]))
print("EN BUYUK BAGIL FARK: %.1e" % enb)
