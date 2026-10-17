# -*- coding: utf-8 -*-
"""index.html hesap çekirdeğinin bağımsız sınaması.

Çalıştırma:  python sina.py        (sinama klasöründe; Node, numpy ve sympy gerekir)

Program (JavaScript) ile bu dosya (Python) aynı fiziği FARKLI yoldan hesaplar:

  programda                               burada
  --------------------------------------  ---------------------------------------------
  Q: kapalı formül E/(1-v12 v21)          Q = inv(S)  (esneklik matrisinin tersi)
  Q-bar: sin/cos'lu kapalı formüller      Q-bar = inv(T) Q R T inv(R)  (matris çarpımı)
  A,B,D: (z1-z0), (z1²-z0²)/2, (z1³-z0³)/3  A,B,D: katman katman Gauss-Legendre sayısal integral
  6x6: ölçekli sistem, kendi Gauss elemesi   6x6: numpy.linalg.solve, ölçeksiz
  malzeme ekseni: açık dönüşüm formülü    malzeme ekseni: T matrisi ile çarpım
  en büyük gösterge: yüzler + işaret      en büyük gösterge: katman başına 4001 noktada
     değişim noktası                         kaba kuvvet tarama

Ayrıca kod karşılaştırması olmayan sınamalar:
  - denge: programın verdiği gerilmeler kalınlık boyunca toplanınca girilen N ve M çıkıyor mu
  - izotrop malzeme: sigma = N/h, sigma = 6M/h², D = E h³ / 12(1-v²)  (ders kitabı formülleri)
  - örnek levha: sympy ile kesirli (yuvarlamasız) kapalı çözüm
  - ilk hasar yükü uygulanınca gösterge 1 oluyor mu
"""
import json
import os
import subprocess
import sys

import numpy as np
import sympy as sp

sys.stdout.reconfigure(encoding="utf-8")
KLASOR = os.path.dirname(os.path.abspath(__file__))
satirlar = []
basarisiz = []


def yaz(*a):
    s = " ".join(str(x) for x in a)
    satirlar.append(s)
    print(s)


def kontrol(ad, kosul, ayrinti=""):
    yaz(("  GEÇTİ  " if kosul else "  KALDI  ") + ad + (("  | " + ayrinti) if ayrinti else ""))
    if not kosul:
        basarisiz.append(ad)


def programi_calistir(hesap=(), sayi=()):
    girdi = json.dumps({"hesap": list(hesap), "sayi": list(sayi)})
    r = subprocess.run(["node", os.path.join(KLASOR, "sina_cekirdek.js")], input=girdi.encode("utf-8"),
                       capture_output=True, check=True)
    return json.loads(r.stdout.decode("utf-8"))


# ------------------------------------------------------------------ bağımsız hesap
def T_matris(aci):
    r = np.deg2rad(aci)
    c, s = np.cos(r), np.sin(r)
    return np.array([[c * c, s * s, 2 * c * s], [s * s, c * c, -2 * c * s], [-c * s, c * s, c * c - s * s]])


R = np.diag([1.0, 1.0, 2.0])
GL_X, GL_W = np.polynomial.legendre.leggauss(3)  # 5. dereceye kadar tam


def tsai_hill(s1, s2, t12, g):
    X = g["Xt"] if s1 >= 0 else g["Xc"]
    Y = g["Yt"] if s2 >= 0 else g["Yc"]
    return (s1 / X) ** 2 - s1 * s2 / X ** 2 + (s2 / Y) ** 2 + (t12 / g["S"]) ** 2


def bagimsiz(g, tarama=4001):
    S = np.array([[1 / g["E1"], -g["nu12"] / g["E1"], 0], [-g["nu12"] / g["E1"], 1 / g["E2"], 0], [0, 0, 1 / g["G12"]]])
    Q = np.linalg.inv(S)
    kal = [k["t"] * k["n"] for k in g["katmanlar"]]
    h = sum(kal)
    z = np.concatenate([[-h / 2], -h / 2 + np.cumsum(kal)])
    Qb = []
    for k in g["katmanlar"]:
        T = T_matris(k["aci"])
        Qb.append(np.linalg.inv(T) @ Q @ R @ T @ np.linalg.inv(R))
    A = np.zeros((3, 3)); B = np.zeros((3, 3)); D = np.zeros((3, 3))
    for i in range(len(kal)):
        orta, yari = (z[i] + z[i + 1]) / 2, (z[i + 1] - z[i]) / 2
        for x, w in zip(GL_X, GL_W):
            zz = orta + yari * x
            A += Qb[i] * w * yari
            B += Qb[i] * zz * w * yari
            D += Qb[i] * zz * zz * w * yari
    K = np.block([[A, B], [B, D]])
    cz = np.linalg.solve(K, np.array(g["yuk"], dtype=float))
    e0, kp = cz[:3], cz[3:]
    kat = []
    for i, k in enumerate(g["katmanlar"]):
        T = T_matris(k["aci"])

        def nokta(zz):
            sxy = Qb[i] @ (e0 + zz * kp)
            s12 = T @ sxy
            return sxy, s12, tsai_hill(s12[0], s12[1], s12[2], g)

        tara = max(nokta(zz)[2] for zz in np.linspace(z[i], z[i + 1], tarama))
        kat.append({"ust": nokta(z[i]), "alt": nokta(z[i + 1]), "enb": tara, "Qb": Qb[i]})
    enb = max(k["enb"] for k in kat)
    return {"Q": Q, "h": h, "z": z, "A": A, "B": B, "D": D, "e0": e0, "kp": kp, "kat": kat, "enb": enb,
            "carpan": (1 / np.sqrt(enb)) if enb > 0 else None}


def bagil_fark(a, b, olcek=None):
    a, b = np.asarray(a, float), np.asarray(b, float)
    o = olcek if olcek is not None else max(np.max(np.abs(a)), np.max(np.abs(b)), 1e-300)
    return float(np.max(np.abs(a - b)) / o)


def karsilastir(g, p):
    """Program sonucu p ile bağımsız hesabı karşılaştırır; en büyük bağıl farkları döndürür."""
    b = bagimsiz(g)
    f = {}
    h = b["h"]
    Ae = np.max(np.abs(b["A"]))
    f["Q"] = bagil_fark(p["Q"], b["Q"])
    f["A"] = bagil_fark(p["A"], b["A"], Ae)
    f["B"] = bagil_fark(p["B"], b["B"], Ae * h)
    f["D"] = bagil_fark(p["D"], b["D"], Ae * h * h)
    eo = max(np.max(np.abs(b["e0"])), np.max(np.abs(b["kp"])) * h, 1e-300)
    f["gerinim"] = bagil_fark(p["eps0"], b["e0"], eo)
    f["egrilik"] = bagil_fark(np.array(p["kappa"]) * h, b["kp"] * h, eo)
    so = max(max(np.max(np.abs(k["ust"][0])), np.max(np.abs(k["alt"][0])), np.max(np.abs(k["ust"][1])), np.max(np.abs(k["alt"][1]))) for k in b["kat"])
    so = max(so, 1e-300)
    f["Qbar"] = max(bagil_fark(pk["Qb"], bk["Qb"], np.max(np.abs(b["Q"]))) for pk, bk in zip(p["katmanlar"], b["kat"]))
    f["gerilme_xy"] = max(max(bagil_fark(pk[y]["sxy"], bk[y][0], so) for y in ("ust", "alt")) for pk, bk in zip(p["katmanlar"], b["kat"]))
    f["gerilme_12"] = max(max(bagil_fark(pk[y]["s12"], bk[y][1], so) for y in ("ust", "alt")) for pk, bk in zip(p["katmanlar"], b["kat"]))
    if b["enb"] > 0:
        f["gosterge_yuz"] = max(max(abs(pk[y]["th"]["gosterge"] - bk[y][2]) for y in ("ust", "alt")) for pk, bk in zip(p["katmanlar"], b["kat"])) / b["enb"]
        # programın katman en büyüğü taramadakinden küçük olmamalı, çok da büyük olmamalı
        f["gosterge_enb_eksik"] = max(max(0.0, bk["enb"] - pk["gosterge"]) for pk, bk in zip(p["katmanlar"], b["kat"])) / b["enb"]
        f["gosterge_enb_fazla"] = max(max(0.0, pk["gosterge"] - bk["enb"]) for pk, bk in zip(p["katmanlar"], b["kat"])) / b["enb"]
        f["carpan"] = abs(p["carpan"] - b["carpan"]) / b["carpan"]
        # denge: programın gerilmelerinden N ve M
        N = np.zeros(3); M = np.zeros(3)
        for pk in p["katmanlar"]:
            s0, s1 = np.array(pk["ust"]["sxy"]), np.array(pk["alt"]["sxy"])
            z0, z1 = pk["zUst"], pk["zAlt"]
            sm, zm = (s0 + s1) / 2, (z0 + z1) / 2
            N += (s0 + s1) / 2 * (z1 - z0)
            M += (s0 * z0 + 4 * sm * zm + s1 * z1) / 6 * (z1 - z0)  # Simpson: ikinci derece için tam
        yuk = np.array(g["yuk"], float)
        yo = max(np.max(np.abs(yuk[:3])), np.max(np.abs(yuk[3:])) / h)
        f["denge_N"] = float(np.max(np.abs(N - yuk[:3])) / yo)
        f["denge_M"] = float(np.max(np.abs(M - yuk[3:])) / h / yo)
    return f, b


# ------------------------------------------------------------------ 1. örnek levha
GPa, MPa, mm = 1e9, 1e6, 1e-3
ORNEK = {"E1": 53.78 * GPa, "E2": 17.93 * GPa, "nu12": 0.25, "G12": 8.62 * GPa,
         "Xt": 1035 * MPa, "Xc": 1035 * MPa, "Yt": 27.6 * MPa, "Yc": 138 * MPa, "S": 41.4 * MPa,
         "katmanlar": [{"aci": 0, "t": 0.127 * mm, "n": 1}, {"aci": 90, "t": 0.127 * mm, "n": 10}, {"aci": 0, "t": 0.127 * mm, "n": 1}],
         "yuk": [1, 0, 0, 0, 0, 0]}

yaz("=" * 78)
yaz("1) ÖRNEK LEVHA  [0 / 90x10 / 0], E-cam/epoksi, Nx = 1 N/m")
yaz("=" * 78)
p = programi_calistir([ORNEK])["hesap"][0]
kontrol("program hata vermedi", p["hatalar"] == [])
f, b = karsilastir(ORNEK, p)
for ad, v in f.items():
    kontrol("program = bağımsız numpy hesabı: %-18s" % ad, v < 1e-9, "bağıl fark %.2e" % v)

# sympy: kesirli sayılarla kapalı çözüm (simetrik çapraz dizilim, yalnız Nx)
r = sp.Rational
E1, E2, nu, G = r(5378, 100) * 10 ** 9, r(1793, 100) * 10 ** 9, r(1, 4), r(862, 100) * 10 ** 9
t = r(127, 1000) / 1000
payda = 1 - nu * nu * E2 / E1
Q11, Q22, Q12 = E1 / payda, E2 / payda, nu * E2 / payda
A11 = 2 * t * Q11 + 10 * t * Q22
A22 = 2 * t * Q22 + 10 * t * Q11
A12 = 12 * t * Q12
ex, ey = sp.symbols("ex ey")
coz = sp.solve([A11 * ex + A12 * ey - 1, A12 * ex + A22 * ey], [ex, ey])
ex_, ey_ = coz[ex], coz[ey]
s1_0, s2_0 = Q11 * ex_ + Q12 * ey_, Q12 * ex_ + Q22 * ey_          # 0° katman: 1 = x
s1_90, s2_90 = Q11 * ey_ + Q12 * ex_, Q12 * ey_ + Q22 * ex_        # 90° katman: 1 = y, 2 = x
Xt, Yt, Yc = r(1035) * 10 ** 6, r(276, 10) * 10 ** 6, r(138) * 10 ** 6
th0 = (s1_0 / Xt) ** 2 - s1_0 * s2_0 / Xt ** 2 + (s2_0 / Yt) ** 2
th90 = (s1_90 / Xt) ** 2 - s1_90 * s2_90 / Xt ** 2 + (s2_90 / Yt) ** 2   # s1_90 < 0 ama Xc = Xt
carpan_sp = 1 / sp.sqrt(sp.Max(th0, th90))
yaz("  sympy kapalı çözüm: eps_x = %.10e, eps_y = %.10e" % (float(ex_), float(ey_)))
yaz("  sympy: 0° katman sigma1 = %.6f Pa, sigma2 = %.6f Pa" % (float(s1_0), float(s2_0)))
yaz("  sympy: 90° katman sigma1 = %.6f Pa, sigma2 = %.6f Pa" % (float(s1_90), float(s2_90)))
yaz("  sympy: gösterge 0° = %.6e, 90° = %.6e, yük çarpanı = %.4f" % (float(th0), float(th90), float(carpan_sp)))
kontrol("sympy: eps_x", abs(p["eps0"][0] - float(ex_)) / float(ex_) < 1e-11)
kontrol("sympy: eps_y", abs(p["eps0"][1] - float(ey_)) / abs(float(ey_)) < 1e-11)
kontrol("sympy: 0° katman sigma1, sigma2", bagil_fark(p["katmanlar"][0]["ust"]["s12"][:2], [float(s1_0), float(s2_0)]) < 1e-11)
kontrol("sympy: 90° katman sigma1, sigma2", bagil_fark(p["katmanlar"][1]["ust"]["s12"][:2], [float(s1_90), float(s2_90)]) < 1e-11)
kontrol("sympy: yük çarpanı", abs(p["carpan"] - float(carpan_sp)) / float(carpan_sp) < 1e-11)
kontrol("elle denge: sigma_x(0°)·2t + sigma_x(90°)·10t = 1 N/m", abs(float(s1_0 * 2 * t + s2_90 * 10 * t) - 1) < 1e-12)
kontrol("kritik katman 2. katman (90°), biçim matris", p["kritikKatman"] == 1 and p["mod"]["tur"] == "matris", p["mod"]["ad"] + " / " + p["mod"]["ayrinti"])

yaz("")
yaz("  PROGRAMIN ÖRNEK SONUÇLARI")
yaz("  h = %.4f mm" % (p["h"] * 1e3))
yaz("  Q (GPa): Q11 = %.4f, Q12 = %.4f, Q22 = %.4f, Q66 = %.4f" % (p["Q"][0][0] / GPa, p["Q"][0][1] / GPa, p["Q"][1][1] / GPa, p["Q"][2][2] / GPa))
yaz("  A (N/m): A11 = %.5e, A12 = %.5e, A22 = %.5e, A66 = %.5e" % (p["A"][0][0], p["A"][0][1], p["A"][1][1], p["A"][2][2]))
yaz("  B: en büyük |Bij| = %.3e (simetrik dizilim: sıfır olmalı)" % np.max(np.abs(p["B"])))
yaz("  D (N·m): D11 = %.5f, D12 = %.5f, D22 = %.5f, D66 = %.5f" % (p["D"][0][0], p["D"][0][1], p["D"][1][1], p["D"][2][2]))
yaz("  eps0 = [%.5e, %.5e, %.3e], kappa = %s" % (p["eps0"][0], p["eps0"][1], p["eps0"][2], p["kappa"]))
for k in p["katmanlar"]:
    yaz("  katman %d (%g°): sigma_x = %.3f, sigma_y = %.3f Pa | sigma1 = %.3f, sigma2 = %.3f Pa | gösterge = %.5e | çarpan = %.1f | %s"
        % (k["no"], k["aci"], k["ust"]["sxy"][0], k["ust"]["sxy"][1], k["ust"]["s12"][0], k["ust"]["s12"][1], k["gosterge"], k["carpan"], k["mod"]["ad"]))
yaz("  İLK HASAR: çarpan = %.2f  ->  Nx = %.2f N/m = %.3f kN/m (Nx/h = %.3f MPa), katman %d, %s (%s)"
    % (p["carpan"], p["carpan"], p["carpan"] / 1e3, p["carpan"] / p["h"] / 1e6, p["kritikKatman"] + 1, p["mod"]["ad"], p["mod"]["ayrinti"]))

# ilk hasar yükü uygulanınca gösterge 1 olmalı
g2 = dict(ORNEK); g2["yuk"] = [p["carpan"], 0, 0, 0, 0, 0]
p2 = programi_calistir([g2])["hesap"][0]
kontrol("ilk hasar yükü girilince en büyük gösterge = 1", abs(p2["gostergeEnb"] - 1) < 1e-12, "%.15f" % p2["gostergeEnb"])
kontrol("o yükte 90° katmanda sigma2 = Yt'ye çok yakın (27,6 MPa)", abs(p2["katmanlar"][1]["ust"]["s12"][1] / MPa - 27.6) < 0.05,
        "%.4f MPa" % (p2["katmanlar"][1]["ust"]["s12"][1] / MPa))

# ------------------------------------------------------------------ 2. izotrop malzeme
yaz("")
yaz("=" * 78)
yaz("2) İZOTROP MALZEME (ders kitabı formülleri)  E = 70 GPa, v = 0,3, h = 2 mm, karışık açılar")
yaz("=" * 78)
E, v, hh = 70e9, 0.3, 2e-3
IZO = {"E1": E, "E2": E, "nu12": v, "G12": E / (2 * (1 + v)), "Xt": 3e8, "Xc": 3e8, "Yt": 3e8, "Yc": 3e8, "S": 3e8 / np.sqrt(3),
       "katmanlar": [{"aci": 37, "t": 0.25e-3, "n": 2}, {"aci": -12.5, "t": 0.5e-3, "n": 1}, {"aci": 90, "t": 0.2e-3, "n": 5}]}
pN = programi_calistir([dict(IZO, yuk=[1000, 0, 0, 0, 0, 0]), dict(IZO, yuk=[0, 0, 0, 10, 0, 0])])["hesap"]
sx = [k[y]["sxy"][0] for k in pN[0]["katmanlar"] for y in ("ust", "alt")]
kontrol("Nx altında her yerde sigma_x = N/h", max(abs(s - 1000 / hh) for s in sx) / (1000 / hh) < 1e-12, "N/h = %.1f Pa" % (1000 / hh))
kontrol("D11 = E h³ / 12(1-v²)", abs(pN[1]["D"][0][0] - E * hh ** 3 / 12 / (1 - v * v)) / pN[1]["D"][0][0] < 1e-12, "%.6f N·m" % pN[1]["D"][0][0])
kontrol("A11 = E h / (1-v²)", abs(pN[1]["A"][0][0] - E * hh / (1 - v * v)) / pN[1]["A"][0][0] < 1e-12)
kontrol("B = 0", np.max(np.abs(pN[1]["B"])) == 0)
ust, alt = pN[1]["katmanlar"][0]["ust"]["sxy"][0], pN[1]["katmanlar"][-1]["alt"]["sxy"][0]
kontrol("Mx altında dış yüzlerde sigma_x = ∓6M/h²", abs(alt - 6 * 10 / hh ** 2) / alt < 1e-12 and abs(ust + 6 * 10 / hh ** 2) / alt < 1e-12,
        "alt %.4e, üst %.4e Pa" % (alt, ust))
kontrol("Mx altında kappa_x = M/D(1-v²)... = 12M/(E h³)", abs(pN[1]["kappa"][0] - 12 * 10 / (E * hh ** 3)) / pN[1]["kappa"][0] < 1e-12)
# izotropta Tsai-Hill (X=Y, S=X/√3) = von Mises: tek eksenli çekmede çarpan = X/sigma
kontrol("izotropta çarpan = X / (N/h)  (Tsai-Hill -> von Mises)", abs(pN[0]["carpan"] - 3e8 / (1000 / hh)) / pN[0]["carpan"] < 1e-12, "%.4f" % pN[0]["carpan"])

# ------------------------------------------------------------------ 3. rastgele levhalar
yaz("")
yaz("=" * 78)
yaz("3) RASTGELE 400 LEVHA (simetrik olmayan dizilimler, 6 yük bileşeni birden, Xt ≠ Xc, Yt ≠ Yc)")
yaz("=" * 78)
rng = np.random.default_rng(12)
durumlar = []
for _ in range(400):
    E1 = rng.uniform(20, 250) * GPa
    E2 = rng.uniform(3, 20) * GPa
    g = {"E1": E1, "E2": E2, "nu12": rng.uniform(0.05, 0.4), "G12": rng.uniform(2, 10) * GPa,
         "Xt": rng.uniform(500, 2500) * MPa, "Xc": rng.uniform(300, 1500) * MPa, "Yt": rng.uniform(20, 80) * MPa,
         "Yc": rng.uniform(100, 250) * MPa, "S": rng.uniform(30, 120) * MPa,
         "katmanlar": [{"aci": float(rng.choice([0, 90, 45, -45, 30, -60, 180, 270, round(rng.uniform(-180, 180), 2)])),
                        "t": float(rng.uniform(0.05, 0.5)) * mm, "n": int(rng.integers(1, 6))} for _ in range(int(rng.integers(1, 9)))],
         "yuk": [float(x) for x in np.concatenate([rng.normal(0, 1e5, 3), rng.normal(0, 60, 3)]) * (rng.random(6) < 0.75)]}
    if not any(g["yuk"]):
        g["yuk"][0] = 1e4
    durumlar.append(g)
sonuclar = programi_calistir(durumlar)["hesap"]
enb_fark, ara_kritik = {}, 0
for g, p in zip(durumlar, sonuclar):
    assert p["hatalar"] == [], p["hatalar"]
    f, b = karsilastir(g, p)
    for ad, vv in f.items():
        enb_fark[ad] = max(enb_fark.get(ad, 0), vv)
    if p["kritikNokta"]["yuz"] == "ara":
        ara_kritik += 1
for ad, vv in enb_fark.items():
    # 'fazla': tarama 4001 noktalı olduğu için sivri tepeyi biraz kaçırabilir; onun eşiği gevşek
    esik = 1e-5 if ad == "gosterge_enb_fazla" else 1e-8
    kontrol("400 levhada en büyük bağıl fark: %-20s" % ad, vv < esik, "%.2e (eşik %.0e)" % (vv, esik))
yaz("  (en büyük göstergenin katman içinde, işaret değişim noktasında çıktığı levha sayısı: %d)" % ara_kritik)

# ------------------------------------------------------------------ 3b. katman içi tepe
# Yukarıdaki 400 levhada en büyük gösterge hep bir yüzde çıkabilir; o zaman programdaki
# "işaret değişim noktası" dalı hiç sınanmamış olur. Burada o dalı zorlayan, kasıtlı tuhaf
# dayanımlı (Xt çok küçük, Xc çok büyük) levhalar üretilir ve kaba kuvvet taramayla karşılaştırılır.
yaz("")
yaz("=" * 78)
yaz("3b) KATMAN İÇİ TEPE — kasıtlı tuhaf dayanımlar (Xt << Xc), eğilme ağırlıklı 3000 levha")
yaz("=" * 78)
rng = np.random.default_rng(7)
durumlar = []
for _ in range(3000):
    g = {"E1": rng.uniform(20, 250) * GPa, "E2": rng.uniform(3, 20) * GPa, "nu12": rng.uniform(0.05, 0.4), "G12": rng.uniform(2, 10) * GPa,
         "Xt": rng.uniform(20, 80) * MPa, "Xc": rng.uniform(1500, 3000) * MPa, "Yt": rng.uniform(300, 600) * MPa,
         "Yc": rng.uniform(300, 600) * MPa, "S": rng.uniform(300, 600) * MPa,
         "katmanlar": [{"aci": float(round(rng.uniform(-90, 90), 1)), "t": float(rng.uniform(0.1, 0.5)) * mm, "n": int(rng.integers(1, 4))}
                       for _ in range(int(rng.integers(1, 4)))],
         "yuk": [float(x) for x in np.concatenate([rng.normal(0, 3e4, 3), rng.normal(0, 30, 3)])]}
    durumlar.append(g)
sonuclar = programi_calistir(durumlar)["hesap"]
icte = [(g, p) for g, p in zip(durumlar, sonuclar) if p["kritikNokta"]["yuz"] == "ara"]
yaz("  en büyük göstergenin katman içinde çıktığı levha sayısı: %d / %d" % (len(icte), len(durumlar)))
kontrol("dal en az bir kez çalıştı", len(icte) > 0)
eksik = fazla = kazanc = 0.0
for g, p in icte:
    b = bagimsiz(g, tarama=40001)
    eksik = max(eksik, max(max(0.0, bk["enb"] - pk["gosterge"]) for pk, bk in zip(p["katmanlar"], b["kat"])) / b["enb"])
    fazla = max(fazla, max(max(0.0, pk["gosterge"] - bk["enb"]) for pk, bk in zip(p["katmanlar"], b["kat"])) / b["enb"])
    yuz_enb = max(max(pk["ust"]["th"]["gosterge"], pk["alt"]["th"]["gosterge"]) for pk in p["katmanlar"])
    kazanc = max(kazanc, p["gostergeEnb"] / yuz_enb - 1)
kontrol("içteki tepe: program taramadan küçük değil", eksik < 1e-9, "%.2e" % eksik)
kontrol("içteki tepe: program taramadan en çok tarama adımı kadar büyük", fazla < 1e-3, "%.2e (40001 noktalı tarama)" % fazla)
yaz("  (bu levhalarda yalnız yüzlere bakılsaydı gösterge en kötü durumda %.1f kat küçük bulunurdu)" % (kazanc + 1))
eksik = 0.0
for g, p in zip(durumlar[:300], sonuclar[:300]):
    f, b = karsilastir(g, p)
    eksik = max(eksik, f["gosterge_enb_eksik"], f["gerilme_12"], f["denge_N"], f["denge_M"])
kontrol("ilk 300 levhada gerilme, denge ve en büyük gösterge", eksik < 1e-8, "%.2e" % eksik)

# ------------------------------------------------------------------ 4. anlamsız girdi
yaz("")
yaz("=" * 78)
yaz("4) ANLAMSIZ GİRDİ — program çökmemeli, hatayı söylemeli")
yaz("=" * 78)


def degis(**kw):
    g = json.loads(json.dumps(ORNEK))
    g.update(kw)
    return g


kotu = [("sıfır kalınlık", degis(katmanlar=[{"aci": 0, "t": 0, "n": 1}])),
        ("boş tablo", degis(katmanlar=[])),
        ("negatif modül", degis(E1=-5e9)),
        ("sıfır modül", degis(G12=0)),
        ("negatif kalınlık", degis(katmanlar=[{"aci": 0, "t": -1e-4, "n": 2}])),
        ("kesirli tabaka adedi", degis(katmanlar=[{"aci": 0, "t": 1e-4, "n": 2.5}])),
        ("sıfır tabaka adedi", degis(katmanlar=[{"aci": 0, "t": 1e-4, "n": 0}])),
        ("fiziksel olmayan Poisson oranı", degis(nu12=1.9)),
        ("sıfır dayanım", degis(Yt=0)),
        ("sayı olmayan yük (null)", degis(yuk=[None, 0, 0, 0, 0, 0]))]
for (ad, _), s in zip(kotu, programi_calistir([g for _, g in kotu])["hesap"]):
    kontrol(ad, len(s["hatalar"]) > 0 and "A" not in s, "; ".join(h["mesaj"] for h in s["hatalar"]))
s0 = programi_calistir([degis(yuk=[0, 0, 0, 0, 0, 0])])["hesap"][0]
kontrol("sıfır yük: hata değil, gösterge 0, çarpan sonsuz", s0["hatalar"] == [] and s0["gostergeEnb"] == 0 and s0["carpan"] is None)

# ------------------------------------------------------------------ 5. sayı okuma
yaz("")
yaz("=" * 78)
yaz("5) SAYI OKUMA — virgül ve nokta")
yaz("=" * 78)
ornek = [("53,78", 53.78), ("53.78", 53.78), (" 0,127 ", 0.127), ("-12,5", -12.5), ("1e3", 1000.0), ("1,5E-2", 0.015), (",5", 0.5), ("7.", 7.0), ("+2", 2.0), ("−3", -3.0)]
kotu_s = ["", "abc", "1.035,5", "1,2,3", "12 mm", "--3", "1e", "NaN", "Infinity"]
rs = programi_calistir(sayi=[a for a, _ in ornek] + kotu_s)["sayi"]
for (a, d), r_ in zip(ornek, rs):
    kontrol("“%s” -> %g" % (a, d), r_["tamam"] and abs(r_["deger"] - d) < 1e-15)
for a, r_ in zip(kotu_s, rs[len(ornek):]):
    kontrol("“%s” reddedildi" % a, not r_["tamam"], r_.get("hata", ""))

yaz("")
yaz("=" * 78)
yaz("SONUÇ: " + ("BÜTÜN SINAMALAR GEÇTİ" if not basarisiz else "%d SINAMA KALDI: %s" % (len(basarisiz), basarisiz)))
yaz("=" * 78)
# Sınanan dosya LEVHA_DOSYA, çıktı dosyası SINA_CIKTI ortam değişkeniyle seçilir (verilmezse ../index.html).
satirlar.insert(0, "Sınanan dosya: " + os.environ.get("LEVHA_DOSYA", "../index.html"))
with open(os.path.join(KLASOR, os.environ.get("SINA_CIKTI", "sina_sonuc.txt")), "w", encoding="utf-8") as fh:
    fh.write("\n".join(satirlar) + "\n")
sys.exit(1 if basarisiz else 0)
