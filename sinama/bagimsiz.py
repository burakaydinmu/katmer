# Bolum 12 - programdan BAGIMSIZ sinama hesabi (numpy). Program yazilirken ajan bu dosyayi gormedi.
# Klasik laminasyon teorisi + Tsai-Hill; dort sinama durumu. Cikti: referans.json
import json, numpy as np

def Qmat(E1, E2, nu12, G12):
    nu21 = nu12 * E2 / E1
    d = 1 - nu12 * nu21
    return np.array([[E1/d, nu12*E2/d, 0], [nu12*E2/d, E2/d, 0], [0, 0, G12]])

def T_sigma(th):
    c, s = np.cos(np.radians(th)), np.sin(np.radians(th))
    return np.array([[c*c, s*s, 2*c*s], [s*s, c*c, -2*c*s], [-c*s, c*s, c*c - s*s]])

R = np.diag([1, 1, 2.0])

def Qbar(Q, th):
    T = T_sigma(th)
    return np.linalg.inv(T) @ Q @ R @ T @ np.linalg.inv(R)

def coz(mal, katlar, yuk):
    """katlar: [(aci, kalinlik_m)] ustten alta; yuk: [Nx,Ny,Nxy,Mx,My,Mxy]"""
    Q = Qmat(mal["E1"], mal["E2"], mal["nu12"], mal["G12"])
    h = sum(t for _, t in katlar)
    z = [-h/2]
    for _, t in katlar: z.append(z[-1] + t)
    A = np.zeros((3, 3)); B = np.zeros((3, 3)); D = np.zeros((3, 3))
    for k, (th, t) in enumerate(katlar):
        Qb = Qbar(Q, th)
        A += Qb * (z[k+1] - z[k]); B += Qb * (z[k+1]**2 - z[k]**2) / 2; D += Qb * (z[k+1]**3 - z[k]**3) / 3
    ABD = np.block([[A, B], [B, D]])
    ek = np.linalg.solve(ABD, np.array(yuk, float))
    e0, kap = ek[:3], ek[3:]
    sonuc = []; fi_max = 0; ilk = None
    for k, (th, t) in enumerate(katlar):
        for yuz, zz in (("alt", z[k]), ("ust", z[k+1])):
            sxy = Qbar(Q, th) @ (e0 + zz * kap)
            s12 = T_sigma(th) @ sxy
            X = mal["Xt"] if s12[0] >= 0 else mal["Xc"]
            Y = mal["Yt"] if s12[1] >= 0 else mal["Yc"]
            fi = s12[0]**2/X**2 - s12[0]*s12[1]/X**2 + s12[1]**2/Y**2 + s12[2]**2/mal["S"]**2
            sonuc.append(dict(kat=k+1, aci=th, yuz=yuz, z_mm=zz*1e3, sxy=sxy.tolist(), s12=s12.tolist(), FI=fi))
            if fi > fi_max: fi_max, ilk = fi, (k+1, yuz)
    return dict(h_mm=h*1e3, A=A.tolist(), B=B.tolist(), D=D.tolist(), e0=e0.tolist(), kap=kap.tolist(),
                noktalar=sonuc, FI_max=fi_max, yuk_carpani=1/np.sqrt(fi_max), ilk_hasar=ilk)

cam = dict(E1=53.78e9, E2=17.93e9, nu12=0.25, G12=8.62e9, Xt=1035e6, Xc=1035e6, Yt=27.6e6, Yc=138e6, S=41.4e6)
t = 0.127e-3
durumlar = {
 "1-odev [0/90x10/0], Nx=1": (cam, [(0, t)] + [(90, t)]*10 + [(0, t)], [1, 0, 0, 0, 0, 0]),
 "2-yari-izotrop [0/45/-45/90]s, Nx=1000": (cam, [(a, t) for a in (0, 45, -45, 90, 90, -45, 45, 0)], [1000, 0, 0, 0, 0, 0]),
 "3-simetrik-olmayan [0/90], Nx=1000, Mx=1": (cam, [(0, t), (90, t)], [1000, 0, 0, 1, 0, 0]),
 "4-aci-kat [30/-30]s, Nx=5000, Nxy=2000, My=0.5": (cam, [(30, t), (-30, t), (-30, t), (30, t)], [5000, 0, 2000, 0, 0.5, 0]),
}
cikti = {ad: coz(*v) for ad, v in durumlar.items()}
json.dump(cikti, open("referans.json", "w", encoding="utf-8"), indent=1)
for ad, r in cikti.items():
    print(ad); print("  A11 %.6g  A12 %.6g  A66 %.6g  B11 %.6g  D11 %.6g" % (r["A"][0][0], r["A"][0][1], r["A"][2][2], r["B"][0][0], r["D"][0][0]))
    print("  e0", ["%.6g" % v for v in r["e0"]], "kap", ["%.6g" % v for v in r["kap"]])
    print("  FI_max %.6g  yuk carpani %.6g  ilk hasar %s" % (r["FI_max"], r["yuk_carpani"], r["ilk_hasar"]))
r = cikti["1-odev [0/90x10/0], Nx=1"]
print("odev: sigma_x 0-kat %.4f Pa, 90-kat %.4f Pa" % (r["noktalar"][0]["sxy"][0], r["noktalar"][2]["sxy"][0]))
q = cikti["2-yari-izotrop [0/45/-45/90]s, Nx=1000"]["A"]
print("yari-izotrop denetimi: A11-A22 = %.3g, A66-(A11-A12)/2 = %.3g" % (q[0][0]-q[1][1], q[2][2]-(q[0][0]-q[0][1])/2))
