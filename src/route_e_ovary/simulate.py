import gzip, os, sys, numpy as np, scipy.sparse as sp
from scipy.io import mmwrite
out, scen = sys.argv[1], sys.argv[2]
rng = np.random.default_rng(11)
M = {
 "fc": "eya br Fas3 CadN Fas2 Sox14 mirr pnt".split(),
 "germ": "vas nos bam orb osk Rbp9 corolla".split(),
 "escort": "ptc hh tj Wnt4 bnl".split(),
 "muscle": "Mhc Mlc1 up zip Act57B".split(),
 "mit": "stg CycB polo aurB Det pav PCNA CycA Incenp".split(),
 "germ_early": "bam nos vas Sxl".split(),
 "germ_late": "orb osk Rbp9 corolla".split(),
 "early": "CadN Fas2 ct SPARC".split(),
 "mid": "br mirr mid Yp1 Yp2 Yp3".split(),
 "late": "Sox14 Vha16-1 Pvf1 dpp cv-2".split(),
 "amp": "Vm26Aa Vm26Ab Vm32E Vm34Ca Cp36 Cp38 Cp15 Cp16 dec-1 psd".split(),
}
prog = ("GstD1 GstD2 GstD3 GstD5 GstD9 GstE1 GstE3 GstE6 GstE7 GstE9 GstE12 GstE14 GstS1 GstO1 GstO2 GstO3 GstT1 GstZ1 GstZ2 Mgstl "
 "Gclc Gclm Gss1 Gss2 Gr1 Grx1 Grx1t PHGPx GILT1 Trxr-1 Trxr-2 Trx-2 TrxT dhd Jafrac1 Jafrac2 Prx3 Prx5 Prx2540-1 Prx2540-2 Srx "
 "Sod1 Sod2 Sod3 Cat Zw Pgd Men Men-b Idh Idh3a Nmdmc Tkt Taldo cnc Keap1 Cyp6a2 Cyp6a8 Ugt35a Ero1L Duox Nox Cyp4e3 Mav Uro Aox1 Pxd Pxt "
 "ND-75 ND-49 ND-30 ND-42 SdhA SdhB UQCR-C1 UQCR-C2 Cyt-c-p COX4 COX5A ATPsynbeta blw Xbp1 Hsc70-3 Ire1 PEK crc Pdi ERp60 Calr Sec61alpha "
 "Sec61beta Srp54 RpL3 RpL4 RpL5 RpL7 RpL8 RpL10 RpL11 RpL13 RpS2 RpS3 RpS4 RpS6 RpS8 RpS9 RpS11 RpS15 Gem dup Cdc6 Orc1 Mcm2 Mcm3 dpa "
 "Mcm5 Mcm6 Mcm7 CycE E2f1 Myc p53 Ku80 Irbp lok mei-41 tefu rad50 Xrp1 spn-A mus304 Hsp70Aa Hsp83 Hsp68 Hsp27 Hsp26 Hsp23 DnaJ-1 Hsf "
 "Hr38 sr kay Jra puc CrebA Act5C alphaTub84B RpL32 Gapdh1 Tbp eIF4A Pgk Ef1alpha48D").split()
named = sorted(set(sum(M.values(), []) + prog))
genes = named + [f"CG{i}" for i in range(4000)] + ["mt:CoI", "mt:CoII", "mt:ND1"]
gi = {g: i for i, g in enumerate(genes)}
base = np.exp(rng.normal(1.2, 1.1, len(genes)))
types = [("mitotic FC", 900, ["fc", "mit", "early"]), ("endocycling FC", 1100, ["fc", "mid", "late"]),
         ("late endocycling FC", 500, ["fc", "late", "amp"]),
         ("germarium germline", 350, ["germ", "mit", "germ_early"]),
         ("nurse cells", 500, ["germ", "germ_late"]),
         ("escort", 350, ["escort"]), ("muscle", 300, ["muscle"])]
gst = [g for g in prog if g.startswith("Gst")]
nadph = "Zw Pgd Men Idh Taldo Tkt".split()
rows, labels = [], []
for lib in ("GSM4381234_ovary1", "GSM4381235_ovary2"):
    for name, n, sets in types:
        mu = base.copy()
        for s in sets:
            for g in M[s]:
                mu[gi[g]] *= 25
        if "endocycling" in name:
            for g in ("Yp1", "Yp2", "Yp3"):
                mu[gi[g]] *= 60                      # yolk genes dominate the library in these cells
            if scen == "effect":
                for g in gst + nadph:
                    mu[gi[g]] *= 2.3
                mu[gi["Jafrac1"]] *= 0.45            # falls in follicle cells
                for g in ("Hsc70-3", "Sec61beta", "Sec61alpha", "Pdi", "Ero1L"):
                    mu[gi[g]] *= 3.0                 # ER load: follicle cells ONLY
        if name == "nurse cells" and scen == "effect":
            mu[gi["Jafrac1"]] *= 0.45                # the SAME fall in the germline: tracks endocycling
            mu[gi["Jafrac2"]] *= 1.02
        depth = 4000 if "endocycling" in name else 2500
        mu = mu / mu.sum() * depth
        for _ in range(n):
            c = rng.poisson(mu * rng.uniform(0.6, 1.5))
            rows.append(sp.csr_matrix(c)); labels.append((lib, name))
    d = os.path.join(out, "geo", "series", "GSE146nnn", "GSE146040", "suppl")
X = sp.vstack(rows).tocsr()
os.makedirs(d, exist_ok=True)
for lib in sorted({l for l, _ in labels}):
    m = np.array([l == lib for l, _ in labels])
    sub = X[np.flatnonzero(m)]
    import io as _io
    b = _io.BytesIO(); mmwrite(b, sp.coo_matrix(sub.T.astype(int)), field="integer")
    open(f"{d}/{lib}_matrix.mtx.gz", "wb").write(gzip.compress(b.getvalue()))
    with gzip.open(f"{d}/{lib}_barcodes.tsv.gz", "wt") as h:
        for i in range(int(m.sum())): h.write(f"CELL{i:05d}-1\n")
    with gzip.open(f"{d}/{lib}_features.tsv.gz", "wt") as h:
        for j, g in enumerate(genes): h.write(f"FBgn{2000000+j}\t{g}\tGene Expression\n")
print("ok", scen, X.shape)
