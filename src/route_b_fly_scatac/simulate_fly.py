"""Fly-format simulated data: fragments.tsv.gz + singlecell.csv per run, and a GTF."""
import gzip, os, sys, numpy as np, pandas as pd
out, scenario = sys.argv[1], sys.argv[2]
rng = np.random.default_rng(int(sys.argv[3]) if len(sys.argv) > 3 else (7 if scenario == "effect" else 8))
os.makedirs(f"{out}/adult_brain", exist_ok=True)
prog = {"cnc": "cnc Keap1 Gclc Gclm GstD1 GstD2 GstE1 GstS1 Zw Pgd Men",
        "trx": "Trxr-1 Trxr-2 Trx-2 Jafrac1 Jafrac2 Prx3 Prx5 PHGPx",
        "anti": "Sod1 Sod2 Cat Eip71CD SelR Ogg1", "iron": "Fer1HCH Fer2LCH Irp1A Tsf1 Mvl",
        "mito": "ND-75 ND-49 ND-30 ND-42 SdhA SdhB UQCR-C1 UQCR-C2 Cyt-c-p COX4 COX5A ATPsynbeta blw",
        "upr": "Hsc70-3 Xbp1 crc PEK Ire1 Atf6 Pdi ERp60 Calr", "hs": "Hsp70Aa Hsp83 Hsp68 Hsp27 Hsp26 Hsc70-4 DnaJ-1 Hsf",
        "dna": "p53 mei-41 tefu lok spn-A mus304 Xrp1 rpr hid", "ieg": "Hr38 sr kay Jra puc CrebA"}
markers = {"neuron": "elav nSyb brp fne", "glia": "repo", "surface": "moody Mdr65", "astro": "alrm Eaat1", "ens": "Eaat2", "cortex": "wrapper"}
prog["lic"] = "Gem dup Cdc6 Orc1 Mcm2 Mcm3 dpa Mcm5 Mcm6 Mcm7 CycE E2f1 Myc CycD"
prog["hripd"] = "Ku80 Irbp Cht8 CG3344 CG4734 Fanci lok rad50 Xrp1 p53"
named = sorted({g for v in list(prog.values()) + list(markers.values()) for g in v.split()} | {"Gem", "dup"})
genes = named + [f"CG{i:05d}" for i in range(1500)]
chroms = ["2L", "2R", "3L", "3R", "X", "4"]
rows, pos = [], {}
for i, g in enumerate(rng.permutation(genes)):
    c = chroms[i % 6] if i % 6 != 5 else "2L"
    start = 10000 + (i // 6) * 6000 + int(rng.integers(0, 1000))
    pos[g] = (c, start, start + 2000)
    rows.append(f'{c}\tFlyBase\tgene\t{start + 1}\t{start + 2000}\t.\t+\t.\tgene_id "FBgn{i:07d}"; gene_name "{g}"; gene_biotype "protein_coding";\n')
with gzip.open(f"{out}/Drosophila_melanogaster.gtf.gz", "wt") as h:
    h.write("#!genome-build BDGP6\n"); h.writelines(rows)
gidx = {g: i for i, g in enumerate(genes)}
base = rng.gamma(0.7, 1.0, len(genes))
types = {"neuron": 0.72, "glia": 0.08, "astro": 0.07, "cortex": 0.06, "surface": 0.04, "ens": 0.03}
tf = {t: rng.lognormal(0, 0.3, len(genes)) for t in types}
for t, ms in markers.items():
    for g in ms.split():
        tf[t if t in tf else "neuron"][gidx[g]] *= 8
for g in markers["glia"].split():
    for t in ("astro", "surface", "ens"): tf[t][gidx[g]] *= 8
n_cells = 1800
recovered = n_cells + 200
mult = recovered / 1000 * 0.008
frag_lines, barcodes = [], []
truth = []
tnames = list(types); tp = np.array(list(types.values()))
for c in range(n_cells):
    t = rng.choice(tnames, p=tp)
    share = types[t]
    dbl = rng.random() < mult * share / (1 - mult + mult * share)
    if dbl: copies = 4
    elif t == "surface": copies = int(rng.choice([8, 16], p=[0.6, 0.4]))
    elif t in ("neuron", "astro", "cortex") and rng.random() < 0.06: copies = 4
    else: copies = 2
    poly = copies > 2 and not dbl
    w = base * tf[t] * np.exp(rng.normal(0, 0.15, len(genes)))
    for name, gs in prog.items():
        act = rng.normal(0, 0.15)
        if scenario == "effect" and poly and name in ("cnc", "lic"): act += 0.6
        for g in gs.split(): w[gidx[g]] *= np.exp(act)
    w = w / w.sum()
    per_copy = int(1500 * np.exp(rng.normal(0, 0.35)))
    bc = "".join(rng.choice(list("ACGT"), 16)) + "-1"
    barcodes.append(bc)
    truth.append((bc, t, copies, dbl))
    for k in range(copies):
        n_in = int(per_copy * 0.6); n_bg = per_copy - n_in
        g = rng.choice(len(genes), n_in, p=w)
        ch = np.array([pos[genes[i]][0] for i in g]); st = np.array([pos[genes[i]][1] for i in g]) + rng.integers(-1000, 2000, n_in)
        bgc = rng.choice(chroms[:4], n_bg); bgs = rng.integers(20000, 10000 + 260 * 6000, n_bg)
        ch = np.r_[ch, bgc]; st = np.r_[st, bgs]
        en = st + np.clip(60 + rng.exponential(110, len(st)).astype(int), 50, 900)
        df = pd.DataFrame({"c": ch, "s": st, "e": en}).sort_values(["c", "s"])
        keep = np.ones(len(df), bool)
        for cc, part in df.groupby("c"):
            s, e = part.s.to_numpy(), part.e.to_numpy()
            ok = np.r_[True, s[1:] >= np.maximum.accumulate(e)[:-1]]
            keep[df.index.get_indexer(part.index)] = ok
        df = df[keep]
        frag_lines.append(df.assign(b=bc))
    for _ in range(rng.poisson(3)):
        cc = rng.choice(chroms[:4]); s0 = int(rng.integers(20000, 1500000))
        d = int(rng.integers(3, 6))
        frag_lines.append(pd.DataFrame({"c": [cc] * d, "s": s0 + rng.integers(0, 30, d), "e": s0 + 150 + rng.integers(0, 30, d), "b": bc}))
f = pd.concat(frag_lines, ignore_index=True)
f["c"] = "chr" + f["c"]          # test chrom-name handling
f["n"] = 1
order = {f"chr{c}": i for i, c in enumerate(chroms)}
f = f.sort_values(["c", "s"], key=lambda s: s.map(order) if s.name == "c" else s)
with gzip.open(f"{out}/adult_brain/fragments.tsv.gz", "wt") as h:
    h.write("# id=sim\n")
    f[["c", "s", "e", "b", "n"]].to_csv(h, sep="\t", header=False, index=False)
noise = ["".join(rng.choice(list("ACGT"), 16)) + "-1" for _ in range(200)]
sc = pd.DataFrame({"barcode": ["NO_BARCODE"] + barcodes + noise, "total": 1,
                   "is__cell_barcode": [0] + [1] * len(barcodes) + [0] * 200})
sc.to_csv(f"{out}/adult_brain/singlecell.csv", index=False)
pd.DataFrame(truth, columns=["barcode", "type", "copies", "doublet"]).to_csv(f"{out}/truth.csv", index=False)
print("fragments:", len(f))
