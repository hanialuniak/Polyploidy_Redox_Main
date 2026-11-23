"""Mock GSE153165: series matrix, a count table (FBgn rows, title columns) or per-sample files, and a GTF."""
import gzip, os, sys, tarfile, io, numpy as np, pandas as pd
out, scenario, layout = sys.argv[1], sys.argv[2], sys.argv[3]
rng = np.random.default_rng(3 if scenario == "effect" else 4)
os.makedirs(out, exist_ok=True)
named = ("cnc Keap1 Gclc Gclm GstD1 GstD2 GstE1 GstS1 Zw Pgd Men Trxr-1 Trxr-2 Trx-2 Jafrac1 Jafrac2 Prx3 Prx5 PHGPx "
         "Sod1 Sod2 Cat Eip71CD SelR Ogg1 Fer1HCH Fer2LCH Irp1A Tsf1 Mvl ND-75 ND-49 ND-30 ND-42 SdhA SdhB UQCR-C1 UQCR-C2 "
         "Cyt-c-p COX4 COX5A ATPsynbeta blw Gem dup Cdc6 Orc1 Mcm2 Mcm3 dpa Mcm5 Mcm6 Mcm7 CycE E2f1 Myc CycD Ku80 Irbp Cht8 "
         "CG3344 CG4734 Fanci lok rad50 Xrp1 p53 Hsc70-3 Xbp1 crc PEK Ire1 Atf6 Pdi ERp60 Calr Hsp70Aa Hsp83 Hsp68 Hsp27 "
         "Hsp26 Hsc70-4 DnaJ-1 Hsf mei-41 tefu spn-A mus304 rpr hid Hr38 sr kay Jra puc CrebA").split()
genes = list(dict.fromkeys(named)) + [f"CG9{i:04d}" for i in range(8000)]
ids = [f"FBgn{i:07d}" for i in range(len(genes))]
with gzip.open(f"{out}/Drosophila_melanogaster.gtf.gz", "wt") as h:
    for i, (g, gid) in enumerate(zip(genes, ids)):
        h.write(f'2L\tFlyBase\tgene\t{1000*i+1}\t{1000*i+800}\t.\t+\t.\tgene_id "{gid}"; gene_name "{g}";\n')
samples = [(r, a, s, k) for a in (0, 2, 7, 21) for s in ("F", "M") for r in ("OL", "CB", "VNC") for k in (1, 2, 3)]
if layout in ("xlsx", "shuffled"):
    titles = [f"{a}d{'Female' if s == 'F' else 'Male'}{r}{k}" for r, a, s, k in samples]
else:
    titles = [f"{r}_{a}d_{s}_rep{k}" for r, a, s, k in samples]
gsm = [f"GSM46{i:05d}" for i in range(len(samples))]
base = rng.lognormal(2, 1.5, len(genes))
gi = {g: i for i, g in enumerate(genes)}
cols = {}
rng_region = {r: np.random.default_rng({"OL": 11, "CB": 12, "VNC": 13}[r]).normal(0, 0.4, len(genes)) for r in ("OL", "CB", "VNC")}
poly = {"OL": {0: .01, 2: .02, 7: .08, 21: .22}, "CB": {0: .01, 2: .01, 7: .02, 21: .04}, "VNC": {0: .01, 2: .01, 7: .02, 21: .03}}
for (r, a, s, k), t in zip(samples, titles):
    lf = rng.normal(0, 0.15, len(genes)) + rng_region[r]
    lf += -0.02 * np.log2(a + 1)                      # general ageing
    for g in "Ku80 Irbp Cht8 CG3344 CG4734 Fanci lok rad50 Xrp1".split():   # OL-specific DNA damage response
        lf[gi[g]] += 0.15 * np.log2(a + 1) * (r == "OL")
    if scenario == "effect":
        for g in "cnc Keap1 Gclc Gclm GstD1 GstD2 GstE1 GstS1 Zw Pgd Men".split():
            lf[gi[g]] += 3.0 * poly[r][a]
    mu = base * np.exp(lf); mu = mu / mu.sum() * 2e7
    cols[t] = rng.negative_binomial(20, 20 / (20 + mu))
counts = pd.DataFrame(cols); counts.insert(0, "gene_id", ids)
if layout in ("xlsx", "shuffled"):
    m = counts.drop(columns="gene_id").astype(float) * 1.03
    m.columns = [f"Sample_{72006 + i}" for i in range(m.shape[1])]
    if layout == "shuffled":
        m = m[list(np.random.default_rng(9).permutation(m.columns))]
        m.columns = [f"Sample_{72006 + i}" for i in range(m.shape[1])]
    m.insert(0, "id", genes)
    m.to_excel(f"{out}/GSE153165_polyploidy_depth_normalized_counts.xlsx", index=False, sheet_name="Sheet")
elif layout == "table":
    counts.to_csv(f"{out}/GSE153165_raw_counts.txt.gz", sep="\t", index=False)
else:
    with tarfile.open(f"{out}/GSE153165_RAW.tar", "w") as tar:
        for g_, t in zip(gsm, titles):
            data = pd.DataFrame({"gene": ids, "count": counts[t]}).to_csv(sep="\t", index=False).encode()
            info = tarfile.TarInfo(f"{g_}_{t}.counts.txt"); info.size = len(data); tar.addfile(info, io.BytesIO(data))
with gzip.open(f"{out}/GSE153165_series_matrix.txt.gz", "wt") as h:
    h.write('!Series_title\t"Polyploidy in the adult Drosophila melanogaster brain"\n')
    h.write("!Sample_title\t" + "\t".join(f'"{t}"' for t in titles) + "\n")
    h.write("!Sample_geo_accession\t" + "\t".join(f'"{g}"' for g in gsm) + "\n")
    tissue = {"OL": "optic lobe", "CB": "central brain", "VNC": "ventral nerve cord"}
    if layout in ("xlsx", "shuffled"):
        h.write("!Sample_source_name_ch1\t" + "\t".join('"Drosophila brain"' for _ in samples) + "\n")
        h.write("!Sample_characteristics_ch1\t" + "\t".join('"genotype: Canton S"' for _ in samples) + "\n")
        print("ok", layout, scenario); sys.exit()
    h.write("!Sample_characteristics_ch1\t" + "\t".join(f'"tissue: {tissue[r]}"' for r, a, s, k in samples) + "\n")
    h.write("!Sample_characteristics_ch1\t" + "\t".join(f'"age: {a} days"' for r, a, s, k in samples) + "\n")
    h.write("!Sample_characteristics_ch1\t" + "\t".join(f'"Sex: {"male" if s == "M" else "female"}"' for r, a, s, k in samples) + "\n")
print("ok", layout, scenario)
