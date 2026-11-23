import gzip, io, os, sys, tarfile, numpy as np
out = sys.argv[1]; rng = np.random.default_rng(5)
named = ("GstD1 GstD2 GstE1 GstE3 GstS1 GstO1 GstT1 GstZ1 Gclc Gclm Gss1 Gr1 PHGPx Trxr-1 Jafrac1 Prx3 Prx5 Sod1 Sod2 Cat Zw Pgd Men Idh "
         "cnc Keap1 Cyp6a2 Ero1L Duox Nox ND-75 ND-49 SdhA COX4 blw Xbp1 Hsc70-3 Ire1 PEK crc Pdi ERp60 Calr RpL3 RpL4 RpS2 RpS3 RpS6 "
         "p53 Ku80 lok Xrp1 rad50 Hsp70Aa Hsp83 Hsp68 Hsp27 Hsf Hr38 sr kay Jra puc Act5C alphaTub84B RpL32 Gapdh1 Tbp eIF4A "
         "CycB polo aurB stg Cdk1 CycA Det pav PCNA RnrS Mcm2 Mcm5 Orc1 dup").split()
genes = named + [f"CG{i}" for i in range(5000)]
fb = [f"FBgn{1000000 + i}" for i in range(len(genes))]
os.makedirs(out, exist_ok=True)
with gzip.open(f"{out}/dm.gtf.gz", "wt") as h:
    for f, g in zip(fb, genes):
        h.write(f'2L\tFlyBase\tgene\t1\t100\t.\t+\t.\tgene_id "{f}"; gene_name "{g}";\n')
base = rng.gamma(1.0, 300, len(genes)) + 5
mit = {"CycB", "polo", "aurB", "stg", "Cdk1", "CycA", "Det", "pav"}
gst = {g for g in genes if g.startswith("Gst")}
def make(gse, groups, test, layout):
    titles, gsms, cols = [], [], []
    for gi, (grp, n) in enumerate(groups):
        tissue = rng.normal(0, 0.5, len(genes))
        for k in range(n):
            mu = base * np.exp(tissue)
            if grp in test:
                mu = mu * np.where([g in mit for g in genes], 0.15, 1) * np.where([g in gst for g in genes], 2.0, 1)
            cols.append(rng.negative_binomial(20, 20 / (20 + mu * rng.uniform(0.7, 1.4))))
            titles.append(f"{grp} replicate {k + 1}" if layout == "xlsx" else f"{grp}_rep{k + 1}"); gsms.append(f"GSM34{len(gsms) + 10:04d}")
    d = f"{out}/geo/series/GSE121nnn/{gse}"; os.makedirs(d + "/matrix", exist_ok=True); os.makedirs(d + "/suppl", exist_ok=True)
    with gzip.open(f"{d}/matrix/{gse}_series_matrix.txt.gz", "wt") as h:
        h.write("!Sample_title\t" + "\t".join(f'"{t}"' for t in titles) + "\n")
        h.write("!Sample_geo_accession\t" + "\t".join(f'"{g}"' for g in gsms) + "\n")
        h.write("!Sample_source_name_ch1\t" + "\t".join(f'"{t.split("_")[0]}"' for t in titles) + "\n")
        h.write("!series_matrix_table_begin\n\"ID_REF\"\t" + "\t".join(gsms) + "\n!series_matrix_table_end\n")
    if layout == "xlsx":
        import pandas as pd
        codes = ["SG" if t.startswith("salivary") else "BD" for t in titles]
        df = pd.DataFrame({"FlyBase ID": fb, "Gene": genes, "Length": rng.integers(500, 5000, len(fb))})
        for k, (c, col) in enumerate(zip(codes, cols)): df[f"{c}{k % 3 + 1}"] = col
        b = io.BytesIO(); df.to_excel(b, index=False); open(f"{d}/suppl/{gse}_GSF1939_SG_RNAseq_raw_counts.xlsx.gz", "wb").write(gzip.compress(b.getvalue()))
        de = pd.DataFrame({"gene": fb, "baseMean": 5.0, "log2FoldChange": 0.1, "pvalue": 0.5, "padj": 0.9})
        b = io.BytesIO(); de.to_excel(b, index=False); open(f"{d}/suppl/{gse}_SG_vs_BD.xlsx.gz", "wb").write(gzip.compress(b.getvalue()))
    elif layout == "table":
        with gzip.open(f"{d}/suppl/{gse}_counts.txt.gz", "wt") as h:
            h.write("gene_id\t" + "\t".join(titles) + "\n")
            for i, f in enumerate(fb): h.write(f + "\t" + "\t".join(str(c[i]) for c in cols) + "\n")
    else:
        with tarfile.open(f"{d}/suppl/{gse}_RAW.tar", "w") as tf:
            for g, t, c in zip(gsms, titles, cols):
                body = gzip.compress(("".join(f"{f}\t{v}\n" for f, v in zip(fb, c)) + "__no_feature\t100\n").encode())
                ti = tarfile.TarInfo(f"{g}_{t}.counts.txt.gz"); ti.size = len(body); tf.addfile(ti, io.BytesIO(body))
make("GSE121956", [("salivary gland", 3), ("brain and disc", 3)], {"salivary gland"}, "xlsx")
make("GSE121955", [("CycA dsRNA", 3), ("Myb dsRNA", 3), ("GFP dsRNA", 3)], {"CycA dsRNA", "Myb dsRNA"}, "tar")
print("ok")
