import gzip, os, sys, numpy as np
out, effect = sys.argv[1], sys.argv[2] == "effect"
rng = np.random.default_rng(3)
cyc = "Mki67 Top2a Ccnb1 Cdk1 Ccna2 Bub1 Aurkb Plk1 Cdc20 Ccnb2".split()
pc = "Glul Cyp2e1 Cyp1a2 Oat Lect2 Axin2 Rhbg Cyp7a1".split(); pp = "Cps1 Ass1 Arg1 Hal Sds Pck1 Gls2 Hsd17b13".split()
gst = "Gsta1 Gsta3 Gsta4 Gstm1 Gstm2 Gstm3 Gstp1 Gsto1 Gstt1 Gstz1 Mgst1".split()
other = "Gclc Gclm Gss Gsr Gpx1 Gpx4 Txn1 Txnrd1 Prdx1 Prdx2 Prdx3 Prdx5 Prdx6 Sod1 Sod2 Cat G6pd Pgd Me1 Idh1 Nqo1 Hmox1 Maoa Xdh Acox1".split()
genes = cyc + pc + pp + gst + other + [f"Gene{i}" for i in range(6000)]
ploidy = [2, 4, 8] * 5
animal = [i // 3 + 1 for i in range(15)]
titles = [f"{p}c_{a}" for p, a in zip(ploidy, animal)]
gsm = [f"GSM38629{38 + i}" for i in range(15)]
base = rng.normal(8, 2, len(genes)); lp = np.log2(ploidy)
M = base[:, None] + rng.normal(0, 0.25, (len(genes), 15)) + rng.normal(0, 0.4, (len(genes), 5))[:, np.array(animal) - 1] + np.array([0, .3, .6] * 5)[None, :]
zone = np.array([0.6, 0, -0.3] * 5) + rng.normal(0, 0.2, 15)
for i, g in enumerate(genes):
    if g in cyc: M[i] += (np.array(ploidy) == 8) * 1.0
    if g in pc: M[i] += zone
    if g in pp: M[i] -= zone
    if effect and g in gst: M[i] -= 0.35 * (lp - 1)
probes = [f"A_44_P{100000 + i}" for i in range(len(genes))]
mdir = f"{out}/series/GSE132nnn/GSE132409/matrix"; os.makedirs(mdir, exist_ok=True)
with gzip.open(f"{mdir}/GSE132409_series_matrix.txt.gz", "wt") as h:
    h.write('!Series_title\t"Transcriptome of adult rat hepatocytes with different ploidy statuses"\n')
    h.write('!Series_platform_id\t"GPL25947"\n')
    h.write("!Sample_title\t" + "\t".join(f'"{t}"' for t in titles) + "\n")
    h.write("!Sample_geo_accession\t" + "\t".join(f'"{g}"' for g in gsm) + "\n")
    h.write('!Sample_source_name_ch1\t' + "\t".join('"liver"' for _ in gsm) + "\n")
    h.write("!series_matrix_table_begin\n")
    h.write('"ID_REF"\t' + "\t".join(f'"{g}"' for g in gsm) + "\n")
    for p, row in zip(probes, M):
        h.write(f'"{p}"\t' + "\t".join(f"{2 ** v:.2f}" for v in row) + "\n")
    h.write("!series_matrix_table_end\n")
adir = f"{out}/platforms/GPL25nnn/GPL25947/annot"; os.makedirs(adir, exist_ok=True)
with gzip.open(f"{adir}/GPL25947.annot.gz", "wt") as h:
    h.write("^Annotation\n!Annotation_platform = GPL25947\n!platform_table_begin\n")
    h.write("ID\tGene title\tGene symbol\tGene ID\n")
    for p, g in zip(probes, genes):
        h.write(f"{p}\tx\t{g}\t1\n")
    h.write("!platform_table_end\n")
print("ok")
