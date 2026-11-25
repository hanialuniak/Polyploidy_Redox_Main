"""Write mock EBI Single Cell Expression Atlas files in the same shape as E-GEOD-234602.

    python3 simulate.py /tmp/mock effect     # salivary-gland-specific ER rise planted
    python3 simulate.py /tmp/mock null       # no gland-specific effect
"""
import gzip, os, sys

import numpy as np

out, mode = sys.argv[1], (sys.argv[2] if len(sys.argv) > 2 else "effect")
acc = "E-GEOD-234602"
d = os.path.join(out, "sc_experiments", acc)
os.makedirs(d, exist_ok=True)
rng = np.random.default_rng(0)

REAL = ("Sec61beta Sec61alpha SrpRalpha Spase25 p24-1 zetaCOP Sec23 Sec31 Srp54 Sar1 fkh sage CrebA hth "
        "Hsc70-3 Pdi ERp60 Calr Xbp1 Ire1 PEK crc Ero1L Jafrac2 Jafrac1 Prx2540-1 Prx2540-2 PHGPx Trxr-1 "
        "Trx-2 Srx Gclc Gclm Gss1 Gr1 Grx1 Sod1 Sod2 Sod3 Cat Prx3 Prx5 Trxr-2 Duox Nox Su(var)3-3 dhd TrxT "
        "Rrp1 Ogg1 CycB polo aurB stg Cdk1 CycA pav tum GstD1 GstD2 GstE1 GstS1 Zw Pgd Men Idh Tkt Taldo "
        "cnc Keap1 ND-75 SdhA COX4 blw RpL3 RpL4 RpS2 RpS6 geminin dup Cdc6 Orc1 Mcm2 Mcm5 CycE E2f1 Myc "
        "p53 Ku80 lok tefu rad50 Xrp1 Hsp70Aa Hsp83 Hsp27 Hsf Hr38 sr kay puc Act5C RpL32 Gapdh1 Tbp Pgk").split()
# FlyBase ids, plus a symbol map written into a small GTF so the notebook can convert
sym2id = {g: f"FBgn{i:07d}" for i, g in enumerate(REAL, start=1)}
bg = {f"bg{i}": f"FBgn{900000 + i:07d}" for i in range(1200)}
sym2id.update(bg)
genes = list(sym2id)
ids = [sym2id[g] for g in genes]

GROUPS = [("salivary gland", "embryonic stage 10-12", 400),
          ("salivary gland", "embryonic stage 13-16", 450),
          ("tracheal cell", "embryonic stage 10-12", 1500),
          ("tracheal cell", "embryonic stage 13-16", 1600),
          ("central nervous system", "embryonic stage 10-12", 1800),
          ("central nervous system", "embryonic stage 13-16", 1700)]

# baseline per gene, then group effects
base = rng.uniform(1.5, 6.0, len(genes))
ER = "Hsc70-3 Pdi ERp60 Calr Xbp1 Ire1 PEK crc Ero1L Sec61alpha Sec61beta Srp54".split()
CREBA = "Sec61beta Sec61alpha SrpRalpha Spase25 p24-1 zetaCOP Sec23 Sec31 Srp54 Sar1".split()
MITOTIC = "CycB polo aurB stg Cdk1 CycA pav tum".split()
gi = {g: i for i, g in enumerate(genes)}

cells, meta = [], []
for ct, stage, n in GROUPS:
    eff = np.zeros(len(genes))
    late = "13-16" in stage
    sg = ct == "salivary gland"
    if late:
        eff += rng.normal(0.05, 0.1, len(genes))          # embryo-wide maturation
        for g in MITOTIC:
            eff[gi[g]] -= 0.8                              # everything slows mitosis late
    if sg and late and mode == "effect":
        for g in CREBA:
            eff[gi[g]] += 1.6                              # CrebA secretory ramp, gland only
        for g in ER:
            eff[gi[g]] += 1.3
        eff[gi["Ero1L"]] += 1.3
        eff[gi["Jafrac2"]] += 0.05                         # ER peroxiredoxin does NOT follow
        eff[gi["Trxr-1"]] += 0.05
        eff[gi["Trx-2"]] += 0.05
        eff[gi["Cat"]] += 1.2
    if sg:
        for g in ("fkh", "sage", "CrebA"):
            eff[gi[g]] += 2.5
    mu = 2 ** (base + eff)
    for _ in range(n):
        cells.append(rng.poisson(mu / mu.sum() * rng.integers(600, 2500)))
        meta.append((ct, stage))

M = np.array(cells).T                                      # genes x cells
names = [f"CELL{i:06d}" for i in range(M.shape[1])]
with open(os.path.join(d, f"{acc}.aggregated_filtered_counts.mtx"), "w") as h:
    nz = np.flatnonzero(M)
    h.write("%%MatrixMarket matrix coordinate integer general\n")
    h.write(f"{M.shape[0]} {M.shape[1]} {len(nz)}\n")
    for k in nz:
        i, j = divmod(int(k), M.shape[1])
        h.write(f"{i + 1} {j + 1} {int(M[i, j])}\n")
open(os.path.join(d, f"{acc}.aggregated_filtered_counts.mtx_rows"), "w").write(
    "".join(f"{i}\t{i}\n" for i in ids))
# EBI also ships a decorated rows file: id in column 1, gene symbol in column 2
# deliberately shuffled: the decorated file is a lookup, not a parallel array. Reading it positionally
# renames every gene after its neighbour, which is exactly the bug this guards against.
_dec = sorted(genes, key=lambda g: sym2id[g])[::-1]
open(os.path.join(d, f"{acc}.aggregated_filtered_counts.decorated.mtx_rows"), "w").write(
    "".join(f"{sym2id[g]}\t{g}\n" for g in _dec))
open(os.path.join(d, f"{acc}.aggregated_filtered_counts.mtx_cols"), "w").write(
    "".join(f"{c}\n" for c in names))
with open(os.path.join(d, f"ExpDesign-{acc}.tsv"), "w") as h:
    # the real file uses normalised names with underscores, and carries ontology columns alongside
    h.write("id\tdevelopmental_stage\tinferred_cell_type_-_authors_labels\t"
            "inferred_cell_type_-_authors_labels_ontology\tdevelopmental_stage_ontology\n")
    for c, (ct, stage) in zip(names, meta):
        h.write(f"{c}\t{stage}\t{ct}\tFBbt:0000000\tFBdv:0000000\n")

gtf = os.path.join(out, "Drosophila_melanogaster.gtf.gz")
with gzip.open(gtf, "wt") as h:
    for g, i in sym2id.items():
        h.write(f'2L\tFlyBase\tgene\t1\t2\t.\t+\t.\tgene_id "{i}"; gene_name "{g}";\n')
print(f"mock ({mode}) written to {d}; gtf at {gtf}")
