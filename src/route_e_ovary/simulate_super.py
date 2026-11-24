"""Mock GEO laid out like a SuperSeries: nothing at series level, files under each sample."""
import gzip, io, os, shutil, sys
src, out = sys.argv[1], sys.argv[2]        # src = a folder already built by simulate.py
old = os.path.join(src, "geo", "series", "GSE146nnn", "GSE146040", "suppl")
gsms = sorted({n.split("_")[0] for n in os.listdir(old) if n.startswith("GSM")})
sdir = os.path.join(out, "geo", "series", "GSE136nnn", "GSE136162")
os.makedirs(os.path.join(sdir, "suppl"), exist_ok=True)
os.makedirs(os.path.join(sdir, "matrix"), exist_ok=True)
with gzip.open(os.path.join(sdir, "matrix", "GSE136162_series_matrix.txt.gz"), "wt") as h:
    h.write('!Series_title\t"A single-cell atlas of the adult Drosophila ovary"\n')
    h.write("!Sample_geo_accession\t" + "\t".join(f'"{g}"' for g in gsms) + "\n")
    h.write("!series_matrix_table_begin\n!series_matrix_table_end\n")
for g in gsms:
    d = os.path.join(out, "geo", "samples", g[:-3] + "nnn", g, "suppl")
    os.makedirs(d, exist_ok=True)
    for n in sorted(os.listdir(old)):
        if n.startswith(g):
            shutil.copy(os.path.join(old, n), os.path.join(d, n))
print("built SuperSeries mock for", gsms)
