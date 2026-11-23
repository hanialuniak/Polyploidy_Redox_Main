#@title Cell 15 - Save and download the results
# The file names record the settings, so runs with different settings do not overwrite each other
tag = "_".join([{"All lines": "allTP53", "TP53-normal lines only": "TP53normal", "Adjust for TP53 status": "TP53adjusted"}[TP53_MODE],
                "prolif" if ADJUST_PROLIFERATION else "noprolif",
                "brain" if LINEAGE == "CNS/Brain only" else "alltissues"])
outputs = {f"route_a_programmes_{tag}.csv": results, f"route_a_findings_{tag}.csv": findings,
           f"route_a_genes_{tag}.csv": gene_results, f"route_a_geminin_{tag}.csv": geminin}
for name, frame in outputs.items():
    frame.to_csv(os.path.join(DATA_DIR, name), index=False)
print("Saved to your Drive folder:", ", ".join(outputs))
try:
    from google.colab import files
    for name in list(outputs):
        files.download(os.path.join(DATA_DIR, name))
except Exception:
    print("Open the Files panel on the left to download the tables and charts.")
