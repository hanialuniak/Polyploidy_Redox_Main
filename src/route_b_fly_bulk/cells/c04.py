#@title Cell 4 - Normalise and choose expressed genes
X = expr[info.column].astype(float)
integer = np.allclose(X.to_numpy(), np.round(X.to_numpy()))
if integer:
    X = np.log2(X / X.sum() * 1e6 + 1)
    print("Counts found: converted to log2 counts per million")
else:
    X = np.log2(X + 1) if X.max().max() > 50 else X
    print("Normalised values found (FPKM/TPM or similar): log2-transformed")
X.columns = range(X.shape[1])
expressed = X[(X > 1).mean(axis=1) >= 0.5]
print(f"Expressed genes: {len(expressed):,} of {len(X):,}")
