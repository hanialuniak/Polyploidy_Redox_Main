#@title Cell 12 - Is the germline signal real, or follicle mRNA leaking into germline droplets?
# Cell 11 called 700-odd genes "shared with the germline". The list cannot be taken at face value, and the
# reason is visible in it: the top of the shared list is Yp1, Yp2, Yp3, Vm34Ca, Vm26Aa, Vm26Ab. Vitelline
# membrane genes are made by FOLLICLE cells, not by nurse cells. They have no business changing in the
# germline at all, so their presence says the germline column is partly picking up follicle transcripts.
#
# How that happens: free mRNA from lysed cells is captured along with every cell, so a transcript that a
# neighbouring follicle cell is producing in enormous quantity appears at a low level in germline droplets
# too. Nurse cells sit inside egg chambers whose follicle cells are pouring out eggshell and yolk message,
# while germarial cysts do not - so ambient mRNA alone produces a germline "change" in the same direction as
# the follicle one. That is exactly the artefact that would manufacture shared calls.
#
# The check: a real germline change should be in a gene the germline actually expresses. Two requirements,
# and the Vm genes are the built-in test of whether the check works - if they are not flagged, it is not.
MIN_GL_DETECT = 0.25       #@param {type:"number"}
MAX_LEVEL_GAP = 2.0        #@param {type:"number"}

_gm = np.flatnonzero(((cells.gl_role == GMIT) & low_g).to_numpy())
_ge = np.flatnonzero(((cells.gl_role == GENDO) & low_g).to_numpy())
gl_order = [g for g in Lg.index if g in gpos]
gl_cols = [gpos[g] for g in gl_order]
gdet_m = pd.Series(np.asarray((X[_gm][:, gl_cols] > 0).mean(axis=0)).ravel(), index=gl_order)
gdet_e = pd.Series(np.asarray((X[_ge][:, gl_cols] > 0).mean(axis=0)).ravel(), index=gl_order)

gl_level = Lg.mean(axis=1)
fol_level = Lp.mean(axis=1)
gap = (fol_level - gl_level).reindex(gl_order)   # positive = much higher in follicle cells than germline

chk = pd.DataFrame({"gl_detect_germarium": gdet_m, "gl_detect_nurse": gdet_e,
                    "gl_level": gl_level.reindex(gl_order), "follicle_level": fol_level.reindex(gl_order),
                    "level_gap": gap})
chk["gl_detect_max"] = chk[["gl_detect_germarium", "gl_detect_nurse"]].max(axis=1)
chk["bleed_through_risk"] = np.where(
    chk.gl_detect_max < MIN_GL_DETECT, "germline barely expresses it",
    np.where(chk.level_gap > MAX_LEVEL_GAP, f"{2 ** MAX_LEVEL_GAP:.0f}x lower in germline than follicle", "ok"))

de2 = de.join(chk[["gl_detect_max", "level_gap", "bleed_through_risk"]])
de2["shared_confirmed"] = (de2.lineage.astype(str).str.startswith("BOTH")) & (de2.bleed_through_risk == "ok")

# Does the check catch the genes we already know are follicle-made? If not, stop reading here.
KNOWN_FOLLICLE = [g for g in "Vm34Ca Vm26Aa Vm26Ab Vm32E Cp36 Cp38 Cp15 Cp16 dec-1 psd".split() if g in de2.index]
if KNOWN_FOLLICLE:
    print("Validation of the check itself. These are eggshell genes, transcribed by follicle cells and not by "
          "nurse cells, so every one of them should be flagged:")
    display(de2.loc[KNOWN_FOLLICLE, ["log2FC", "germline_log2FC", "gl_detect_max", "level_gap",
                                     "bleed_through_risk"]].round(2))
    caught = (de2.loc[KNOWN_FOLLICLE, "bleed_through_risk"] != "ok").mean()
    print(f"Flagged: {caught:.0%} of them.")
    if caught < 0.6:
        print("*** The check is NOT catching genes known to be follicle-made. Do not use the verdicts below; "
              "lower MIN_GL_DETECT or raise MAX_LEVEL_GAP and look at the numbers for these genes by hand. ***")

shared0 = int(de2.lineage.astype(str).str.startswith("BOTH").sum())
shared1 = int(de2.shared_confirmed.sum())
print(f"\nShared with the germline before this check: {shared0:,} genes. After it: {shared1:,}.")
print(f"Removed: {shared0 - shared1:,} genes the germline barely expresses, or expresses far below the "
      f"follicle level, where ambient mRNA is the simpler explanation.")

keep2 = de2[de2.shared_confirmed & (de2.libraries_agree.fillna(True))]
keep2 = keep2.reindex(keep2.log2FC.abs().sort_values(ascending=False).index)
cols2 = [c for c in ["log2FC", "germline_log2FC", "detect_endo", "gl_detect_max", "level_gap", "expr_rank_pct",
                     "in_set"] if c in keep2.columns]
print(f"\n=== Shared with the germline AND expressed in the germline: {len(keep2):,} genes ===")
print("This is the shortlist for changes that track the endocycle rather than secretory differentiation.")
display(keep2.head(50)[cols2].round(2))

# Where the curated redox genes land once the check is applied.
WATCH = [g for g in "Jafrac1 Jafrac2 Ero1L Gclc Trx-2 Trxr-1 Rrp1 Ogg1 mre11 Cat Hsc70-3 Sec61beta geminin".split()
         if g in de2.index]
print("\nThe genes this project has been following, with the same check applied:")
display(de2.loc[WATCH, [c for c in ["log2FC", "germline_log2FC", "gl_detect_max", "level_gap",
                                    "bleed_through_risk", "lineage"] if c in de2.columns]].round(2))
print("""
Read the direction of the bias before reading anything else. Ambient follicle mRNA pushes the germline column
TOWARDS the follicle one, so it manufactures agreement and cannot manufacture disagreement. A redox gene that
still reads flat in the germline is therefore a stronger negative after this check than before it, not a weaker
one. What the check can overturn is the shared list, not the follicle-only calls.
""")
de2.sort_values("log2FC").to_csv(os.path.join(GDIR, "route_e_unbiased_de_checked.csv"))
print("Saved to", os.path.join(GDIR, "route_e_unbiased_de_checked.csv"))
