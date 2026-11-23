# Putting this on GitHub

## 1. Create the repository

On [github.com/new](https://github.com/new): name it `polyploidy-redox`, add no README,
no .gitignore and no licence — they are already here.

## 2. Upload

**In the browser**, which needs no tools installed: open the empty repository, click
**uploading an existing file**, and drag the whole unzipped folder in. GitHub keeps the folder
structure. Commit.

**Or from a terminal:**

```bash
cd polyploidy-redox
git init
git add .
git commit -m "Analysis of polyploidy and redox state across public datasets"
git branch -M main
git remote add origin https://github.com/YOUR-USERNAME/polyploidy-redox.git
git push -u origin main
```

## 3. Two small edits

- `CITATION.cff` — replace `USERNAME` in the repository URL.
- `README.md` — nothing needed, but check the data-source table still matches what you used.

## 4. Decide: public or private

Private is sensible until the work is published. Supervisors and collaborators can be added
under **Settings → Collaborators**. You can switch it to public later without losing history.

## Notes

- **No data is committed.** `.gitignore` excludes count matrices, GEO archives, annotation
  files and anything large. The notebooks fetch or point to their own data, so the repository
  stays a few megabytes.
- **Notebooks are committed with their outputs stripped** where possible. If you commit a
  notebook you have run in Colab, it will carry its output cells — useful as a record, but it
  makes diffs noisy. To clear them: **Edit → Clear all outputs** in Colab before downloading.
- **Opening a notebook from GitHub in Colab:** replace `github.com` with
  `colab.research.google.com/github` in the URL, so anyone can run it without downloading.
- If you later want a DOI for a thesis or paper, connect the repository to
  [Zenodo](https://zenodo.org) and make a release. `CITATION.cff` is already there for that.
