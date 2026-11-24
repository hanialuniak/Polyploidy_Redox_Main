#!/usr/bin/env python3
"""Rebuild every notebook in notebooks/ from its cell sources in src/."""
import pathlib
import shutil
import subprocess
import sys

ROOT = pathlib.Path(__file__).resolve().parent
NOTEBOOKS = ROOT / "notebooks"
failed = []
for build in sorted(ROOT.glob("src/*/build*.py")) + sorted(ROOT.glob("src/*/*/build*.py")):
    where = build.parent
    print(f"--- {build.relative_to(ROOT)}")
    r = subprocess.run([sys.executable, build.name], cwd=where, capture_output=True, text=True)
    print((r.stdout or r.stderr).strip()[:300])
    if r.returncode != 0:
        failed.append(str(build.relative_to(ROOT)))
        continue
    for nb in where.glob("*.ipynb"):
        shutil.copy(nb, NOTEBOOKS / nb.name)
        print(f"    -> notebooks/{nb.name}")
    for pdf in where.glob("*.pdf"):
        shutil.copy(pdf, ROOT / "docs" / pdf.name)
        print(f"    -> docs/{pdf.name}")
if failed:
    print("\nFAILED:", ", ".join(failed))
    sys.exit(1)
print("\nAll notebooks rebuilt.")
