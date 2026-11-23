#@title Cell 5 - Gene programmes (fly gene symbols)
PROGRAMMES = {
    "REDOX: cnc (NRF2) pathway": "cnc Keap1 Gclc Gclm GstD1 GstD2 GstE1 GstS1 Zw Pgd Men",
    "REDOX: thioredoxin/peroxiredoxin": "Trxr-1 Trxr-2 Trx-2 Jafrac1 Jafrac2 Prx3 Prx5 PHGPx",
    "REDOX: antioxidant enzymes": "Sod1 Sod2 Cat Eip71CD SelR Ogg1",
    "REDOX: iron": "Fer1HCH Fer2LCH Irp1A Tsf1 Mvl",
    "REDOX: mitochondrial respiration": "ND-75 ND-49 ND-30 ND-42 SdhA SdhB UQCR-C1 UQCR-C2 Cyt-c-p COX4 COX5A ATPsynbeta blw",
    "EXPLORE: licensing and cell cycle re-entry": "geminin dup Cdc6 Orc1 Mcm2 Mcm3 dpa Mcm5 Mcm6 Mcm7 CycE E2f1 Myc CycD",
    "EXPLORE: p53 head DNA damage response": "Ku80 Irbp Cht8 CG3344 CG4734 Fanci lok rad50 Xrp1 p53",
    "CONTROL: unfolded protein response": "Hsc70-3 Xbp1 crc PEK Ire1 Atf6 Pdi ERp60 Calr",
    "CONTROL: heat shock": "Hsp70Aa Hsp83 Hsp68 Hsp27 Hsp26 Hsc70-4 DnaJ-1 Hsf",
    "CONTROL: DNA damage": "p53 mei-41 tefu lok spn-A mus304 Xrp1 rpr hid",
    "CONTROL: immediate-early genes": "Hr38 sr kay Jra puc CrebA",
}
# Some genes are listed under a different name in the annotation, so each one is looked up under its synonyms too.
# FlyBase calls Geminin "geminin" (CG3183); "Gem" is only a synonym, which is why earlier runs could not find it.
ALIASES = {"geminin": ["Gem", "CG3183"], "dup": ["Cdt1"], "dpa": ["Mcm4"], "Zw": ["G6PD", "G6pd"],
           "Ero1L": ["Ero1", "CG1333"], "PHGPx": ["CG12013"], "Gss1": ["Gss"], "Gr1": ["Gr", "GR"],
           "Trxr-1": ["Trxr1"], "Jafrac1": ["Prx2"], "Men-b": ["Men-b", "CG5889"]}


def resolve(g, present):
    if g in present:
        return g
    for a in ALIASES.get(g, []):
        if a in present:
            return a
    return None


_present = set(expressed.index)
_renamed, _resolved = {}, {}
for n, v in PROGRAMMES.items():
    out = []
    for g in v.split():
        r = resolve(g, _present)
        if r:
            out.append(r)
            if r != g:
                _renamed[g] = r
    _resolved[n] = out
PROGRAMMES = _resolved
if _renamed:
    print("Found under a synonym:", ", ".join(f"{k} -> {v}" for k, v in _renamed.items()))
for n, v in PROGRAMMES.items():
    print(f"{n}: {len(v)} expressed genes")
PROGRAMMES = {n: v for n, v in PROGRAMMES.items() if len(v) >= 4}
