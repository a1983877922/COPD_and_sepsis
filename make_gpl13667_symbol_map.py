# -*- coding: utf-8 -*-
"""
Build an offline probe -> gene symbol mapping table for GPL13667 (Affymetrix HG-U219,
used by GSE65682).
Data sources (both official Bioconductor packages, local builds, no internet needed):
  hgu219.db    inst/extdata/hgu219.sqlite   : probe_id -> entrez gene_id
  org.Hs.eg.db inst/extdata/org.Hs.eg.sqlite: entrez _id  -> gene symbol
Output: GPL13667_symbol_map.csv (two columns: probe, symbol; tab-separated, one row per probe)
Multi-gene probes have their symbols joined with " /// " (matching the style of the Gene
Symbol column in the NCBI SOFT file; the R-side load_cohort fallback keeps the first symbol).
"""
import sqlite3, csv, os

HERE   = os.path.dirname(os.path.abspath(__file__))
HGU219 = os.path.join(HERE, "hgu219_tmp", "hgu219.db", "inst", "extdata", "hgu219.sqlite")
ORGEG  = os.path.join(HERE, "org.Hs.eg.db", "inst", "extdata", "org.Hs.eg.sqlite")
OUT    = os.path.join(HERE, "GPL13667_symbol_map.csv")

assert os.path.exists(HGU219), f"missing hgu219.sqlite: {HGU219}"
assert os.path.exists(ORGEG),  f"missing org.Hs.eg.sqlite: {ORGEG}"

# 1) probe -> [entrez ids]
con = sqlite3.connect(HGU219)
rows = con.execute("SELECT probe_id, gene_id FROM probes").fetchall()
con.close()
probe2gene = {}
for probe, gid in rows:
    probe2gene.setdefault(probe, []).append(str(gid))

# 2) entrez -> symbol
# Note: in org.Hs.eg.db 3.23, gene_info._id is an internal ID, not the NCBI Entrez ID;
#     Entrez lives in the genes table (gene_id column), so join gene_info (internal _id -> symbol).
con = sqlite3.connect(ORGEG)
cur = con.cursor()
gi = cur.execute(
    "SELECT g.gene_id, i.symbol FROM genes g "
    "JOIN gene_info i ON g._id = i._id WHERE i.symbol IS NOT NULL"
).fetchall()
con.close()
eg2sym = {str(gid): sym for gid, sym in gi if sym}

# 3) probe -> symbol
out = [("probe", "symbol")]
n_unmapped = 0
for probe in probe2gene:
    syms = [eg2sym[g] for g in probe2gene[probe] if g in eg2sym]
    syms = list(dict.fromkeys(syms))          # dedupe, keep order
    if not syms:
        n_unmapped += 1
        continue
    out.append((probe, " /// ".join(syms)))

with open(OUT, "w", newline="", encoding="utf-8") as f:
    w = csv.writer(f, delimiter="\t", lineterminator="\n")
    w.writerows(out)

print(f"unique probes: {len(probe2gene)} | mapped to symbol: {len(out)-1} | no symbol (control/unmapped): {n_unmapped}")
print("output:", OUT)

# 4) spot-check that key genes have probes on the array
for g in ["SOCS3", "STAT1", "STAT3", "GBP1", "IFI44L", "JAK1", "MX1"]:
    hits = [p for p, s in out[1:] if s.split(" /// ")[0] == g]
    print(f"  {g}: {len(hits)} probes  e.g. {hits[:3]}")
