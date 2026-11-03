"""
Like-charge regions (LCRs) of yeast FG-Nups (main text Methods 'Sequence analysis', Fig. 1b).

Definitions
  positive residues: K, R      negative residues: D, E      (sequences read from the N-terminus)
  A positive LCR starts at K or R and runs up to, but not including, the next D or E. It may hold
  further K/R and any number of uncharged residues. A negative LCR is the same with signs swapped.
  The last LCR of a segment ends at the segment end. LCR length = number of residues in it.
Segments
  The charge-rich segment is the C-terminal part (anchored to the scaffold), after the boundary
  found from the cumulative charge profile (ESI Fig. S1); the charge-poor segment is the N-terminal
  free end, up to the boundary. Boundaries (residue number from the N-terminus) are given below.
Input : yeast_FG_Nups.fasta (NCBI sequences; accession numbers in the headers)
Output: lcr_lengths.csv (every LCR) and lcr_histogram.csv (Fig. 1b, fractions of all pooled LCRs)
Usage : python3 lcr_analysis.py yeast_FG_Nups.fasta
"""
import sys, numpy as np, pandas as pd

BOUNDARY = {"Nup116": 774, "Nup100": 609, "Nup49": 275, "Nup57": 286, "Nsp1": 243}   # central-channel FG-Nups
BINS = [(1, 1), (2, 2), (3, 3), (4, 4), (5, 5), (6, 10), (11, 20), (21, 50), (51, 100), (101, 300), (301, 600)]

def read_fasta(f):
    seqs, name = {}, None
    for line in open(f):
        line = line.strip()
        if line.startswith(">"): name = line[1:].split()[0]; seqs[name] = ""
        elif name: seqs[name] += line
    return seqs

def lcrs(seg):
    """list of (sign, length) for every LCR in a segment"""
    out, sign, start = [], None, None
    for i, c in enumerate(seg):
        s = 1 if c in "KR" else (-1 if c in "DE" else 0)
        if s == 0: continue
        if sign is None: sign, start = s, i
        elif s != sign: out.append((sign, i - start)); sign, start = s, i
    if sign is not None: out.append((sign, len(seg) - start))          # last LCR ends at the segment end
    return out

seqs = read_fasta(sys.argv[1] if len(sys.argv) > 1 else "yeast_FG_Nups.fasta")
rows = []
for name, b in BOUNDARY.items():
    s = seqs[name]
    for seg, label in ((s[b:], "charge-rich"), (s[:b], "charge-poor")):
        for sign, L in lcrs(seg):
            rows.append((name, label, "positive" if sign > 0 else "negative", L))
d = pd.DataFrame(rows, columns=["FG_Nup", "segment", "sign", "length"])
d.to_csv("lcr_lengths.csv", index=False)

N = len(d)
lab = ["%d" % a if a == z else "%d-%d" % (a, z) for a, z in BINS]
d["bin"] = pd.cut(d.length, [a - 0.5 for a, _ in BINS] + [BINS[-1][1] + 0.5], labels=lab)
h = d.groupby(["segment", "sign", "bin"], observed=False).size().unstack("bin") / N   # all bars sum to 1
h.to_csv("lcr_histogram.csv")
rich = d[d.segment == "charge-rich"]; poor = d[d.segment == "charge-poor"]
print(f"N = {N} LCRs")
print(f"charge-rich: {100 * (rich.length <= 20).mean():.0f}% of LCRs are 20 residues or shorter")
print(f"charge-poor: {100 * (poor.length > 50).mean():.0f}% of LCRs are longer than 50 residues")
print(h.round(3).to_string())
