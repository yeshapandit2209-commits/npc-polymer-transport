"""
Cumulative charge along yeast FG-Nup sequences: locates the charge-rich segment.

For each FG-Nup, counts charged residues (K, R, D, E), positive (K, R) and negative (D, E) residues
cumulatively from the N-terminus (free end). The slope of the cumulative count is the local charge
density. The boundary between the charge-poor free end and the charge-rich segment at the anchored
C-terminus is the point where the slope rises sharply. The boundaries used here are listed in BOUNDARY;
Nup159 and Nup60 have no boundary because their charges are spread along the whole chain.

Input : yeast_FG_Nups.fasta
Output: cumulative_charge.csv   residue, count of charged / positive / negative residues, per FG-Nup
        charge_segments.csv     length and charge fraction of each segment, slope ratio at the boundary
        cumulative_charge.png
Usage : python3 cumulative_charge.py yeast_FG_Nups.fasta
"""
import sys, numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt

BOUNDARY = {"Nup42": 365, "Nup116": 774, "Nup100": 609, "Nup49": 275, "Nup57": 286, "Nsp1": 243}
ORDER = ["Nup42", "Nup116", "Nup100", "Nup49", "Nup57", "Nsp1", "Nup159", "Nup60"]

def read_fasta(f):
    seqs, name = {}, None
    for line in open(f):
        line = line.strip()
        if line.startswith(">"): name = line[1:].split()[0]; seqs[name] = ""
        elif name: seqs[name] += line
    return seqs

seqs = read_fasta(sys.argv[1] if len(sys.argv) > 1 else "yeast_FG_Nups.fasta")
rows, seg = [], []
fig, axs = plt.subplots(2, 4, figsize=(13, 6))
for ax, n in zip(axs.ravel(), ORDER):
    s = seqs[n]; i = np.arange(1, len(s) + 1)
    pos = np.cumsum([c in "KR" for c in s]); neg = np.cumsum([c in "DE" for c in s]); tot = pos + neg
    rows += [(n, k, t, p, q) for k, t, p, q in zip(i, tot, pos, neg)]
    ax.plot(i, tot, "k", label="charged (K, R, D, E)"); ax.plot(i, pos, "C0", label="positive (K, R)"); ax.plot(i, neg, "C3", label="negative (D, E)")
    if n in BOUNDARY:
        b = BOUNDARY[n]; ax.axvline(b, color="0.4", ls="--")
        f_poor = tot[b - 1] / b; f_rich = (tot[-1] - tot[b - 1]) / (len(s) - b)
        seg.append((n, len(s), b, len(s) - b, (len(s) - b) / len(s), f_poor, f_rich, f_rich / f_poor))
        ax.set_title("%s  (%.0f%% | %.0f%% charged)" % (n, 100 * f_poor, 100 * f_rich))
    else:
        seg.append((n, len(s), 0, len(s), 1.0, np.nan, tot[-1] / len(s), np.nan))
        ax.set_title("%s  (%.0f%% charged)" % (n, 100 * tot[-1] / len(s)))
    ax.set_xlabel("Residue (from N-terminus)"); ax.set_ylabel("Cumulative count")
axs[0, 0].legend(fontsize=7, frameon=False); fig.tight_layout(); fig.savefig("cumulative_charge.png", dpi=200)
pd.DataFrame(rows, columns=["FG_Nup", "residue", "n_charged", "n_positive", "n_negative"]).to_csv("cumulative_charge.csv", index=False)
S = pd.DataFrame(seg, columns=["FG_Nup", "length", "boundary", "charge_rich_length", "charge_rich_fraction",
                               "charged_fraction_poor", "charged_fraction_rich", "slope_ratio"])
S.to_csv("charge_segments.csv", index=False); print(S.round(3).to_string(index=False))
