# Self-contained demo: fold designed hairpins with Nussinov and with stacking rewards.
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Arc

from folding import fold, pairs_from_dotbracket, pair_f1

TARGET = "((((((((.....))))))))....((((((......))))))....(((((....)))))"


def design(db, rng):
    seq = ["A"] * len(db)
    for i, j in pairs_from_dotbracket(db):
        a, b = [("G", "C"), ("C", "G"), ("A", "U"), ("U", "A")][rng.choice(4, p=[0.35, 0.35, 0.15, 0.15])]
        seq[i], seq[j] = a, b
    for k, c in enumerate(db):
        if c == ".":
            seq[k] = rng.choice(["A", "A", "C", "U"])
    return "".join(seq)


def arcs(ax, pairs, n, colour, label, above=True):
    for i, j in pairs:
        w = j - i
        ax.add_patch(Arc(((i + j) / 2, 0), w, w, theta1=0 if above else 180, theta2=180 if above else 360,
                         color=colour, lw=1.2))
    ax.plot([], [], color=colour, label=label)


def main():
    os.makedirs("figures", exist_ok=True)
    os.makedirs("results", exist_ok=True)
    rng = np.random.default_rng(2)
    true_pairs = pairs_from_dotbracket(TARGET)
    rows = []
    for rep in range(12):
        seq = design(TARGET, rng)
        shuffled = "".join(rng.permutation(list(seq)))
        for name, stacking in [("Nussinov (max pairs)", False), ("stacking-aware", True)]:
            db, pairs, _ = fold(seq, stacking=stacking)
            _, sh_pairs, _ = fold(shuffled, stacking=stacking)
            rows.append({"design": rep, "method": name, "f1_vs_designed": round(pair_f1(pairs, true_pairs), 3),
                         "pairs_designed_seq": len(pairs), "pairs_shuffled_seq": len(sh_pairs)})
    df = pd.DataFrame(rows)
    table = df.groupby("method")[["f1_vs_designed", "pairs_designed_seq", "pairs_shuffled_seq"]].mean().round(2).reset_index()
    table.to_csv("results/summary.csv", index=False)

    seq = design(TARGET, np.random.default_rng(0))
    db_n, p_n, dp = fold(seq)
    db_s, p_s, _ = fold(seq, stacking=True)
    with open("results/structures.txt", "w") as fh:
        fh.write(f"sequence        {seq}\ndesigned        {TARGET}\nNussinov        {db_n}\nstacking-aware  {db_s}\n")
    fig, (ax, ax2) = plt.subplots(1, 2, figsize=(17, 5.2), gridspec_kw={"width_ratios": [2.2, 1]})
    n = len(seq)
    arcs(ax, true_pairs, n, "#444444", "designed structure", above=True)
    arcs(ax, p_n, n, "#c0504d", f"Nussinov (F1 {pair_f1(p_n, true_pairs):.2f})", above=False)
    arcs(ax, [], n, "#2a7f62", f"stacking-aware (F1 {pair_f1(p_s, true_pairs):.2f}, see results/)")
    span = max(j - i for i, j in true_pairs + p_n) / 2 + 1
    ax.set_xlim(-1, n)
    ax.set_ylim(-span, span)
    ax.axhline(0, color="k", lw=0.5)
    ax.set_yticks([])
    ax.set_xlabel("nucleotide position")
    ax.set_title("Arc diagram: designed structure (top) vs max-pairs prediction (bottom)")
    ax.legend(loc="upper right", fontsize=8)
    x = np.arange(len(table))
    ax2.bar(x - 0.2, table["pairs_designed_seq"], 0.4, label="designed sequence", color="#2a7f62")
    ax2.bar(x + 0.2, table["pairs_shuffled_seq"], 0.4, label="same letters, shuffled", color="#bbbbbb")
    for i, f in enumerate(table["f1_vs_designed"]):
        ax2.text(i, 1, f"F1 {f:.2f}", ha="center", fontsize=9)
    ax2.set_xticks(x, table["method"])
    ax2.set_ylabel("base pairs predicted (mean of 12 designs)")
    ax2.set_title("Shuffled RNA 'folds' too")
    ax2.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig("figures/demo.png", dpi=120)
    print(table.to_string(index=False))


if __name__ == "__main__":
    main()
