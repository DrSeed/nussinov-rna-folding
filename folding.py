# RNA secondary structure by dynamic programming: Nussinov and a stacking-aware variant.
import numpy as np

PAIRS = {("G", "C"), ("C", "G"), ("A", "U"), ("U", "A"), ("G", "U"), ("U", "G")}


def can_pair(a, b):
    return (a, b) in PAIRS


def _fold(seq, min_loop, score):
    n = len(seq)
    dp = np.zeros((n, n))
    for span in range(min_loop + 1, n):
        for i in range(n - span):
            j = i + span
            best = max(dp[i + 1, j], dp[i, j - 1])
            if can_pair(seq[i], seq[j]):
                best = max(best, dp[i + 1, j - 1] + score(seq, i, j))
            for k in range(i + 1, j):
                best = max(best, dp[i, k] + dp[k + 1, j])
            dp[i, j] = best
    return dp


def _traceback(seq, dp, min_loop, score):
    pairs, stack = [], [(0, len(seq) - 1)]
    while stack:
        i, j = stack.pop()
        if j - i <= min_loop:
            continue
        if np.isclose(dp[i, j], dp[i + 1, j]):
            stack.append((i + 1, j))
        elif np.isclose(dp[i, j], dp[i, j - 1]):
            stack.append((i, j - 1))
        elif can_pair(seq[i], seq[j]) and np.isclose(dp[i, j], dp[i + 1, j - 1] + score(seq, i, j)):
            pairs.append((i, j))
            stack.append((i + 1, j - 1))
        else:
            for k in range(i + 1, j):
                if np.isclose(dp[i, j], dp[i, k] + dp[k + 1, j]):
                    stack += [(i, k), (k + 1, j)]
                    break
    return sorted(pairs)


def _one(seq, i, j):
    return 1.0


def _stacking(seq, i, j):
    # Reward pairs that stack on a neighbouring pair, and G-C over A-U over G-U.
    strength = {"GC": 3.0, "CG": 3.0, "AU": 2.0, "UA": 2.0, "GU": 1.0, "UG": 1.0}[seq[i] + seq[j]]
    stacked = i + 1 < j - 1 and can_pair(seq[i + 1], seq[j - 1])
    return strength if stacked else strength - 2.5


def fold(seq, min_loop=3, stacking=False):
    """Return (dot-bracket structure, list of (i, j) pairs, DP matrix)."""
    seq = seq.upper().replace("T", "U")
    score = _stacking if stacking else _one
    dp = _fold(seq, min_loop, score)
    pairs = _traceback(seq, dp, min_loop, score)
    db = ["."] * len(seq)
    for i, j in pairs:
        db[i], db[j] = "(", ")"
    return "".join(db), pairs, dp


def pairs_from_dotbracket(db):
    stack, pairs = [], []
    for k, c in enumerate(db):
        if c == "(":
            stack.append(k)
        elif c == ")":
            pairs.append((stack.pop(), k))
    return sorted(pairs)


def pair_f1(pred, true):
    pred, true = set(pred), set(true)
    if not pred or not true:
        return 0.0
    tp = len(pred & true)
    precision, recall = tp / len(pred), tp / len(true)
    return 0.0 if tp == 0 else 2 * precision * recall / (precision + recall)
