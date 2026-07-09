"""cyclic_snf_table.py — generator of the "cyclic Smith-normal-form table
for n <= 12" promised in the paper's Code-and-data-availability paragraph.

Added to the public repository 2026-07-09 (the table was promised at
publication but not shipped until then — see CLAIM_LEDGER.md GAP-1).

For the cyclic Latin square L_ij = (i+j) mod n + 1 and each order
n = 2..12, computes the exact Smith normal form over Z of the
(n-1)x(n-1) difference matrix A (last column as reference) and checks
the paper's statement (Sec. 3 / SNF section): the invariant factors are

    (1, n, n, ..., n)          with n-2 copies of n,

hence |det(A)| = n^{n-2} and |det(E_std)| = n * |det(A)| = n^{n-1},
matching prop:cyclic. Writes results/certified/cyclic_snf_table.json.

Exact integer arithmetic throughout; deterministic (no randomness).
"""
import sys as _sys
try:
    _sys.stdout.reconfigure(encoding="utf-8")
    _sys.stderr.reconfigure(encoding="utf-8")
except Exception:
    pass

import json
import os


def cyclic_latin(n):
    return [[(i + j) % n + 1 for j in range(n)] for i in range(n)]


def compute_A(L):
    n = len(L)
    return [[L[i][j] - L[i][n - 1] for j in range(n - 1)]
            for i in range(n - 1)]


def smith_normal_form(M):
    """Invariant factors of an integer matrix (exact, in-place algorithm)."""
    A = [row[:] for row in M]
    rows, cols = len(A), len(A[0])
    factors = []
    t = 0
    while t < min(rows, cols):
        # find nonzero entry of minimal absolute value in A[t:][t:]
        pivot = None
        for i in range(t, rows):
            for j in range(t, cols):
                if A[i][j] != 0 and (pivot is None
                                     or abs(A[i][j]) < abs(A[pivot[0]][pivot[1]])):
                    pivot = (i, j)
        if pivot is None:
            break
        i0, j0 = pivot
        A[t], A[i0] = A[i0], A[t]
        for row in A:
            row[t], row[j0] = row[j0], row[t]
        # eliminate; restart if a smaller remainder appears
        while True:
            done = True
            for i in range(t + 1, rows):
                q = A[i][t] // A[t][t]
                if q != 0:
                    for j in range(t, cols):
                        A[i][j] -= q * A[t][j]
                if A[i][t] != 0:
                    A[t], A[i] = A[i], A[t]
                    done = False
            for j in range(t + 1, cols):
                q = A[t][j] // A[t][t]
                if q != 0:
                    for i in range(t, rows):
                        A[i][j] -= q * A[i][t]
                if A[t][j] != 0:
                    for row in A:
                        row[t], row[j] = row[j], row[t]
                    done = False
            if done:
                break
        # enforce divisibility d_t | A[i][j] for the rest
        redo = False
        for i in range(t + 1, rows):
            for j in range(t + 1, cols):
                if A[i][j] % A[t][t] != 0:
                    for k in range(t, cols):
                        A[t][k] += A[i][k]
                    redo = True
                    break
            if redo:
                break
        if redo:
            continue
        factors.append(abs(A[t][t]))
        t += 1
    return factors


def main():
    results = {}
    all_ok = True
    print("=== Cyclic SNF table, n = 2..12 ===")
    for n in range(2, 13):
        A = compute_A(cyclic_latin(n))
        inv = smith_normal_form(A)
        expected = [1] + [n] * (n - 2)
        det_abs = 1
        for d in inv:
            det_abs *= d
        ok = (inv == expected) and (det_abs == n ** (n - 2))
        all_ok = all_ok and ok
        status = "PASS" if ok else "FAIL"
        print(f"  [{status}] n={n:2d}: invariant factors = {inv}"
              f" (|det A| = n^{{n-2}} = {n ** (n - 2)})")
        results[str(n)] = {
            "invariant_factors": inv,
            "expected": expected,
            "abs_det_A": det_abs,
            "abs_det_E_std": n * det_abs,
            "matches_paper": ok,
        }

    out_dir = os.path.join(os.path.dirname(__file__), "..", "results",
                           "certified")
    os.makedirs(out_dir, exist_ok=True)
    out_path = os.path.join(out_dir, "cyclic_snf_table.json")
    with open(out_path, "w") as f:
        json.dump({
            "table": "cyclic Smith-normal-form table, n <= 12",
            "generator": "scripts/cyclic_snf_table.py",
            "convention": "A[i][j] = L[i][j] - L[i][n-1], (n-1)x(n-1)",
            "results": results,
            "all_match_paper": all_ok,
        }, f, indent=2)
    print(f"\nTable written: {out_path}")
    print("ALL ROWS MATCH PAPER" if all_ok else "MISMATCH FOUND")
    return 0 if all_ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
