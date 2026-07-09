"""jm_n8_corpus.py — regeneration of the "10,000-sample n=8 table"
promised in the paper's Code-and-data-availability paragraph.

PROVENANCE NOTE (2026-07-09): the original corpus behind the paper's
rem:n8-sampling ("approximately 11% satisfy 32 | det(E_std) but
64 does not divide det(E_std)") did not survive as a dataset or script;
this script REGENERATES a corpus to the paper's stated specification
(10,000 Jacobson-Matthews samples at n = 8, 5*n^2 = 320 mixing steps
each, seed 0) using the same JM sampler implementation as
scripts/jm_paper_verify.py. The paper's ~11% figure is explicitly
approximate ("not claimed as a rigorous proportion"); this regeneration
checks it, it does not reproduce the original byte-for-byte. See
CLAIM_LEDGER.md GAP-1.

Also locates the first SHARPNESS WITNESS (v2(det E_std) = 5, i.e.
32 | det(E_std) but 64 does not divide), certifying that the n^2/2
bound is attained at n = 8. By thm:sharp-v2, every sharp sample
necessarily has dim ker_F2(A mod 2) = 1 (asserted in code for the
recorded witness; proved in-paper for all). NOTE: the paper's
proof describes its witness as "the first JM sample (seed 0) satisfying
dim ker_F2 = 1"; that selection is specific to the original sampler
implementation and is NOT portably reproducible — with this
implementation the first dim-ker-1 sample happens not to be sharp
(recorded in the summary). The mathematical content (a sharpness
witness exists and is verified exactly) is what this script certifies.

Writes results/certified/jm_n8_corpus_summary.json.
"""
import sys as _sys
try:
    _sys.stdout.reconfigure(encoding="utf-8")
    _sys.stderr.reconfigure(encoding="utf-8")
except Exception:
    pass

import json
import os
import random
import time

N = 8
N_SAMPLES = 10_000
STEPS = 5 * N * N  # 320
SEED = 0


def jm_sample(n, steps):
    """Jacobson-Matthews sampler (same implementation as jm_paper_verify.py)."""
    cube = [[[0] * n for _ in range(n)] for _ in range(n)]
    for i in range(n):
        for j in range(n):
            cube[i][j][(i + j) % n] = 1
    proper = True
    ir = ic = ik = 0
    count = 0
    while count < steps or not proper:
        if proper:
            i = random.randint(0, n - 1)
            j = random.randint(0, n - 1)
            row_ij = cube[i][j]
            k_cur = 0
            while row_ij[k_cur] != 1:
                k_cur += 1
            k = random.randint(0, n - 2)
            if k >= k_cur:
                k += 1
            count += 1
        else:
            i, j, k = ir, ic, ik
            cands = [x for x in range(n) if cube[i][j][x] == 1]
            k_cur = cands[random.randint(0, len(cands) - 1)]
        jc = [jj for jj in range(n) if jj != j and cube[i][jj][k] == 1]
        j_p = jc[random.randint(0, len(jc) - 1)]
        ic_ = [ii for ii in range(n) if ii != i and cube[ii][j][k] == 1]
        i_p = ic_[random.randint(0, len(ic_) - 1)]
        cube[i][j][k] += 1
        cube[i][j][k_cur] -= 1
        cube[i][j_p][k] -= 1
        cube[i][j_p][k_cur] += 1
        cube[i_p][j][k] -= 1
        cube[i_p][j][k_cur] += 1
        cube[i_p][j_p][k] += 1
        cube[i_p][j_p][k_cur] -= 1
        if cube[i_p][j_p][k_cur] == -1:
            proper = False
            ir, ic, ik = i_p, j_p, k_cur
        else:
            proper = True
    L = [[0] * n for _ in range(n)]
    for i in range(n):
        for j in range(n):
            for k in range(n):
                if cube[i][j][k] == 1:
                    L[i][j] = k + 1
                    break
    return L


def compute_A(L):
    n = len(L)
    return [[L[i][j] - L[i][n - 1] for j in range(n - 1)]
            for i in range(n - 1)]


def bareiss_det(M):
    n = len(M)
    A = [row[:] for row in M]
    sign = 1
    prev = 1
    for i in range(n):
        if A[i][i] == 0:
            found = False
            for k in range(i + 1, n):
                if A[k][i] != 0:
                    A[i], A[k] = A[k], A[i]
                    sign *= -1
                    found = True
                    break
            if not found:
                return 0
        for j in range(i + 1, n):
            for k in range(i + 1, n):
                A[j][k] = (A[j][k] * A[i][i] - A[j][i] * A[i][k]) // prev
            A[j][i] = 0
        prev = A[i][i]
    return sign * A[n - 1][n - 1]


def kernel_dim_F2(M):
    B = [[x % 2 for x in row] for row in M]
    rows, cols = len(B), len(B[0])
    pr = 0
    for col in range(cols):
        pivot = None
        for row in range(pr, rows):
            if B[row][col] == 1:
                pivot = row
                break
        if pivot is None:
            continue
        B[pr], B[pivot] = B[pivot], B[pr]
        for row in range(rows):
            if row != pr and B[row][col] == 1:
                B[row] = [(a + b) % 2 for a, b in zip(B[row], B[pr])]
        pr += 1
    return cols - pr


def v2(x):
    if x == 0:
        return None  # singular; excluded from valuation stats
    x = abs(x)
    v = 0
    while x % 2 == 0:
        x //= 2
        v += 1
    return v


def main():
    random.seed(SEED)
    t0 = time.time()
    v2_hist = {}
    n_singular = 0
    n_sharp = 0  # v2(det E_std) = 5, i.e. 32 | det, 64 does not divide
    n_bound_violation = 0
    first_ker1_sample = None   # implementation-dependent (see docstring)
    first_sharp_witness = None  # the certified sharpness witness

    print(f"=== JM corpus: n={N}, {N_SAMPLES} samples, {STEPS} steps, "
          f"seed {SEED} ===")
    for s in range(N_SAMPLES):
        L = jm_sample(N, STEPS)
        A = compute_A(L)
        det_A = bareiss_det(A)
        det_E = N * det_A
        val = v2(det_E)
        if val is None:
            n_singular += 1
        else:
            v2_hist[val] = v2_hist.get(val, 0) + 1
            if val < 5:
                n_bound_violation += 1
            if val == 5:
                n_sharp += 1
                if first_sharp_witness is None:
                    k = kernel_dim_F2(A)
                    assert k == 1, (
                        f"sharp sample with dim ker = {k} contradicts "
                        f"thm:sharp-v2")
                    first_sharp_witness = {
                        "sample_index": s,
                        "L": L,
                        "det_A": det_A,
                        "det_E_std": det_E,
                        "v2_det_E_std": val,
                        "dim_ker_F2": k,
                    }
        if first_ker1_sample is None and val is not None:
            k = kernel_dim_F2(A)
            if k == 1:
                first_ker1_sample = {
                    "sample_index": s,
                    "v2_det_E_std": val,
                    "is_sharp": val == 5,
                }
        if (s + 1) % 1000 == 0:
            print(f"  {s + 1} samples, sharp so far: {n_sharp} "
                  f"({100.0 * n_sharp / (s + 1):.1f}%)")

    n_valued = N_SAMPLES - n_singular
    frac_sharp = n_sharp / n_valued if n_valued else 0.0
    elapsed = time.time() - t0

    print(f"\n  singular (det = 0): {n_singular}")
    print(f"  v2(det E_std) histogram: "
          f"{dict(sorted(v2_hist.items()))}")
    print(f"  bound violations (v2 < 5): {n_bound_violation}")
    print(f"  sharp samples (32 | det, 64 does not divide): {n_sharp} "
          f"/ {n_valued} = {100.0 * frac_sharp:.2f}%  (paper: ~11%)")
    if first_sharp_witness is not None:
        w = first_sharp_witness
        print(f"  first sharpness witness: index {w['sample_index']}, "
              f"det(E_std) = {w['det_E_std']}, dim ker_F2 = 1")
    if first_ker1_sample is not None:
        w = first_ker1_sample
        print(f"  first dim-ker-1 sample: index {w['sample_index']}, "
              f"v2(det E_std) = {w['v2_det_E_std']}, sharp: {w['is_sharp']}"
              f"  (implementation-dependent; see docstring)")

    out_dir = os.path.join(os.path.dirname(__file__), "..", "results",
                           "certified")
    os.makedirs(out_dir, exist_ok=True)
    out_path = os.path.join(out_dir, "jm_n8_corpus_summary.json")
    with open(out_path, "w") as f:
        json.dump({
            "table": "n=8 Jacobson-Matthews corpus (REGENERATED; see "
                     "docstring of scripts/jm_n8_corpus.py)",
            "params": {"n": N, "samples": N_SAMPLES, "steps": STEPS,
                       "seed": SEED},
            "singular": n_singular,
            "v2_det_E_std_histogram":
                {str(k): v for k, v in sorted(v2_hist.items())},
            "bound_violations_v2_below_5": n_bound_violation,
            "sharp_count": n_sharp,
            "sharp_fraction": frac_sharp,
            "paper_claim": "approximately 11% (rem:n8-sampling; "
                           "explicitly non-rigorous)",
            "first_sharpness_witness": first_sharp_witness,
            "first_dim_ker_1_sample": first_ker1_sample,
            "note": "the paper's 'first dim-ker-1 sample is the witness' "
                    "selection is sampler-implementation-specific; this "
                    "regeneration certifies the mathematical content (a "
                    "verified sharpness witness exists) instead",
        }, f, indent=2)
    print(f"\nSummary written: {out_path}")
    print(f"Total time: {elapsed:.1f}s")

    ok = n_bound_violation == 0 and first_sharp_witness is not None
    print("ALL CHECKS PASS" if ok else "CHECK OUTCOME: see summary "
          "(bound violation or no sharpness witness found)")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
