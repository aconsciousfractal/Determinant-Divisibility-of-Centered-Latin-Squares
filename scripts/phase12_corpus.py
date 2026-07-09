#!/usr/bin/env python3
"""
Phase 12 — Counterexample Corpus & Feature Analysis

Research-phase script, added to the public repository 2026-07-09: this is
the GENERATOR of the tracked `results/phase12/` datasets (the paper's
"4,000-sample n=10 switch-chain table"). NOTE the actual seed is
SEED = 20260412 below; the paper's availability paragraph says
"seed 42", which is incorrect — see CLAIM_LEDGER.md GAP-1.

Pipeline:
  1. Switch-chain MCMC on R(10,5) — sample doubly balanced binary matrices
  2. F2-rank check — classify full-rank (counterexample pattern) vs deficient
  3. Feature extraction (Levels 1-3) for every sample
  4. König lift for full-rank patterns → verify det(A) odd
  5. Save per-sample matrices + features

Targets: ≥500 full-rank patterns + ≥2000 rank-deficient controls
At 17.8% full-rank rate, ~3500 total samples suffice.
"""

import numpy as np
from collections import Counter, defaultdict
from fractions import Fraction
import json, time, sys, os

# ═══════════════════════════════════════════════════════════════════
# CONSTANTS
# ═══════════════════════════════════════════════════════════════════
N, K = 10, 5
N_SAMPLES = 4000          # expect ~710 full-rank, ~3290 deficient
BURN_IN = 50_000
THIN = 200
SEED = 2026_04_12
RNG = np.random.default_rng(SEED)

OUT_DIR = os.path.join(os.path.dirname(__file__), '..', 'results', 'phase12')
os.makedirs(OUT_DIR, exist_ok=True)

# ═══════════════════════════════════════════════════════════════════
# SWITCH-CHAIN MCMC
# ═══════════════════════════════════════════════════════════════════

def make_start():
    """Circulant circ({0,1,2,3,4}) as starting matrix."""
    P = np.zeros((N, N), dtype=np.int8)
    for i in range(N):
        for d in range(K):
            P[i, (i + d) % N] = 1
    return P

def switch(P, rng):
    """One switch step on the binary matrix (preserves row/col sums)."""
    i, j = rng.choice(N, 2, replace=False)
    c, d = rng.choice(N, 2, replace=False)
    a, b, e, f = P[i, c], P[i, d], P[j, c], P[j, d]
    if a == 1 and f == 1 and b == 0 and e == 0:
        P[i, c], P[i, d], P[j, c], P[j, d] = 0, 1, 1, 0
    elif b == 1 and e == 1 and a == 0 and f == 0:
        P[i, c], P[i, d], P[j, c], P[j, d] = 1, 0, 0, 1

# ═══════════════════════════════════════════════════════════════════
# F2 LINEAR ALGEBRA
# ═══════════════════════════════════════════════════════════════════

def rank_f2(P):
    """F2-rank of an integer matrix."""
    M = P.copy() % 2
    n_rows, n_cols = M.shape
    r = 0
    for c in range(n_cols):
        piv = None
        for row in range(r, n_rows):
            if M[row, c] % 2:
                piv = row; break
        if piv is None:
            continue
        M[[r, piv]] = M[[piv, r]]
        pv = M[r].copy()
        for row in range(n_rows):
            if row != r and M[row, c] % 2:
                M[row] = (M[row] + pv) % 2
        r += 1
    return r

def nullspace_f2(P):
    """Return basis vectors of ker(P^T) over F2."""
    # We want left null-space of P, i.e. ker(P^T)
    M = P.T.copy() % 2
    n_rows, n_cols = M.shape
    pivots = []
    r = 0
    for c in range(n_cols):
        piv = None
        for row in range(r, n_rows):
            if M[row, c] % 2:
                piv = row; break
        if piv is None:
            continue
        M[[r, piv]] = M[[piv, r]]
        pv = M[r].copy()
        for row in range(n_rows):
            if row != r and M[row, c] % 2:
                M[row] = (M[row] + pv) % 2
        pivots.append(c)
        r += 1
    free = [c for c in range(n_cols) if c not in pivots]
    basis = []
    for fc in free:
        v = np.zeros(n_cols, dtype=np.int8)
        v[fc] = 1
        for pi, pc in enumerate(pivots):
            v[pc] = M[pi, fc] % 2
        basis.append(v)
    return basis

# ═══════════════════════════════════════════════════════════════════
# KÖNIG LIFT: balanced P -> Latin square
# ═══════════════════════════════════════════════════════════════════

def _aug(u, adj, mr, vis):
    for v in adj[u]:
        if not vis[v]:
            vis[v] = True
            if mr[v] < 0 or _aug(mr[v], adj, mr, vis):
                mr[v] = u; return True
    return False

def perfect_matching(adj, n):
    mr = [-1] * n
    for u in range(n):
        _aug(u, adj, mr, [False] * n)
    ml = [-1] * n
    for j in range(n):
        if mr[j] >= 0:
            ml[mr[j]] = j
    return ml if all(x >= 0 for x in ml) else None

def decompose_matchings(P_list):
    n = len(P_list)
    rem = [r[:] for r in P_list]
    ms = []
    while any(rem[0]):
        adj = [[j for j in range(n) if rem[i][j]] for i in range(n)]
        m = perfect_matching(adj, n)
        assert m is not None, "decomposition failed"
        ms.append(m)
        for i in range(n):
            rem[i][m[i]] -= 1
    return ms

def lift_to_latin(P_np):
    """P is numpy array. Returns list-of-lists Latin square, symbols 1..n."""
    n = P_np.shape[0]
    P_list = P_np.tolist()
    Pc = [[1 - P_list[i][j] for j in range(n)] for i in range(n)]
    mo = decompose_matchings(P_list)
    me = decompose_matchings(Pc)
    osym = list(range(1, n + 1, 2))
    esym = list(range(2, n + 1, 2))
    L = [[0] * n for _ in range(n)]
    for k, m in enumerate(mo):
        for i in range(n):
            L[i][m[i]] = osym[k]
    for k, m in enumerate(me):
        for i in range(n):
            L[i][m[i]] = esym[k]
    return L

# ═══════════════════════════════════════════════════════════════════
# EXACT INTEGER DETERMINANT (Bareiss via Fraction)
# ═══════════════════════════════════════════════════════════════════

def det_exact(M):
    n = len(M)
    A = [[Fraction(M[i][j]) for j in range(n)] for i in range(n)]
    s = 1
    for c in range(n):
        pv = next((r for r in range(c, n) if A[r][c]), None)
        if pv is None:
            return 0
        if pv != c:
            A[c], A[pv] = A[pv], A[c]; s = -s
        for r in range(c + 1, n):
            if A[r][c]:
                f = A[r][c] / A[c][c]
                for j in range(n):
                    A[r][j] -= f * A[c][j]
    d = s
    for i in range(n):
        d *= A[i][i]
    return int(d)

# ═══════════════════════════════════════════════════════════════════
# FEATURE EXTRACTION
# ═══════════════════════════════════════════════════════════════════

def extract_features(P):
    """Compute Level 1-3 features for a 10x10 doubly balanced matrix P."""
    n = P.shape[0]
    feat = {}

    # ── Level 1: Pattern-level ──
    rk = rank_f2(P)
    feat['rank'] = int(rk)

    col_weights = sorted(int(P[:, j].sum()) for j in range(n))
    feat['col_weights'] = col_weights

    # Overlap matrix (over Z)
    ov_mat = (P @ P.T)
    ov_upper = []
    for i in range(n):
        for j in range(i + 1, n):
            ov_upper.append(int(ov_mat[i, j]))
    feat['min_ov'] = min(ov_upper)
    feat['max_ov'] = max(ov_upper)
    feat['n_ov1'] = ov_upper.count(1)
    feat['n_ov2'] = ov_upper.count(2)
    feat['n_ov3'] = ov_upper.count(3)
    feat['n_ov4'] = ov_upper.count(4)
    feat['var_ov'] = float(np.var(ov_upper))
    feat['ov_histogram'] = [ov_upper.count(v) for v in range(K + 1)]

    # sigma_min (over R)
    sv = np.linalg.svd(P.astype(float), compute_uv=False)
    feat['sigma_min'] = float(sv[-1])
    feat['sigma_2nd'] = float(sv[-2]) if len(sv) > 1 else 0.0

    # ── Level 2: Code-level (on B = A mod 2) ──
    # B_{ij} = P_{ij} XOR P_{i, n-1}  (Lemma 4.2)
    last_col = P[:, n - 1].copy()
    B = np.zeros((n - 1, n - 1), dtype=np.int8)
    for i in range(n - 1):
        for j in range(n - 1):
            B[i, j] = (P[i, j] ^ last_col[i]) % 2

    # Hamming distances between rows of B
    hamming = []
    for i in range(n - 1):
        for j in range(i + 1, n - 1):
            hamming.append(int(np.sum(B[i] != B[j])))
    feat['d_min'] = min(hamming) if hamming else 0
    feat['d_max'] = max(hamming) if hamming else 0
    feat['d_mean'] = float(np.mean(hamming)) if hamming else 0.0

    # Weight enumerator of row code of B
    row_weights = [int(np.sum(B[i])) for i in range(n - 1)]
    feat['B_row_weights'] = sorted(row_weights)

    # Spectrum of BB^T over Z (eigenvalues)
    G = B.astype(np.float64) @ B.astype(np.float64).T
    eigs = sorted(np.linalg.eigvalsh(G))
    feat['gram_eigs'] = [round(float(e), 6) for e in eigs]
    feat['gram_eig_min'] = float(eigs[0])

    # Null-space analysis for rank-deficient
    if rk < n:
        ns = nullspace_f2(P)
        if ns:
            feat['null_dim'] = len(ns)
            feat['null_weights'] = sorted(int(np.sum(v)) for v in ns)
            feat['null_min_weight'] = min(feat['null_weights'])
        else:
            feat['null_dim'] = 0
            feat['null_weights'] = []
            feat['null_min_weight'] = None
    else:
        feat['null_dim'] = 0
        feat['null_weights'] = []
        feat['null_min_weight'] = None

    # ── Level 3: Graph-level ──
    # Overlap graph: nodes = rows, edge weight = ov(i,j)
    # Laplacian: L = D - W  where W_{ij} = ov(i,j), D_{ii} = sum_j W_{ij}
    W = ov_mat.astype(np.float64).copy()
    np.fill_diagonal(W, 0)
    D = np.diag(W.sum(axis=1))
    Lap = D - W
    lap_eigs = sorted(np.linalg.eigvalsh(Lap))
    feat['fiedler'] = float(lap_eigs[1])  # 2nd smallest eigenvalue
    feat['lap_eig3'] = float(lap_eigs[2]) if len(lap_eigs) > 2 else 0.0

    # Binary overlap graph: edge iff ov(i,j) >= 1
    # Cheeger-like: algebraic connectivity of unweighted version
    W_bin = (ov_mat > 0).astype(np.float64)
    np.fill_diagonal(W_bin, 0)
    D_bin = np.diag(W_bin.sum(axis=1))
    Lap_bin = D_bin - W_bin
    lap_bin_eigs = sorted(np.linalg.eigvalsh(Lap_bin))
    feat['fiedler_bin'] = float(lap_bin_eigs[1])

    return feat

# ═══════════════════════════════════════════════════════════════════
# MAIN
# ═══════════════════════════════════════════════════════════════════

def main():
    print("=" * 72)
    print("  PHASE 12 — COUNTEREXAMPLE CORPUS & FEATURE ANALYSIS")
    print("=" * 72)
    print(f"  N={N}, K={K}, samples={N_SAMPLES}, burn_in={BURN_IN}, thin={THIN}")
    print(f"  Seed={SEED}")
    print(f"  Output: {OUT_DIR}")
    print()

    # ── Step 1: Burn-in ──
    t0 = time.time()
    print(f"  Burn-in ({BURN_IN} steps)...", end='', flush=True)
    P = make_start()
    for _ in range(BURN_IN):
        switch(P, RNG)
    print(f" done ({time.time()-t0:.1f}s)")

    # ── Step 2: Sample + extract features ──
    full_rank = []   # list of (matrix, features)
    deficient = []
    rank_dist = Counter()

    t1 = time.time()
    for s in range(N_SAMPLES):
        for _ in range(THIN):
            switch(P, RNG)
        Pc = P.copy()
        feat = extract_features(Pc)
        rank_dist[feat['rank']] += 1

        record = {
            'matrix': Pc.tolist(),
            'features': feat,
        }

        if feat['rank'] == N:
            full_rank.append(record)
        else:
            deficient.append(record)

        done = s + 1
        if done % 200 == 0 or done == N_SAMPLES:
            elapsed = time.time() - t1
            rate = done / elapsed
            n_fr = len(full_rank)
            print(f"\r  Sampling: {done}/{N_SAMPLES}  "
                  f"full-rank={n_fr} ({100*n_fr/done:.1f}%)  "
                  f"{rate:.1f} samples/s  ", end='', flush=True)

    print()
    elapsed_sample = time.time() - t1
    print(f"\n  Sampling complete: {time.time()-t0:.1f}s total")
    print(f"  Full-rank: {len(full_rank)}")
    print(f"  Deficient: {len(deficient)}")
    print(f"  Rank distribution: {dict(sorted(rank_dist.items()))}")

    # ── Step 3: König lift for full-rank patterns ──
    print(f"\n  Lifting {len(full_rank)} full-rank patterns to Latin squares...")
    t2 = time.time()
    n_ce = 0
    for idx, rec in enumerate(full_rank):
        P_np = np.array(rec['matrix'], dtype=np.int8)
        L = lift_to_latin(P_np)

        # Verify parity preservation: L mod 2 == P
        L_np = np.array(L, dtype=np.int64)
        assert np.array_equal(L_np % 2, P_np), f"Parity mismatch at idx {idx}"

        # Compute A_{ij} = L_{ij} - L_{i,n-1}  (column-reduced matrix)
        A = [[L[i][j] - L[i][N - 1] for j in range(N - 1)] for i in range(N - 1)]
        d = det_exact(A)
        rec['features']['det_A'] = d
        rec['features']['det_A_odd'] = (d % 2 != 0)
        rec['features']['det_A_div5'] = (d % 5 == 0)
        det_E = N * d
        rec['features']['det_E'] = det_E
        rec['features']['n2_divides'] = (det_E % (N * N) == 0)
        rec['features']['is_counterexample'] = not rec['features']['n2_divides']

        if rec['features']['is_counterexample']:
            n_ce += 1

        # Store Latin square (compact: just the matrix)
        rec['latin_square'] = L

        if (idx + 1) % 100 == 0 or idx + 1 == len(full_rank):
            print(f"\r    Lifted {idx+1}/{len(full_rank)}, "
                  f"counterexamples={n_ce}  ", end='', flush=True)

    print()
    print(f"  Counterexamples: {n_ce}/{len(full_rank)} "
          f"({100*n_ce/len(full_rank) if full_rank else 0:.1f}%)")
    print(f"  Lift time: {time.time()-t2:.1f}s")

    # ── Step 4: Summary statistics ──
    print("\n" + "=" * 72)
    print("  SUMMARY STATISTICS")
    print("=" * 72)

    # Feature comparison: full-rank vs false-positive vs true-negative
    fr_records = [r['features'] for r in full_rank]
    # False positives: min_ov >= 1 but rank < N
    fp_records = [r['features'] for r in deficient if r['features']['min_ov'] >= 1]
    # True negatives: min_ov == 0
    tn_records = [r['features'] for r in deficient if r['features']['min_ov'] == 0]

    print(f"\n  Populations:")
    print(f"    Full-rank (counterexample patterns): {len(fr_records)}")
    print(f"    False-positive (min_ov>=1, rank<N):  {len(fp_records)}")
    print(f"    True-negative (min_ov=0):            {len(tn_records)}")

    # Compare key features
    def stats(values, name):
        if not values:
            return f"    {name}: (empty)"
        a = np.array(values)
        return (f"    {name}: mean={np.mean(a):.4f}  "
                f"std={np.std(a):.4f}  "
                f"min={np.min(a):.4f}  max={np.max(a):.4f}")

    for label, recs in [("Full-rank", fr_records),
                         ("False-positive", fp_records),
                         ("True-negative", tn_records)]:
        if not recs:
            continue
        print(f"\n  ── {label} (n={len(recs)}) ──")
        print(stats([r['sigma_min'] for r in recs], "sigma_min"))
        print(stats([r['fiedler'] for r in recs], "fiedler  "))
        print(stats([r['fiedler_bin'] for r in recs], "fiedler_b"))
        print(stats([r['d_min'] for r in recs], "d_min    "))
        print(stats([r['var_ov'] for r in recs], "var_ov   "))
        print(stats([r['n_ov1'] for r in recs], "n_ov1    "))
        print(stats([r['gram_eig_min'] for r in recs], "gram_emin"))
        nw = [r.get('null_min_weight') for r in recs if r.get('null_min_weight') is not None]
        if nw:
            print(stats(nw, "null_wmin"))

    # ── H3: Fiedler value as discriminator ──
    print(f"\n  ── H3: FIEDLER VALUE AS FULL-RANK DISCRIMINATOR ──")
    if fr_records and fp_records:
        fr_f = np.array([r['fiedler'] for r in fr_records])
        fp_f = np.array([r['fiedler'] for r in fp_records])
        # Try thresholds
        print(f"  {'Thresh':>7}  {'TP':>5}  {'FP':>5}  {'FN':>5}  {'Prec':>7}  {'Rec':>7}  {'F1':>7}")
        best_f1 = 0
        best_t = 0
        for t_pct in range(0, 101, 5):
            t = np.percentile(np.concatenate([fr_f, fp_f]), t_pct)
            tp = int(np.sum(fr_f >= t))
            fp_count = int(np.sum(fp_f >= t))
            fn = int(np.sum(fr_f < t))
            prec = tp / (tp + fp_count) if (tp + fp_count) > 0 else 0
            rec = tp / (tp + fn) if (tp + fn) > 0 else 0
            f1 = 2 * prec * rec / (prec + rec) if (prec + rec) > 0 else 0
            if f1 > best_f1:
                best_f1 = f1; best_t = t
            print(f"  {t:>7.3f}  {tp:>5}  {fp_count:>5}  {fn:>5}  "
                  f"{prec:>7.3f}  {rec:>7.3f}  {f1:>7.3f}")
        print(f"\n  Best F1={best_f1:.3f} at threshold={best_t:.4f}")

    # ── sigma_min as discriminator (baseline) ──
    print(f"\n  ── BASELINE: sigma_min AS FULL-RANK DISCRIMINATOR ──")
    if fr_records and fp_records:
        fr_s = np.array([r['sigma_min'] for r in fr_records])
        fp_s = np.array([r['sigma_min'] for r in fp_records])
        print(f"  {'Thresh':>7}  {'TP':>5}  {'FP':>5}  {'FN':>5}  {'Prec':>7}  {'Rec':>7}  {'F1':>7}")
        best_f1_s = 0
        best_t_s = 0
        for t in [x * 0.05 for x in range(1, 20)]:
            tp = int(np.sum(fr_s >= t))
            fp_count = int(np.sum(fp_s >= t))
            fn = int(np.sum(fr_s < t))
            prec = tp / (tp + fp_count) if (tp + fp_count) > 0 else 0
            rec = tp / (tp + fn) if (tp + fn) > 0 else 0
            f1 = 2 * prec * rec / (prec + rec) if (prec + rec) > 0 else 0
            if f1 > best_f1_s:
                best_f1_s = f1; best_t_s = t
            print(f"  {t:>7.3f}  {tp:>5}  {fp_count:>5}  {fn:>5}  "
                  f"{prec:>7.3f}  {rec:>7.3f}  {f1:>7.3f}")
        print(f"\n  Best F1={best_f1_s:.3f} at threshold={best_t_s:.4f}")

    # ── H5: null min weight for false positives ──
    print(f"\n  ── H5: NULL SPACE MIN WEIGHT (false positives only) ──")
    if fp_records:
        nw = [r.get('null_min_weight') for r in fp_records
              if r.get('null_min_weight') is not None]
        if nw:
            c = Counter(nw)
            print(f"  Weight distribution: {dict(sorted(c.items()))}")
            print(f"  Mean={np.mean(nw):.2f}, min={min(nw)}, max={max(nw)}")

    # ── Counterexample properties ──
    print(f"\n  ── COUNTEREXAMPLE PROPERTIES ──")
    ce_recs = [r for r in full_rank if r['features'].get('is_counterexample')]
    if ce_recs:
        dets = [r['features']['det_A'] for r in ce_recs]
        print(f"  Total counterexamples: {len(ce_recs)}")
        print(f"  All det(A) odd: {all(d % 2 != 0 for d in dets)}")
        print(f"  All 5|det(A): {all(d % 5 == 0 for d in dets)}")
        print(f"  5|det(A): {sum(1 for d in dets if d % 5 == 0)}/{len(dets)}")
        # v_2 and v_5 distributions
        def vp(n, p):
            if n == 0: return float('inf')
            c, n = 0, abs(n)
            while n % p == 0: n //= p; c += 1
            return c
        v2_vals = sorted(set(vp(d, 2) for d in dets))
        v5_vals = sorted(set(vp(d, 5) for d in dets))
        print(f"  v₂(det(A)) values: {v2_vals}")
        print(f"  v₅(det(A)) values: {v5_vals}")
        print(f"  |det(A)| range: [{min(abs(d) for d in dets)}, {max(abs(d) for d in dets)}]")
        # Column weight profiles
        cw_profiles = Counter(tuple(r['features']['col_weights']) for r in ce_recs)
        print(f"  Distinct column weight profiles: {len(cw_profiles)}")
        for prof, cnt in cw_profiles.most_common(10):
            print(f"    {prof}: {cnt}")

    # ── Step 5: Save ──
    print(f"\n  Saving results to {OUT_DIR}...")

    # Save compact summary (features only, no matrices)
    summary = {
        'params': {'N': N, 'K': K, 'N_SAMPLES': N_SAMPLES,
                   'BURN_IN': BURN_IN, 'THIN': THIN, 'SEED': SEED},
        'counts': {
            'total': N_SAMPLES,
            'full_rank': len(full_rank),
            'deficient': len(deficient),
            'false_positive': len(fp_records),
            'true_negative': len(tn_records),
            'counterexamples': n_ce,
        },
        'rank_distribution': dict(sorted(rank_dist.items())),
    }
    with open(os.path.join(OUT_DIR, 'summary.json'), 'w') as f:
        json.dump(summary, f, indent=2, default=str)

    # Save full-rank corpus (matrices + features + LS)
    # Exclude numpy arrays, convert to lists
    def clean(rec):
        r = {}
        r['matrix'] = rec['matrix']
        r['features'] = {k: (v if not isinstance(v, np.ndarray) else v.tolist())
                         for k, v in rec['features'].items()}
        if 'latin_square' in rec:
            r['latin_square'] = rec['latin_square']
        return r

    with open(os.path.join(OUT_DIR, 'full_rank_corpus.json'), 'w') as f:
        json.dump([clean(r) for r in full_rank], f)
    print(f"    full_rank_corpus.json: {len(full_rank)} records")

    # Save deficient controls (matrices + features, no LS)
    with open(os.path.join(OUT_DIR, 'deficient_controls.json'), 'w') as f:
        json.dump([clean(r) for r in deficient], f)
    print(f"    deficient_controls.json: {len(deficient)} records")

    # Save feature matrix as CSV for quick analysis
    all_recs = full_rank + deficient
    csv_path = os.path.join(OUT_DIR, 'features.csv')
    with open(csv_path, 'w') as f:
        header = ['rank', 'min_ov', 'max_ov', 'n_ov1', 'n_ov2', 'n_ov3',
                  'n_ov4', 'var_ov', 'sigma_min', 'sigma_2nd', 'd_min',
                  'd_max', 'd_mean', 'gram_eig_min', 'fiedler', 'fiedler_bin',
                  'lap_eig3', 'null_dim', 'null_min_weight', 'group']
        f.write(','.join(header) + '\n')
        for rec in all_recs:
            feat = rec['features']
            group = 'FR' if feat['rank'] == N else (
                'FP' if feat['min_ov'] >= 1 else 'TN')
            row = [
                feat['rank'], feat['min_ov'], feat['max_ov'],
                feat['n_ov1'], feat['n_ov2'], feat['n_ov3'], feat['n_ov4'],
                f"{feat['var_ov']:.6f}", f"{feat['sigma_min']:.6f}",
                f"{feat['sigma_2nd']:.6f}", feat['d_min'], feat['d_max'],
                f"{feat['d_mean']:.4f}", f"{feat['gram_eig_min']:.6f}",
                f"{feat['fiedler']:.6f}", f"{feat['fiedler_bin']:.6f}",
                f"{feat['lap_eig3']:.6f}", feat.get('null_dim', 0),
                feat.get('null_min_weight', ''), group,
            ]
            f.write(','.join(str(v) for v in row) + '\n')
    print(f"    features.csv: {len(all_recs)} rows")

    total_time = time.time() - t0
    print(f"\n  Total time: {total_time:.1f}s ({total_time/60:.1f}min)")
    print("  DONE.")


if __name__ == '__main__':
    main()
