"""
sandpile_n6_exhaustive.py — Exhaustive n=6 analysis (9408 reduced LS)
Tasks: 
  - SNF profiles + distribution
  - C7: |K(A)| | |K(LSG)| 
  - C9: parity vs K(A) structure
  - C7 refined conjecture: failure iff prime support mismatch

Output saved to n6_results.txt
"""
import numpy as np
from math import gcd
from functools import reduce
from collections import Counter
import time, sys

# ═══════════════════════════════════════════════════════════════════════
# Core routines
# ═══════════════════════════════════════════════════════════════════════

def bareiss_det(M_np):
    """Exact integer determinant via Bareiss algorithm."""
    n = M_np.shape[0]
    if n == 0: return 1
    A = [[int(M_np[i][j]) for j in range(n)] for i in range(n)]
    sign = 1; prev = 1
    for i in range(n):
        if A[i][i] == 0:
            found = False
            for k in range(i+1, n):
                if A[k][i] != 0:
                    A[i], A[k] = A[k], A[i]; sign *= -1; found = True; break
            if not found: return 0
        for j in range(i+1, n):
            for k in range(i+1, n):
                A[j][k] = (A[j][k]*A[i][i] - A[j][i]*A[i][k]) // prev
            A[j][i] = 0
        prev = A[i][i]
    return sign * A[n-1][n-1]


def get_A(L):
    """Reduced matrix A_{ij} = L_{ij} - L_{i,n-1}."""
    n = len(L)
    return [[L[i][j] - L[i][n-1] for j in range(n-1)] for i in range(n-1)]


def snf_diag(A_int):
    """Smith Normal Form — returns diagonal entries (including zeros)."""
    A = [list(row) for row in A_int]
    m = len(A); n = len(A[0]) if m else 0
    diag = []
    for k in range(min(m, n)):
        found = True
        while found:
            found = False
            min_val = None; mi, mj = k, k
            for i in range(k, m):
                for j in range(k, n):
                    if A[i][j] != 0 and (min_val is None or abs(A[i][j]) < abs(min_val)):
                        min_val = A[i][j]; mi, mj = i, j
            if min_val is None: break
            if mi != k: A[k], A[mi] = A[mi], A[k]
            if mj != k:
                for i in range(m): A[i][k], A[i][mj] = A[i][mj], A[i][k]
            if A[k][k] < 0:
                for j in range(n): A[k][j] = -A[k][j]
            for i in range(k+1, m):
                if A[i][k] != 0:
                    q = A[i][k] // A[k][k]
                    for j in range(n): A[i][j] -= q * A[k][j]
                    if A[i][k] != 0: found = True
            for j in range(k+1, n):
                if A[k][j] != 0:
                    q = A[k][j] // A[k][k]
                    for i in range(m): A[i][j] -= q * A[i][k]
                    if A[k][j] != 0: found = True
            if not found:
                for i in range(k+1, m):
                    for j in range(k+1, n):
                        if A[k][k] != 0 and A[i][j] % A[k][k] != 0:
                            for jj in range(n): A[k][jj] += A[i][jj]
                            found = True; break
                    if found: break
        diag.append(abs(A[k][k]) if A[k][k] != 0 else 0)
    return diag


def snf_nonzero(A_int):
    """SNF invariant factors > 1."""
    d = snf_diag(A_int)
    return tuple(x for x in d if x > 1)


def parity_det(L):
    """Compute det(A) and return its sign/parity: +1, -1, or 0."""
    A = get_A(L)
    A_np = np.array(A, dtype=int)
    d = bareiss_det(A_np)
    if d > 0: return 1, d
    elif d < 0: return -1, d
    else: return 0, 0


def factorize(n):
    """Return dict {prime: exponent}."""
    factors = {}
    d = 2
    temp = abs(n) if n != 0 else 0
    while d * d <= temp:
        while temp % d == 0:
            factors[d] = factors.get(d, 0) + 1
            temp //= d
        d += 1
    if temp > 1:
        factors[temp] = 1
    return factors


# ═══════════════════════════════════════════════════════════════════════
# Efficient reduced Latin square enumeration for n=6
# ═══════════════════════════════════════════════════════════════════════

def enumerate_reduced_ls(n):
    """Enumerate all reduced LS of order n.
    First row = [1,2,...,n], first column = [1,2,...,n]."""
    results = []
    L = [[0]*n for _ in range(n)]
    L[0] = list(range(1, n+1))
    for i in range(n):
        L[i][0] = i + 1
    
    def fill(row, col):
        if row == n:
            results.append([r[:] for r in L])
            return
        next_row, next_col = (row, col+1) if col+1 < n else (row+1, 1)
        
        # Find available values
        used_row = set(L[row][:col])
        used_col = set(L[r][col] for r in range(row))
        available = set(range(1, n+1)) - used_row - used_col
        
        for val in sorted(available):
            L[row][col] = val
            fill(next_row, next_col)
        L[row][col] = 0
    
    fill(1, 1)
    return results


def compute_lsg_sandpile(L):
    """Compute |K(LSG)| for the Latin Square Graph of L.
    LSG: n² vertices (cells), edges when sharing row, col, or symbol."""
    n = len(L)
    N = n * n
    
    # Build adjacency matrix
    adj = [[0]*N for _ in range(N)]
    for i1 in range(n):
        for j1 in range(n):
            v1 = i1 * n + j1
            for i2 in range(n):
                for j2 in range(n):
                    if i1 == i2 and j1 == j2: continue
                    v2 = i2 * n + j2
                    if i1 == i2 or j1 == j2 or L[i1][j1] == L[i2][j2]:
                        adj[v1][v2] = 1
    
    # Laplacian = D - A
    lap = [[0]*N for _ in range(N)]
    for i in range(N):
        deg = sum(adj[i])
        for j in range(N):
            lap[i][j] = -adj[i][j]
        lap[i][i] = deg
    
    # |K(LSG)| = det of any (N-1)×(N-1) cofactor / 1
    # = det(reduced Laplacian)
    M = N - 1
    red_lap = np.array([[lap[i][j] for j in range(M)] for i in range(M)], dtype=int)
    return bareiss_det(red_lap)


# ═══════════════════════════════════════════════════════════════════════
# MAIN
# ═══════════════════════════════════════════════════════════════════════

if __name__ == "__main__":
    out = open("n6_results.txt", "w", encoding="utf-8")
    def pr(s=""):
        print(s); out.write(s + "\n"); out.flush()
    
    pr("=" * 72)
    pr("  EXHAUSTIVE n=6 ANALYSIS")
    pr("=" * 72)
    
    n = 6
    
    # ─── Step 1: Enumerate ───
    pr(f"\n>>> Enumerating reduced LS of order {n}...")
    t0 = time.time()
    all_ls = enumerate_reduced_ls(n)
    t1 = time.time()
    pr(f"    Found: {len(all_ls)} reduced LS in {t1-t0:.1f}s")
    
    if len(all_ls) != 9408:
        pr(f"    WARNING: expected 9408, got {len(all_ls)}")
    
    # ─── Step 2: Compute SNF, det, parity for all ───
    pr(f"\n>>> Computing SNF + det for all {len(all_ls)} LS...")
    t0 = time.time()
    
    snf_profiles = Counter()  # snf_tuple → count
    det_values = Counter()    # abs(det) → count
    parity_map = {}           # (parity, snf_tuple) → count
    det_to_snf = {}           # abs(det) → set of snf_tuples
    singular_count = 0
    
    all_data = []  # (idx, det_val, parity, snf_tuple)
    
    for idx, L in enumerate(all_ls):
        if idx % 1000 == 0 and idx > 0:
            elapsed = time.time() - t0
            pr(f"    ... {idx}/{len(all_ls)} ({elapsed:.1f}s)")
        
        A = get_A(L)
        A_np = np.array(A, dtype=int)
        det_val = bareiss_det(A_np)
        
        snf = snf_nonzero(A)
        parity = 1 if det_val > 0 else (-1 if det_val < 0 else 0)
        abs_det = abs(det_val)
        
        all_data.append((idx, det_val, parity, snf))
        
        snf_profiles[snf] += 1
        det_values[abs_det] += 1
        
        key = (parity, snf)
        parity_map[key] = parity_map.get(key, 0) + 1
        
        det_to_snf.setdefault(abs_det, set()).add(snf)
        
        if det_val == 0:
            singular_count += 1
    
    t1 = time.time()
    pr(f"    Done in {t1-t0:.1f}s")
    
    # ─── Results ───
    pr(f"\n{'='*72}")
    pr(f"  BLOCK 1: SNF Profile Distribution (n={n})")
    pr(f"{'='*72}")
    pr(f"  Total: {len(all_ls)} reduced LS")
    pr(f"  Singular (det=0): {singular_count}")
    pr(f"  Distinct SNF profiles: {len(snf_profiles)}")
    pr(f"  Distinct |det| values: {len(det_values)}")
    
    pr(f"\n  {'SNF':<40s} {'|K|':>10s} {'count':>8s} {'%':>6s}")
    pr(f"  {'-'*40} {'-'*10} {'-'*8} {'-'*6}")
    for snf, cnt in sorted(snf_profiles.items(), key=lambda x: (-x[1], x[0])):
        k_order = reduce(lambda a,b: a*b, snf) if snf else 1
        pct = 100 * cnt / len(all_ls)
        snf_str = " × ".join(f"Z/{d}" for d in snf) if snf else "(trivial)"
        pr(f"  {snf_str:<40s} {k_order:>10d} {cnt:>8d} {pct:>5.1f}%")
    
    pr(f"\n  |det| distribution (top 20):")
    pr(f"  {'|det|':>12s} {'count':>8s} {'factorization':<30s}")
    for dv, cnt in sorted(det_values.items(), key=lambda x: -x[1])[:20]:
        facs = factorize(dv) if dv > 0 else {}
        fac_str = " · ".join(f"{p}^{a}" if a > 1 else str(p) for p, a in sorted(facs.items())) if facs else "0"
        pr(f"  {dv:>12d} {cnt:>8d} {fac_str:<30s}")
    
    # ─── C9: Parity vs K(A) ───
    pr(f"\n{'='*72}")
    pr(f"  BLOCK 2: C9 — Parity vs K(A) (n={n})")
    pr(f"{'='*72}")
    
    # For each SNF profile, check if all instances have same parity
    snf_parities = {}  # snf_tuple → set of parities
    for idx, det_val, parity, snf in all_data:
        snf_parities.setdefault(snf, set()).add(parity)
    
    conflicts = 0
    pr(f"\n  {'SNF':<40s} {'parities':>15s} {'conflict?':>10s}")
    for snf in sorted(snf_parities.keys()):
        pars = snf_parities[snf]
        has_conflict = len(pars) > 1
        if has_conflict: conflicts += 1
        par_str = str(sorted(pars))
        tag = "YES ✗" if has_conflict else "no"
        pr(f"  {str(snf):<40s} {par_str:>15s} {tag:>10s}")
    
    pr(f"\n  Profiles with parity conflict: {conflicts}/{len(snf_parities)}")
    
    # Also check: same |det| → same parity?
    det_parities = {}
    for idx, det_val, parity, snf in all_data:
        det_parities.setdefault(abs(det_val), set()).add(parity)
    
    det_conflicts = sum(1 for v in det_parities.values() if len(v) > 1)
    pr(f"  |det| values with parity conflict: {det_conflicts}/{len(det_parities)}")
    
    # ─── C7: Divisibility by K(LSG) ───
    pr(f"\n{'='*72}")
    pr(f"  BLOCK 3: C7 — |K(A)| | |K(LSG)| (n={n})")
    pr(f"{'='*72}")
    
    # Compute K(LSG) for the cyclic LS (one representative)
    pr(f"\n  Computing |K(LSG)| for cyclic LS (n={n})...")
    cyc_ls = [[(i+j) % n + 1 for j in range(n)] for i in range(n)]
    t0 = time.time()
    K_LSG = abs(compute_lsg_sandpile(cyc_ls))
    t1 = time.time()
    pr(f"  |K(LSG)| = {K_LSG}")
    pr(f"  Computed in {t1-t0:.1f}s")
    
    # Factorize K_LSG
    facs_lsg = factorize(K_LSG)
    fac_str = " · ".join(f"{p}^{a}" if a > 1 else str(p) for p, a in sorted(facs_lsg.items()))
    pr(f"  Factorization: {fac_str}")
    prime_support_lsg = set(facs_lsg.keys())
    pr(f"  Prime support: {sorted(prime_support_lsg)}")
    
    # Check if K(LSG) is the same for a few other LS
    pr(f"\n  Checking if |K(LSG)| varies across LS:")
    sample_indices = [0, 100, 500, 1000, 5000, 9000] if len(all_ls) > 9000 else list(range(min(5, len(all_ls))))
    lsg_values = set()
    for si in sample_indices:
        if si < len(all_ls):
            kl = abs(compute_lsg_sandpile(all_ls[si]))
            lsg_values.add(kl)
            match_tag = "= same" if kl == K_LSG else f"≠ DIFFERENT ({kl})"
            pr(f"    LS {si}: |K(LSG)| {match_tag}")
    
    if len(lsg_values) == 1:
        pr(f"    → |K(LSG)| appears constant across all n={n} LS!")
    else:
        pr(f"    → |K(LSG)| varies! {len(lsg_values)} distinct values")
    
    # Now check C7 for all LS
    pr(f"\n  C7 divisibility check:")
    c7_pass = 0; c7_fail = 0
    fail_orders = Counter()
    
    for idx, det_val, parity, snf in all_data:
        k_order = abs(det_val)
        if k_order == 0:
            continue  # skip singular
        if K_LSG % k_order == 0:
            c7_pass += 1
        else:
            c7_fail += 1
            fail_orders[k_order] += 1
    
    nonsingular = c7_pass + c7_fail
    pr(f"    Non-singular LS: {nonsingular}")
    pr(f"    PASS: {c7_pass} ({100*c7_pass/nonsingular:.1f}%)")
    pr(f"    FAIL: {c7_fail} ({100*c7_fail/nonsingular:.1f}%)")
    
    if c7_fail > 0:
        pr(f"\n    Failing |K(A)| orders:")
        pr(f"    {'|K(A)|':>10s} {'count':>8s} {'factorization':<30s} {'alien primes':<20s}")
        for order, cnt in sorted(fail_orders.items()):
            facs = factorize(order)
            fac_str = " · ".join(f"{p}^{a}" if a > 1 else str(p) for p, a in sorted(facs.items()))
            aliens = set(facs.keys()) - prime_support_lsg
            alien_str = str(sorted(aliens)) if aliens else "none"
            pr(f"    {order:>10d} {cnt:>8d} {fac_str:<30s} {alien_str:<20s}")
    
    # ─── C7 refined conjecture ───
    pr(f"\n  C7 REFINED CONJECTURE: failure iff prime support ⊄ support(|K(LSG)|)")
    # Check: is it exactly the case that C7 fails iff K(A) has an "alien" prime?
    alien_fail = 0; alien_pass = 0; no_alien_fail = 0; no_alien_pass = 0
    for idx, det_val, parity, snf in all_data:
        k_order = abs(det_val)
        if k_order == 0: continue
        facs_k = factorize(k_order)
        has_alien = bool(set(facs_k.keys()) - prime_support_lsg)
        divides = (K_LSG % k_order == 0)
        
        if has_alien and not divides: alien_fail += 1
        elif has_alien and divides: alien_pass += 1
        elif not has_alien and not divides: no_alien_fail += 1
        else: no_alien_pass += 1
    
    pr(f"    alien prime + fails:  {alien_fail}")
    pr(f"    alien prime + passes: {alien_pass}")
    pr(f"    no alien + fails:     {no_alien_fail}")
    pr(f"    no alien + passes:    {no_alien_pass}")
    
    if alien_pass == 0 and no_alien_fail == 0:
        pr(f"    → CONJECTURE CONFIRMED for n={n}: "
           f"C7 fails ⟺ prime support mismatch")
    else:
        pr(f"    → CONJECTURE NEEDS REFINEMENT")
    
    # ─── Summary ───
    pr(f"\n{'='*72}")
    pr(f"  SUMMARY for n={n}")
    pr(f"{'='*72}")
    pr(f"  Reduced LS: {len(all_ls)}")
    pr(f"  Singular: {singular_count}")
    pr(f"  Distinct SNF profiles: {len(snf_profiles)}")
    pr(f"  Distinct |det|: {len(det_values)}")
    pr(f"  Parity determines K(A): {'NO' if conflicts > 0 else 'YES'}")
    pr(f"  C7 pass rate: {c7_pass}/{nonsingular} ({100*c7_pass/nonsingular:.1f}%)")
    pr(f"  C7 failure ⟺ alien prime: {'YES' if (alien_pass == 0 and no_alien_fail == 0) else 'NO'}")
    if n == 6:
        _d36 = sum(1 for _, dv, _, _ in all_data if (6*dv) % 36 == 0)
        _d18not36 = sum(1 for _, dv, _, _ in all_data if dv != 0 and (6*dv) % 18 == 0 and (6*dv) % 36 != 0)
        pr(f"  Reduced LS with 36 | det(E_std): {_d36}")
        pr(f"  Reduced LS with 18 | det(E_std) but 36 ∤ det(E_std): {_d18not36}")

    out.close()
    print(f"\n  Results saved to n6_results.txt")
