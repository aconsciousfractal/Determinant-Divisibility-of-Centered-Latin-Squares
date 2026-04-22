"""
BREAKTHROUGH VERIFICATION: dim ker_F₂(A² mod 2) ≥ 2 as the UNIFIED criterion.

Key insight from explore_dim1_deep.py:
  For dim ker_F₂(A) = 1:
    PASS (n²|det) ⟺ dim ker_F₂(A²) = 2
    FAIL (n²∤det) ⟺ dim ker_F₂(A²) = 1

Since ker(A) ⊆ ker(A²), for dim ker(A) ≥ 2 we automatically get dim ker(A²) ≥ 2.

UNIFIED CONJECTURE: n² | det(E_std) ⟺ dim ker_F₂(A² mod 2) ≥ 2

Equivalent reformulation:
  ker(A²) has dim ≥ 2 iff ∃ w with A·w ≡ 1 (mod 2) [i.e. 1 ∈ im(A mod 2)]
  OR dim ker(A) ≥ 2.

So: n² | det ⟺ 1 ∈ im(A mod 2) OR rank_F₂(A) ≤ n-3.

This script verifies the UNIFIED criterion at n=4 (exhaustive) and n=8 (massive sample).
"""
import sys as _sys
try:
    _sys.stdout.reconfigure(encoding="utf-8")
    _sys.stderr.reconfigure(encoding="utf-8")
except Exception:
    pass
import numpy as np
from collections import Counter, defaultdict

def bareiss(A):
    n = len(A); B = [[int(x) for x in r] for r in A]; sign = 1; prev = 1
    for i in range(n):
        piv = None
        for k in range(i, n):
            if B[k][i] != 0: piv = k; break
        if piv is None: return 0
        if piv != i: B[i], B[piv] = B[piv], B[i]; sign *= -1
        for j in range(i+1, n):
            for k2 in range(i+1, n):
                B[j][k2] = (B[i][i]*B[j][k2] - B[j][i]*B[i][k2]) // prev
            B[j][i] = 0
        prev = B[i][i]
    return sign * B[n-1][n-1]

def v_p(x, p):
    if x == 0: return 999
    v = 0; t = abs(x)
    while t % p == 0: t //= p; v += 1
    return v

def kernel_dim_F2(M):
    n_rows, n_cols = M.shape
    B = M.copy() % 2
    pivot_cols = []; row = 0
    for col in range(n_cols):
        found = None
        for r in range(row, n_rows):
            if B[r, col] % 2 == 1: found = r; break
        if found is None: continue
        B[[row, found]] = B[[found, row]]
        pivot_cols.append(col)
        for r in range(n_rows):
            if r != row and B[r, col] % 2 == 1:
                B[r] = (B[r] + B[row]) % 2
        row += 1
    return n_cols - len(pivot_cols)

def one_in_image_F2(M):
    """Check if 1-vector is in the column space of M mod 2."""
    n_rows, n_cols = M.shape
    # Augmented system: [M | 1] mod 2, check consistency
    B = np.zeros((n_rows, n_cols + 1), dtype=np.int64)
    B[:, :n_cols] = M % 2
    B[:, n_cols] = 1  # target = all-ones
    
    # Row reduce
    row = 0
    for col in range(n_cols):
        found = None
        for r in range(row, n_rows):
            if B[r, col] % 2 == 1: found = r; break
        if found is None: continue
        B[[row, found]] = B[[found, row]]
        for r in range(n_rows):
            if r != row and B[r, col] % 2 == 1:
                B[r] = (B[r] + B[row]) % 2
        row += 1
    
    # Check if any row below pivot rows has [0...0 | 1]
    for r in range(row, n_rows):
        if B[r, n_cols] % 2 == 1:
            return False  # inconsistent
    return True

def swap_chain(L, rng, steps=500):
    n = L.shape[0]; L2 = L.copy()
    for _ in range(steps):
        r = rng.integers(n); c1, c2 = rng.choice(n, 2, replace=False)
        v1, v2_ = L2[r,c1], L2[r,c2]
        cs = [r2 for r2 in range(n) if r2 != r and L2[r2,c1] == v2_ and L2[r2,c2] == v1]
        if cs:
            r2 = rng.choice(cs)
            L2[r,c1], L2[r,c2] = v2_, v1
            L2[r2,c1], L2[r2,c2] = v1, v2_
    return L2

# ═══════════════════════════════════════════════════════════════════════
# EXHAUSTIVE n=4
# ═══════════════════════════════════════════════════════════════════════
def gen_all_reduced_LS(n):
    results = []
    def fill(L, row, col):
        if row == n:
            results.append([r[:] for r in L])
            return
        next_row, next_col = (row, col+1) if col+1 < n else (row+1, 1)
        used_row = set(L[row][:col])
        used_col = set(L[r][col] for r in range(row))
        for v in range(1, n+1):
            if v not in used_row and v not in used_col:
                L[row][col] = v
                fill(L, next_row, next_col)
                L[row][col] = 0
    L = [[0]*n for _ in range(n)]
    for j in range(n): L[0][j] = j+1
    for i in range(n): L[i][0] = i+1
    fill(L, 1, 1)
    return results

print("=" * 72)
print("UNIFIED CRITERION VERIFICATION")
print("Conjecture: n² | det(E) ⟺ dim ker_F₂(A² mod 2) ≥ 2")
print("           ⟺ 1 ∈ im(A mod 2)  OR  dim ker_F₂(A) ≥ 2")
print("=" * 72)

# ── n=4 exhaustive ────────────────────────────────────────────────────
print("\n─── n=4 (EXHAUSTIVE) ───")
n = 4
all_LS = gen_all_reduced_LS(n)
results_n4 = []
for Llist in all_LS:
    L = np.array(Llist, dtype=np.int64)
    A = np.zeros((n-1, n-1), dtype=np.int64)
    for i in range(n-1):
        for j in range(n-1):
            A[i,j] = int(L[i,j]) - int(L[i,n-1])
    dA = bareiss(A.tolist())
    dE = n * dA
    if dE == 0: continue
    
    n2_div = (dE % (n**2) == 0)
    kdim = kernel_dim_F2(A)
    A2 = A @ A
    kdim_A2 = kernel_dim_F2(A2)
    one_im = one_in_image_F2(A)
    
    print(f"  LS: det(A)={dA:8d}, v₂={v_p(dA,2)}, n²|det={'Y' if n2_div else 'N'}, "
          f"dim ker(A)={kdim}, dim ker(A²)={kdim_A2}, 1∈im(A)={one_im}")
    
    # Verify criterion
    criterion = (kdim_A2 >= 2)
    criterion_alt = (kdim >= 2 or one_im)
    assert criterion == n2_div, f"CRITERION FAILS at n=4!"
    assert criterion == criterion_alt, f"Two formulations disagree!"

print(f"  n=4: CRITERION VERIFIED on all {len(all_LS)} reduced LS ✓")

# ── n=8 massive sample ───────────────────────────────────────────────
print("\n─── n=8 (4500 samples, 3 bases) ───")
n = 8
L_cyc = np.array([[(i+j) % n + 1 for j in range(n)] for i in range(n)], dtype=np.int64)
L_nc8 = np.array([
    [1,2,3,4,5,6,7,8],[2,1,4,3,6,5,8,7],[3,4,1,2,7,8,5,6],
    [4,3,2,1,8,7,6,5],[5,6,7,8,1,2,3,4],[6,5,8,7,2,1,4,3],
    [7,8,5,6,3,4,1,2],[8,7,6,5,4,3,2,1]], dtype=np.int64)
L_d4 = np.array([
    [1,2,3,4,5,6,7,8],[2,3,4,1,6,7,8,5],[3,4,1,2,7,8,5,6],
    [4,1,2,3,8,5,6,7],[5,8,7,6,1,4,3,2],[6,5,8,7,2,1,4,3],
    [7,6,5,8,3,2,1,4],[8,7,6,5,4,3,2,1]], dtype=np.int64)
bases = [L_cyc, L_nc8, L_d4]

total = 0
pass_count = 0
fail_count = 0
criterion_match = 0
criterion_mismatch = 0
alt_match = 0

# Cross-tabulation for the unified criterion
xtab = defaultdict(lambda: {'pass': 0, 'fail': 0})

for bi, Lb in enumerate(bases):
    for seed in range(1500):
        rng = np.random.default_rng(seed + bi * 100000)
        L = swap_chain(Lb, rng, 1000)
        
        A = np.zeros((n-1, n-1), dtype=np.int64)
        for i in range(n-1):
            for j in range(n-1):
                A[i,j] = int(L[i,j]) - int(L[i,n-1])
        
        dA = bareiss(A.tolist())
        dE = n * dA
        if dE == 0: continue
        total += 1
        
        n2_div = (dE % (n**2) == 0)
        
        # Unified criterion
        A2 = A @ A
        kdim_A2 = kernel_dim_F2(A2)
        criterion = (kdim_A2 >= 2)
        
        # Alternative formulation
        kdim = kernel_dim_F2(A)
        one_im = one_in_image_F2(A)
        criterion_alt = (kdim >= 2 or one_im)
        
        if criterion == n2_div:
            criterion_match += 1
        else:
            criterion_mismatch += 1
        
        if criterion_alt == n2_div:
            alt_match += 1
        
        if n2_div:
            pass_count += 1
        else:
            fail_count += 1
        
        xtab[kdim_A2]['pass' if n2_div else 'fail'] += 1

print(f"  Total: {total}")
print(f"  n²|det: {pass_count} PASS, {fail_count} FAIL")
print()
print(f"  UNIFIED CRITERION (dim ker(A²) ≥ 2):")
print(f"    Match: {criterion_match}/{total} ({100*criterion_match/total:.2f}%)")
print(f"    Mismatch: {criterion_mismatch}/{total}")
if criterion_mismatch == 0:
    print(f"    ★★★ PERFECT — ZERO mismatches! ★★★")
print()
print(f"  ALT FORMULATION (dim ker(A) ≥ 2 OR 1∈im(A)):")
print(f"    Match: {alt_match}/{total} ({100*alt_match/total:.2f}%)")
if alt_match == total:
    print(f"    ★★★ PERFECT — ZERO mismatches! ★★★")

print()
print(f"  Cross-tabulation: dim ker(A²) vs n²-divisibility")
print(f"  {'dim ker(A²)':>12s}  {'PASS':>8s}  {'FAIL':>8s}")
print(f"  {'─'*12}  {'─'*8}  {'─'*8}")
for k in sorted(xtab):
    print(f"  {k:>12d}  {xtab[k]['pass']:>8d}  {xtab[k]['fail']:>8d}")

# ── n=12 sample ──────────────────────────────────────────────────────
print("\n─── n=12 (300 samples) ───")
n = 12
L_cyc12 = np.array([[(i+j) % n + 1 for j in range(n)] for i in range(n)], dtype=np.int64)

# Z4 x Z3
L_z4z3 = np.zeros((n, n), dtype=np.int64)
for i in range(n):
    for j in range(n):
        i1, i2 = i % 4, i // 4
        j1, j2 = j % 4, j // 4
        s1, s2 = (i1+j1) % 4, (i2+j2) % 3
        L_z4z3[i,j] = s1 + s2*4 + 1

bases12 = [L_cyc12, L_z4z3]
total12 = 0; match12 = 0; mismatch12 = 0

for bi, Lb in enumerate(bases12):
    for seed in range(150):
        rng = np.random.default_rng(seed + bi * 100000)
        L = swap_chain(Lb, rng, 600)
        
        A = np.zeros((n-1, n-1), dtype=np.int64)
        for i in range(n-1):
            for j in range(n-1):
                A[i,j] = int(L[i,j]) - int(L[i,n-1])
        
        dA = bareiss(A.tolist())
        dE = n * dA
        if dE == 0: continue
        total12 += 1
        
        n2_div = (dE % (n**2) == 0)
        A2 = A @ A
        kdim_A2 = kernel_dim_F2(A2)
        criterion = (kdim_A2 >= 2)
        
        if criterion == n2_div:
            match12 += 1
        else:
            mismatch12 += 1

print(f"  Total: {total12}")
print(f"  Criterion match: {match12}/{total12}")
print(f"  Mismatch: {mismatch12}")
if mismatch12 == 0:
    print(f"  ★★★ PERFECT at n=12 too ★★★")

# ── FINAL SUMMARY ────────────────────────────────────────────────────
print()
print("=" * 72)
print("FINAL SUMMARY")
print("=" * 72)
total_all = len(all_LS) + total + total12
match_all = len(all_LS) + criterion_match + match12
print(f"""
THEOREM CANDIDATE (necessary AND sufficient for n ≡ 0 mod 4):

  n² | det(E_std) ⟺ dim ker_F₂((A mod 2)²) ≥ 2

Equivalently:

  n² | det(E_std) ⟺ 1 ∈ im(A mod 2)  OR  rank_F₂(A mod 2) ≤ n-3

EVIDENCE:
  n=4:  exhaustive ({len(all_LS)} reduced LS)        — 0 mismatches
  n=8:  {total} swap-chain samples (3 bases)  — {criterion_mismatch} mismatches
  n=12: {total12} swap-chain samples (2 bases)  — {mismatch12} mismatches
  TOTAL: {match_all}/{total_all} perfect matches

MATHEMATICAL EXPLANATION:
  Write B = A mod 2 (the F₂ reduction).
  - 1 ∈ ker(B) always (for 4|n).
  - If dim ker(B) ≥ 2: at least 2 elementary divisors of A are even
    → v₂(det A) ≥ a+1 where n=2^a·q → n²|det(E).
  - If dim ker(B) = 1 (only 1): need 1 ∈ im(B), meaning ∃ w: Bw = 1.
    This ensures the "lift" of the kernel vector 1 to Z/4Z yields
    a second even elementary divisor in the Smith NF.
    → v₂(det A) ≥ 3 → n²|det(E).
  - If dim ker(B) = 1 AND 1 ∉ im(B): only one even elementary divisor,
    and it contributes only v₂ = 1 to det(A).
    → v₂(det A) = 2 exactly → n² ∤ det(E).
""")
