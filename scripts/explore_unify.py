"""
Three investigations in one:

1. UNIFICATION: Does dim ker(A² mod 2) ≥ 2 also characterize n²|det for n ≡ 2 mod 4?
   (n=6 and n=10)

2. WHY n=12 always has dim ker ≥ 3:
   Structural analysis of of ker_F₂(A mod 2) at n=12.

3. MINIMUM ker dim at each n ≡ 0 mod 4:
   Is there a formula? Does n=2^a have min dim ker = 1, while n=4k (k odd) has higher?
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

# ═══════════════════════════════════════════════════════════════════════
# PART 1: n ≡ 2 mod 4 — does dim ker(A²) ≥ 2 work?
# ═══════════════════════════════════════════════════════════════════════
print("=" * 72)
print("PART 1: Does dim ker(A² mod 2) ≥ 2 work for n ≡ 2 mod 4?")
print("=" * 72)

# n=6 (exhaustive)
print("\n─── n=6 EXHAUSTIVE ───")
n = 6
all_LS6 = gen_all_reduced_LS(n)
print(f"Total reduced LS: {len(all_LS6)}")

xtab = defaultdict(lambda: {'pass': 0, 'fail': 0})
xtab_old = defaultdict(lambda: {'pass': 0, 'fail': 0})
mismatch_new = 0
mismatch_old = 0

for Llist in all_LS6:
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
    
    # Old criterion: dim ker(A) ≥ 1
    old_crit = (kdim >= 1)
    # New criterion: dim ker(A²) ≥ 2
    new_crit = (kdim_A2 >= 2)
    
    if old_crit != n2_div: mismatch_old += 1
    if new_crit != n2_div: mismatch_new += 1
    
    xtab[kdim_A2]['pass' if n2_div else 'fail'] += 1
    xtab_old[kdim]['pass' if n2_div else 'fail'] += 1

print(f"\nOld criterion (dim ker(A) ≥ 1): {mismatch_old} mismatches ({'PERFECT' if mismatch_old==0 else 'FAILS'})")
print(f"New criterion (dim ker(A²) ≥ 2): {mismatch_new} mismatches ({'PERFECT' if mismatch_new==0 else 'FAILS'})")

print(f"\nCross-tab OLD: dim ker(A) vs n²|det")
print(f"  {'dim ker(A)':>12s}  {'PASS':>8s}  {'FAIL':>8s}")
for k in sorted(xtab_old):
    print(f"  {k:>12d}  {xtab_old[k]['pass']:>8d}  {xtab_old[k]['fail']:>8d}")

print(f"\nCross-tab NEW: dim ker(A²) vs n²|det")
print(f"  {'dim ker(A²)':>12s}  {'PASS':>8s}  {'FAIL':>8s}")
for k in sorted(xtab):
    print(f"  {k:>12d}  {xtab[k]['pass']:>8d}  {xtab[k]['fail']:>8d}")

# Check: is dim ker(A²) ≥ 2 EQUIVALENT to dim ker(A) ≥ 1 for n ≡ 2 mod 4?
# i.e., does ker(A) = {0} imply ker(A²) ⊆ {0,1} with dim ≤ 1?
# And ker(A) ≥ 1 implies ker(A²) ≥ 2?
print(f"\nAre the two criteria equivalent at n=6?")
equiv_count = 0
inequiv_count = 0
for Llist in all_LS6:
    L = np.array(Llist, dtype=np.int64)
    A = np.zeros((n-1, n-1), dtype=np.int64)
    for i in range(n-1):
        for j in range(n-1):
            A[i,j] = int(L[i,j]) - int(L[i,n-1])
    dA = bareiss(A.tolist())
    if n * dA == 0: continue
    kdim = kernel_dim_F2(A)
    A2 = A @ A
    kdim_A2 = kernel_dim_F2(A2)
    
    old_c = (kdim >= 1)
    new_c = (kdim_A2 >= 2)
    if old_c == new_c:
        equiv_count += 1
    else:
        inequiv_count += 1
        print(f"  DIFFER: kdim={kdim}, kdim_A2={kdim_A2}")

print(f"  Equivalent: {equiv_count}, Different: {inequiv_count}")

# ═══════════════════════════════════════════════════════════════════════
# PART 2: WHY n=12 has min ker ≥ 3
# ═══════════════════════════════════════════════════════════════════════
print()
print("=" * 72)
print("PART 2: Structural analysis of ker_F₂(A) at various n")
print("=" * 72)

# For each n ≡ 0 mod 4, generate some LS and check min ker dim
for n in [4, 8, 12, 16, 20]:
    L_cyc = np.array([[(i+j) % n + 1 for j in range(n)] for i in range(n)], dtype=np.int64)
    
    if n <= 5:
        # Exhaustive
        all_LS = gen_all_reduced_LS(n)
        kdims = []
        for Llist in all_LS:
            Larr = np.array(Llist, dtype=np.int64)
            A = np.zeros((n-1, n-1), dtype=np.int64)
            for i in range(n-1):
                for j in range(n-1):
                    A[i,j] = int(Larr[i,j]) - int(Larr[i,n-1])
            dA = bareiss(A.tolist())
            if n * dA == 0: continue
            kdims.append(kernel_dim_F2(A))
        print(f"\nn={n} (exhaustive, {len(kdims)} LS):")
        print(f"  dim ker distribution: {dict(sorted(Counter(kdims).items()))}")
        print(f"  min dim ker = {min(kdims)}, max = {max(kdims)}")
    else:
        # Sampling
        kdims = []
        num_samples = 200 if n <= 12 else 50
        for seed in range(num_samples):
            rng = np.random.default_rng(seed)
            L = swap_chain(L_cyc, rng, max(400, n*50))
            A = np.zeros((n-1, n-1), dtype=np.int64)
            for i in range(n-1):
                for j in range(n-1):
                    A[i,j] = int(L[i,j]) - int(L[i,n-1])
            dA = bareiss(A.tolist())
            if n * dA == 0: continue
            kdims.append(kernel_dim_F2(A))
        print(f"\nn={n} ({len(kdims)} samples via swap-chain):")
        print(f"  dim ker distribution: {dict(sorted(Counter(kdims).items()))}")
        print(f"  min dim ker = {min(kdims)}, max = {max(kdims)}")
        
        # Check: what is rank_F₂(A mod 2) for cyclic base?
        A_cyc = np.zeros((n-1, n-1), dtype=np.int64)
        for i in range(n-1):
            for j in range(n-1):
                A_cyc[i,j] = int(L_cyc[i,j]) - int(L_cyc[i,n-1])
        kdim_cyc = kernel_dim_F2(A_cyc)
        print(f"  Cyclic base: dim ker = {kdim_cyc}")

# ═══════════════════════════════════════════════════════════════════════
# PART 3: Theoretical analysis of the ker structure
# ═══════════════════════════════════════════════════════════════════════
print()
print("=" * 72)
print("PART 3: Structural kernel analysis")
print("=" * 72)

# For the cyclic Latin square L_{ij} = (i+j) mod n + 1:
# A_{ij} = L_{ij} - L_{i,n-1} = ((i+j) mod n) - ((i+n-1) mod n)
# = ((i+j) mod n) - ((i-1) mod n)
# Let's compute A mod 2 for cyclic LS and analyze the kernel.

for n in [4, 8, 12, 16, 20]:
    L = np.array([[(i+j) % n + 1 for j in range(n)] for i in range(n)], dtype=np.int64)
    A = np.zeros((n-1, n-1), dtype=np.int64)
    for i in range(n-1):
        for j in range(n-1):
            A[i,j] = int(L[i,j]) - int(L[i,n-1])
    
    B = A % 2  # mod 2 reduction
    kdim = kernel_dim_F2(A)
    
    # Row sums: (A·1)_i = n·((n+1)/2 - L_{i,n-1}). For cyclic: L_{i,n-1} = (i+n-1)%n+1.
    row_sums = A.sum(axis=1)
    row_sums_mod2 = row_sums % 2
    row_sums_mod4 = row_sums % 4
    
    print(f"\nn={n}: dim ker_F₂ of cyclic = {kdim}")
    print(f"  A mod 2 has rank {n-1-kdim} out of {n-1}")
    print(f"  Row sums mod 2: {list(row_sums_mod2[:6])}...")
    print(f"  Row sums mod 4: {list(row_sums_mod4[:6])}...")
    
    # The key structural question: for n=12, why is rank so low?
    # A_{ij} = ((i+j) mod 12) - ((i+11) mod 12) = ((i+j) mod 12) - ((i+11) mod 12)
    # All entries of A mod 2: let's look at the pattern
    if n <= 16:
        print(f"  B = A mod 2 ({n-1}×{n-1}):")
        for i in range(min(n-1, 12)):
            print(f"    {list(B[i,:min(n-1,12)])}")

# Summary: Is there a pattern for min ker dim as function of n?
print()
print("=" * 72)
print("Q: Minimum possible dim ker_F₂ for n ≡ 0 mod 4?")
print("=" * 72)
print("  n=4:  min ker = 1 (from exhaustive)")
print("  n=8:  min ker = 1 (from 4500+ samples)")
print("  n=12: min ker = 3 (from 500+ samples)")
print("  n=16: min ker = 1 (from 200 samples)")
print("  n=20: see above")
print()
print("Hypothesis: min dim ker depends on arithmetic structure of n.")
print("  n=4=2²: min=1")
print("  n=8=2³: min=1")
print("  n=12=4·3: min=3 (extra ker vectors from 3|n?)")
print("  n=16=2⁴: min=1") 
print("  n=20=4·5: prediction → min≥5?")
