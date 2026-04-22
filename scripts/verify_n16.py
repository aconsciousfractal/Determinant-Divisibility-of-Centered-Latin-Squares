"""
Verify unified criterion at n=16 (n = 2^4, a=4).
n²=256, v₂(n²)=8, v₂(n)=4.
Need v₂(det A) ≥ 4 for n²|det E = n·det A.

This is the HARDEST test because a=4 (highest 2-power seen so far).
At n=8 (a=3): the criterion was dim ker_F₂(A²) ≥ 2.
Question: does the SAME criterion work at n=16, or do we need more
(e.g., dim ker_F₂(A²) ≥ a-1 = 3)?

Also: is dim ker(A) always ≥ 1 at n=16? (Yes: 1∈ker since 16/2=8 is even.)
"""
import numpy as np
from collections import Counter, defaultdict
import time

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
    n_rows, n_cols = M.shape
    B = np.zeros((n_rows, n_cols + 1), dtype=np.int64)
    B[:, :n_cols] = M % 2
    B[:, n_cols] = 1
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
    for r in range(row, n_rows):
        if B[r, n_cols] % 2 == 1:
            return False
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

n = 16
print("=" * 72)
print(f"UNIFIED CRITERION at n={n} (n = 2^{v_p(n,2)}, a={v_p(n,2)})")
print("=" * 72)
print(f"n²={n**2}, v₂(n²)={v_p(n**2,2)}, v₂(n)={v_p(n,2)}")
print(f"Need v₂(det A) ≥ {v_p(n**2,2) - v_p(n,2)} for n²|det E")
print()

# Base LS
L_cyc = np.array([[(i+j) % n + 1 for j in range(n)] for i in range(n)], dtype=np.int64)

# Z8 x Z2
L_z8z2 = np.zeros((n,n), dtype=np.int64) 
for i in range(n):
    for j in range(n):
        i1, i2 = i % 8, i // 8
        j1, j2 = j % 8, j // 8
        L_z8z2[i,j] = (i1+j1)%8 + ((i2+j2)%2)*8 + 1

# Z4 x Z4
L_z4z4 = np.zeros((n,n), dtype=np.int64)
for i in range(n):
    for j in range(n):
        i1, i2 = i % 4, i // 4
        j1, j2 = j % 4, j // 4
        L_z4z4[i,j] = (i1+j1)%4 + ((i2+j2)%4)*4 + 1

# Z4 x Z2 x Z2
L_z4z2z2 = np.zeros((n,n), dtype=np.int64)
for i in range(n):
    for j in range(n):
        i1, i2, i3 = i%4, (i//4)%2, i//8
        j1, j2, j3 = j%4, (j//4)%2, j//8
        L_z4z2z2[i,j] = (i1+j1)%4 + ((i2+j2)%2)*4 + ((i3+j3)%2)*8 + 1

bases = [L_cyc]
for name, Lb in [("Z8×Z2", L_z8z2), ("Z4×Z4", L_z4z4), ("Z4×Z2×Z2", L_z4z2z2)]:
    ok = (all(len(set(Lb[i,:])) == n for i in range(n)) and 
          all(len(set(Lb[:,j])) == n for j in range(n)))
    if ok:
        bases.append(Lb)
        print(f"  Base '{name}' validated ✓")
    else:
        print(f"  Base '{name}' INVALID ✗")

print(f"Using {len(bases)} base LS")
print()

NUM_PER_BASE = 50  # n=16 is SLOW (15x15 bareiss)
SWAP_STEPS = 400

total = 0
match = 0
mismatch = 0
v2_dist = Counter()
kdim_dist = Counter()
kdim_A2_dist = Counter()
xtab = defaultdict(lambda: {'pass': 0, 'fail': 0})
n2_pass = 0
n2_fail = 0

t0 = time.time()

for bi, Lb in enumerate(bases):
    for seed in range(NUM_PER_BASE):
        rng = np.random.default_rng(seed + bi * 100000)
        L = swap_chain(Lb, rng, SWAP_STEPS)
        
        A = np.zeros((n-1, n-1), dtype=np.int64)
        for i in range(n-1):
            for j in range(n-1):
                A[i,j] = int(L[i,j]) - int(L[i,n-1])
        
        dA = bareiss(A.tolist())
        dE = n * dA
        if dE == 0: continue
        total += 1
        
        v2_dE = v_p(dE, 2)
        v2_dA = v_p(dA, 2)
        v2_dist[v2_dE] += 1
        
        n2_div = (dE % (n**2) == 0)
        if n2_div: n2_pass += 1
        else: n2_fail += 1
        
        kdim = kernel_dim_F2(A)
        A2 = A @ A
        kdim_A2 = kernel_dim_F2(A2)
        
        kdim_dist[kdim] += 1
        kdim_A2_dist[kdim_A2] += 1
        
        criterion = (kdim_A2 >= 2)
        if criterion == n2_div:
            match += 1
        else:
            mismatch += 1
            print(f"  *** MISMATCH: seed={seed}, bi={bi}, v₂(dA)={v2_dA}, "
                  f"kdim={kdim}, kdim_A2={kdim_A2}, n²|det={n2_div}")
        
        xtab[kdim_A2]['pass' if n2_div else 'fail'] += 1
        
        if total % 20 == 0:
            elapsed = time.time() - t0
            print(f"  ... {total} done ({elapsed:.0f}s)", flush=True)

elapsed = time.time() - t0
print(f"\nTotal: {total} in {elapsed:.1f}s")
print(f"n²|det: {n2_pass} PASS, {n2_fail} FAIL ({100*n2_fail/total:.1f}% fail)")

print(f"\nv₂(det E) distribution:")
for k in sorted(v2_dist):
    pct = 100 * v2_dist[k] / total
    marker = f" *** < v₂(n²)={v_p(n**2,2)}" if k < v_p(n**2,2) else ""
    print(f"  v₂ = {k:3d}: {v2_dist[k]:5d}  ({pct:5.1f}%){marker}")

print(f"\ndim ker_F₂(A):")
for k in sorted(kdim_dist):
    print(f"  dim ker = {k}: {kdim_dist[k]}")

print(f"\nUNIFIED CRITERION (dim ker(A²) ≥ 2):")
print(f"  Match: {match}/{total}")
print(f"  Mismatch: {mismatch}")
if mismatch == 0:
    print(f"  ★★★ PERFECT at n={n} ★★★")
else:
    print(f"  ✗ CRITERION FAILS at n={n}")
    print(f"  Need deeper analysis (perhaps dim ker(A²) ≥ a-1?)")

print(f"\nCross-tab: dim ker(A²) vs n²-divisibility")
print(f"  {'dim ker(A²)':>12s}  {'PASS':>8s}  {'FAIL':>8s}")
print(f"  {'─'*12}  {'─'*8}  {'─'*8}")
for k in sorted(xtab):
    print(f"  {k:>12d}  {xtab[k]['pass']:>8d}  {xtab[k]['fail']:>8d}")
