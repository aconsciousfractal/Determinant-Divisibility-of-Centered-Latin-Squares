"""
Exploration at n=12 (n ≡ 0 mod 4, a=2 since 12=4·3).
Key questions:
  - Does dim ker_F₂ ≥ 2 ⟹ n²|det(E) hold?
  - What is the failure rate?
  - Is the sharp bound v₂(det A) ≥ dim_ker + 1 still valid?
  - For dim_ker=1: is the mod-8 separation still perfect?
  
n=12: n²=144, v₂(n²)=4, v₂(n)=2, need v₂(det A) ≥ 2 for n²|det E.
This is WEAKER than n=8 (needs ≥3). So dim ker ≥ 1 (= just **1** in ker) 
already gives v₂≥2 via conj 29.2 → ALL should pass!

Wait: v₂(12²)=v₂(144)=4, v₂(12)=v₂(12)=2, so v₂(det E)=v₂(n)+v₂(det A).
Need v₂(det E)≥4, so v₂(det A)≥2. 
From row sums: (A·1)_i = (n/2)·(odd), n/2=6, v₂(6)=1.
So 1 ∈ ker(A mod 2) gives v₂(det A) ≥ 1 extra. Combined with base floor:
  base: from row sums each row sum = 6·odd → v₂ contribution...

Let me just compute and see.
"""
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

# ═══════════════════════════════════════════════════════════════════════
n = 12
print("=" * 72)
print(f"EXPLORATION: n={n} (n ≡ {n%4} mod 4)")
print("=" * 72)
print(f"n²={n**2}, v₂(n²)={v_p(n**2,2)}, v₂(n)={v_p(n,2)}")
print(f"Need v₂(det A) ≥ {v_p(n**2,2) - v_p(n,2)} for n²|det E")
print(f"n/2 = {n//2}, v₂(n/2) = {v_p(n//2,2)}")
print(f"Since n/2={n//2} is even (4|n), 1 ∈ ker(A mod 2) automatically.")
print()

# Base Latin squares
L_cyc = np.array([[(i+j) % n + 1 for j in range(n)] for i in range(n)], dtype=np.int64)

# Z6 x Z2 (abelian group of order 12)
L_z6z2 = np.zeros((n, n), dtype=np.int64)
for i in range(n):
    for j in range(n):
        # (i1,i2) + (j1,j2) in Z6 x Z2
        i1, i2 = i % 6, i // 6
        j1, j2 = j % 6, j // 6
        s1, s2 = (i1+j1) % 6, (i2+j2) % 2
        L_z6z2[i,j] = s1 + s2*6 + 1

# Z4 x Z3 (another abelian group)
L_z4z3 = np.zeros((n, n), dtype=np.int64)
for i in range(n):
    for j in range(n):
        i1, i2 = i % 4, i // 4
        j1, j2 = j % 4, j // 4
        s1, s2 = (i1+j1) % 4, (i2+j2) % 3
        L_z4z3[i,j] = s1 + s2*4 + 1

bases = [L_cyc]
for name, Lb in [("Z6×Z2", L_z6z2), ("Z4×Z3", L_z4z3)]:
    ok = (all(len(set(Lb[i,:])) == n for i in range(n)) and 
          all(len(set(Lb[:,j])) == n for j in range(n)))
    if ok:
        bases.append(Lb)
        print(f"  Base '{name}' validated ✓")
    else:
        print(f"  Base '{name}' INVALID ✗")

print(f"Using {len(bases)} base LS")
print()

# Sampling
NUM_PER_BASE = 200  # n=12 is slower
SWAP_STEPS = 600

tab = defaultdict(lambda: {'pass': 0, 'fail': 0, 'v2_min': 999, 'v2s': []})
total = 0
v2_dist = Counter()
dim1_mod = Counter()  # det(A) mod 4 for dim=1

import time
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
        v2_dA = v_p(dA, 2)
        v2_dE = v_p(dE, 2)
        kdim = kernel_dim_F2(A)
        n2_div = (dE % (n**2) == 0)
        
        v2_dist[v2_dE] += 1
        tab[kdim]['v2_min'] = min(tab[kdim]['v2_min'], v2_dA)
        tab[kdim]['v2s'].append(v2_dA)
        if n2_div:
            tab[kdim]['pass'] += 1
        else:
            tab[kdim]['fail'] += 1
        
        if kdim == 1:
            dim1_mod[(abs(dA) % 4, n2_div)] += 1

elapsed = time.time() - t0
print(f"Sampled {total} LS in {elapsed:.1f}s")

# ── Results ────────────────────────────────────────────────────────────
print()
print("─" * 50)
print("v₂(det E_std) distribution")
print("─" * 50)
for k in sorted(v2_dist):
    pct = 100 * v2_dist[k] / total
    marker = f" *** < v₂(n²)={v_p(n**2,2)}" if k < v_p(n**2,2) else ""
    print(f"  v₂ = {k:3d}: {v2_dist[k]:5d}  ({pct:5.1f}%){marker}")

n2_pass = sum(tab[k]['pass'] for k in tab)
n2_fail = sum(tab[k]['fail'] for k in tab)
print(f"\nn² | det E: {n2_pass}/{total} ({100*n2_pass/total:.1f}%)")
print(f"n² ∤ det E: {n2_fail}/{total} ({100*n2_fail/total:.1f}%)")

print()
print("─" * 50)
print("CROSS-TABULATION: dim ker_F₂ vs n²-divisibility")
print("─" * 50)
print(f"  {'dim ker':>8s}  {'PASS':>8s}  {'FAIL':>8s}  {'min v₂(dA)':>12s}  {'Total':>8s}")
print(f"  {'─'*8}  {'─'*8}  {'─'*8}  {'─'*12}  {'─'*8}")
for k in sorted(tab):
    d = tab[k]
    t = d['pass'] + d['fail']
    print(f"  {k:>8d}  {d['pass']:>8d}  {d['fail']:>8d}  {d['v2_min']:>12d}  {t:>8d}")

all_pass_ge2 = all(tab[k]['fail'] == 0 for k in tab if k >= 2)
print(f"\n  Conj 29.1 (dim ker ≥ 2 ⟹ n²|det)? {'YES ✓' if all_pass_ge2 else 'NO ✗'}")

# Sharp bound
print()
print("─" * 50)
print("Sharp bound: v₂(det A) ≥ dim_ker + 1 ?")
print("─" * 50)
for k in sorted(tab):
    d = tab[k]
    bound = k + 1
    status = "✓" if d['v2_min'] >= bound else f"✗ (min={d['v2_min']} < {bound})"
    print(f"  dim ker = {k}: min v₂(dA) = {d['v2_min']} ≥ {bound}?  {status}")

# n=12 specific: v₂(det A) ≥ 2 needed. Since dim ker ≥ 1 always (1∈ker),
# Conj 29.2 gives v₂ ≥ 2 → should ALL pass n²|det!
need = v_p(n**2, 2) - v_p(n, 2)
print(f"\n  For n={n}: need v₂(det A) ≥ {need}.")
print(f"  Since dim ker ≥ 1 always, Conj 29.2 gives v₂ ≥ 2.")
if need <= 2:
    print(f"  → Conj 29.2 predicts ALL n={n} LS should have n²|det!")
    if n2_fail == 0:
        print(f"  ✓ CONFIRMED: 0 failures in {total} samples.")
    else:
        print(f"  ✗ DISPROVED: {n2_fail} failures found!")

# dim=1 mod analysis
if dim1_mod:
    print()
    print("─" * 50)
    print("Dim ker=1: det(A) mod 4")  
    print("─" * 50)
    for (mod4, n2d), cnt in sorted(dim1_mod.items()):
        status = "PASS" if n2d else "FAIL"
        print(f"  |det(A)| mod 4 = {mod4}, {status}: {cnt}")

print()
print("=" * 72)
print("v₂(det A) full distribution by ker dim")
print("=" * 72)
for kdim in sorted(tab):
    sub = tab[kdim]['v2s']
    if not sub: continue
    c = Counter(sub)
    print(f"\n  dim ker = {kdim} ({len(sub)} samples):")
    for v in sorted(c):
        print(f"    v₂(det A) = {v}: {c[v]}  ({100*c[v]/len(sub):.1f}%)")
