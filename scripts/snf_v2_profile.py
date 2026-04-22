"""
snf_v2_profile.py — Analyze 2-adic invariant factor structure via rank mod 2^k.

Key insight from n=4 exhaustive (SNF computed):
  - kdim=1: SNF = [1, 1, 8] for ALL 144 LS → v2 profile (0,0,3)
  - kdim=2: SNF = [1, 4, 4] for ALL 144 LS → v2 profile (0,2,2)

The number of invariant factors with v2 >= t equals:
  (n-1) - rank_Z/2^t(A mod 2^t)

We compute rank_Z/2^t using the "hopping ranks" trick:
  r_t = #{elementary divisors NOT divisible by 2^t}
  This equals the rank of A mod 2^t over Z/2^t.

For practical computation, we use the fact that for a matrix M over Z/2^t:
  rank(M mod 2^t) = r  <=>  exactly r elementary divisors are units mod 2^t
  i.e., the number of invariant factors with v2 < t.

So: #{d_i : v2(d_i) >= t} = (n-1) - rank_Z/2^t(A mod 2^t)
And: #{d_i : v2(d_i) = t} = rank_Z/2^{t+1}(A) - rank_Z/2^t(A) ... NO, reversed:
  = (n-1 - rank_t) - (n-1 - rank_{t+1}) = rank_{t+1} - rank_t  if rank is "number of units"

Actually, the cleanest approach: compute rank(A mod 2) over F_2 (= number of NONZERO
invariant factors mod 2), then rank(A mod 4) over Z/4Z, etc.

But Z/2^k isn't a field, so "rank" needs care. Instead, use:
  #{i : 2^t | d_i} = dim ker(A mod 2^t over Z/2^t) ... not quite.

SIMPLEST: Use the formula
  #{i : v2(d_i) >= t} = #{i : d_i ≡ 0 mod 2^t} = (n-1) - rank_2^t(A)
  where rank_2^t(A) = F_2-rank of (A / gcd-extracted mod 2)... 

Actually the standard result is:
  For integer matrix A with SNF d_1|...|d_m,
  #{i : p^k | d_i} = m - rank_{F_p}(A mod p^k / p^{k-1})... no.

The classical way: for p prime, the p-adic profile of the SNF is determined by:
  rank(A mod p) = #{d_i : p ∤ d_i}  (i.e., # units mod p)
  
Then from the lifted matrix approach:
  A = 2^0 * (odd part) stacked.

SIMPLEST APPROACH: For the LAST invariant factor d_{n-1}:
  v2(d_{n-1}) = v2(det A) - sum of v2 of other factors
  
But let's not overthink this. Let me just compute:
  - det(A) and v2(det(A))
  - kdim = dim ker_F2(A mod 2) = number of elementary divisors divisible by 2
  - rank of A mod 4 over Z/4Z (using SNF-lite)
  - rank of A mod 8 over Z/8Z
And from these, reconstruct the 2-adic profile of the SNF.

The number of elementary divisors divisible by 2^t is:
  f(t) = (n-1) - r_t
where r_t is the "2^t-rank" of A.

r_0 = n-1 (all factors are nonzero if det != 0)
r_1 = (n-1) - kdim = rank_F2(A mod 2)
r_t for t >= 2: hardest to compute without full SNF.

ALTERNATIVE: Just compute the SNF but with a BETTER algorithm.
Use the Kannan-Bachem algorithm or just use Python's exact arithmetic
with careful GCD operations. My previous algorithm might have been stuck
in a loop; let me add a max iterations guard.

Actually, let me just implement a proper modular approach and validate
against the known n=4 data.
"""
import sys as _sys
try:
    _sys.stdout.reconfigure(encoding="utf-8")
    _sys.stderr.reconfigure(encoding="utf-8")
except Exception:
    pass

import numpy as np
import random, time, sys
from collections import Counter, defaultdict
from math import gcd

# ── Helpers (reused) ──────────────────────────────────────────────

def kernel_dim_F2(M_int):
    B = np.array(M_int % 2, dtype=int)
    rows, cols = B.shape; pr = 0
    for col in range(cols):
        found = False
        for row in range(pr, rows):
            if B[row, col] == 1:
                B[[pr, row]] = B[[row, pr]]; found = True; break
        if not found: continue
        for row in range(rows):
            if row != pr and B[row, col] == 1:
                B[row] = (B[row] + B[pr]) % 2
        pr += 1
    return cols - pr

def bareiss(M):
    n = M.shape[0]; A = M.astype(object); sign = 1; prev = 1
    for i in range(n):
        if A[i, i] == 0:
            found = False
            for k in range(i + 1, n):
                if A[k, i] != 0:
                    A[[i, k]] = A[[k, i]]; sign *= -1; found = True; break
            if not found: return 0
        for j in range(i + 1, n):
            for k in range(i + 1, n):
                A[j, k] = (A[j, k] * A[i, i] - A[j, i] * A[i, k]) // prev
            A[j, i] = 0
        prev = A[i, i]
    return sign * A[n - 1, n - 1]

def v2(x):
    if x == 0: return 999
    x = abs(int(x)); v = 0
    while x % 2 == 0: x //= 2; v += 1
    return v

def get_A(L):
    n = L.shape[0]
    return np.array([[L[i,j]-L[i,n-1] for j in range(n-1)] for i in range(n-1)], dtype=int)

def jm_sample(n, steps):
    cube = [[[0]*n for _ in range(n)] for _ in range(n)]
    for i in range(n):
        for j in range(n):
            cube[i][j][(i+j)%n] = 1
    proper = True; ir = ic = ik = 0; count = 0
    while count < steps or not proper:
        if proper:
            i = random.randint(0,n-1); j = random.randint(0,n-1)
            row_ij = cube[i][j]; k_cur = 0
            while row_ij[k_cur] != 1: k_cur += 1
            k = random.randint(0,n-2)
            if k >= k_cur: k += 1
            count += 1
        else:
            i,j,k = ir,ic,ik
            cands = [x for x in range(n) if cube[i][j][x]==1]
            k_cur = cands[random.randint(0,len(cands)-1)]
        jc = [jj for jj in range(n) if jj!=j and cube[i][jj][k]==1]
        j_p = jc[random.randint(0,len(jc)-1)]
        ic_ = [ii for ii in range(n) if ii!=i and cube[ii][j][k]==1]
        i_p = ic_[random.randint(0,len(ic_)-1)]
        cube[i][j][k]+=1; cube[i][j][k_cur]-=1
        cube[i][j_p][k]-=1; cube[i][j_p][k_cur]+=1
        cube[i_p][j][k]-=1; cube[i_p][j][k_cur]+=1
        cube[i_p][j_p][k]+=1; cube[i_p][j_p][k_cur]-=1
        if cube[i_p][j_p][k_cur]==-1:
            proper=False; ir,ic,ik = i_p,j_p,k_cur
        else:
            proper=True
    L = np.zeros((n,n),dtype=int)
    for i in range(n):
        for j in range(n):
            for k in range(n):
                if cube[i][j][k]==1: L[i,j]=k+1; break
    return L

def all_ls_4():
    results = []
    L = [[0]*4 for _ in range(4)]
    def ok(r,c,v):
        for j2 in range(c):
            if L[r][j2]==v: return False
        for i2 in range(r):
            if L[i2][c]==v: return False
        return True
    def solve(pos):
        if pos==16:
            results.append([row[:] for row in L]); return
        r,c = pos//4, pos%4
        for v in range(1,5):
            if ok(r,c,v): L[r][c]=v; solve(pos+1); L[r][c]=0
    solve(0)
    return results

# ── Proper SNF via iterated GCD (Kannan-Bachem style) ────────────

def snf(A_np):
    """Smith Normal Form for integer matrix. Returns invariant factors."""
    M = [[int(A_np[i,j]) for j in range(A_np.shape[1])] for i in range(A_np.shape[0])]
    rows, cols = len(M), len(M[0])
    m = min(rows, cols)
    
    for k in range(m):
        MAX_ITER = 1000
        itr = 0
        while itr < MAX_ITER:
            itr += 1
            # Find smallest nonzero in submatrix M[k:, k:]
            best = None; best_val = None
            for i in range(k, rows):
                for j in range(k, cols):
                    if M[i][j] != 0:
                        if best_val is None or abs(M[i][j]) < best_val:
                            best = (i, j); best_val = abs(M[i][j])
            if best is None:
                break
            
            # Move best to (k, k)
            bi, bj = best
            if bi != k:
                M[k], M[bi] = M[bi], M[k]
            if bj != k:
                for i in range(rows):
                    M[i][k], M[i][bj] = M[i][bj], M[i][k]
            
            if M[k][k] < 0:
                for j in range(cols):
                    M[k][j] = -M[k][j]
            
            # Row eliminate
            changed = False
            for i in range(k+1, rows):
                if M[i][k] != 0:
                    q = M[i][k] // M[k][k]
                    for j in range(cols):
                        M[i][j] -= q * M[k][j]
                    if M[i][k] != 0:
                        changed = True
            
            # Column eliminate
            for j in range(k+1, cols):
                if M[k][j] != 0:
                    q = M[k][j] // M[k][k]
                    for i in range(rows):
                        M[i][j] -= q * M[i][k]
                    if M[k][j] != 0:
                        changed = True
            
            if not changed:
                # Check divisibility
                divides_all = True
                for i in range(k+1, rows):
                    for j in range(k+1, cols):
                        if M[i][j] != 0 and M[i][j] % M[k][k] != 0:
                            divides_all = False
                            # Fix: add row i to row k
                            for jj in range(cols):
                                M[k][jj] += M[i][jj]
                            changed = True
                            break
                    if not divides_all:
                        break
                if divides_all:
                    break  # Done with position k
    
    return [abs(M[i][i]) if i < rows and i < cols and abs(M[i][i]) > 0 else 0 for i in range(m)]

# ── Compound matrix rank (for proving bounds) ────────────────────

def compound_rank_F2(A_np, p_power):
    """Compute rank of A mod p_power over Z/p_power*Z for p=2.
    Uses: #{d_i with v2 < t} = rank_2^t(A)."""
    # For Z/2^t: row reduce with modular arithmetic
    mod = p_power
    B = [[int(A_np[i,j]) % mod for j in range(A_np.shape[1])] for i in range(A_np.shape[0])]
    rows, cols = len(B), len(B[0])
    pr = 0
    for col in range(cols):
        # Find pivot with smallest v2
        best_row = None; best_v = 999
        for row in range(pr, rows):
            if B[row][col] % mod != 0:
                vv = v2(B[row][col])
                if vv < best_v:
                    best_v = vv; best_row = row
        if best_row is None:
            continue
        B[pr], B[best_row] = B[best_row], B[pr]
        # Try to make pivot a unit mod 2^t
        pivot = B[pr][col]
        # Eliminate below
        for row in range(rows):
            if row != pr and B[row][col] % mod != 0:
                # Find multipliers to eliminate
                g = gcd(pivot, B[row][col])
                g = gcd(g, mod)
                a = pivot // g
                b = B[row][col] // g
                for j in range(cols):
                    B[row][j] = (a * B[row][j] - b * B[pr][j]) % mod
        pr += 1
    return pr

# ══════════════════════════════════════════════════════════════════════
random.seed(271828)
print("=" * 70)
print("  SNF 2-ADIC PROFILE — Invariant Factors Analysis")
print("=" * 70)

# ── Part 1: n=4 exhaustive ───────────────────────────────────────────
print("\n>>> n=4 EXHAUSTIVE — SNF + rank mod 2^t")
all4 = all_ls_4()
by_kdim_4 = defaultdict(list)

for ls in all4:
    L = np.array(ls, dtype=int); A = get_A(L)
    det_val = bareiss(A)
    if det_val == 0: continue
    kdim = kernel_dim_F2(A)
    factors = snf(A)
    v2_det = v2(det_val)
    v2_profile = tuple(v2(d) for d in factors)
    by_kdim_4[kdim].append((factors, v2_det, v2_profile))

print(f"  Total non-singular: {sum(len(v) for v in by_kdim_4.values())}")
for kd in sorted(by_kdim_4):
    entries = by_kdim_4[kd]
    profiles = Counter(e[2] for e in entries)
    print(f"\n  kdim={kd}: {len(entries)} LS")
    print(f"    v2 profiles: {dict(sorted(profiles.items()))}")
    # Count number of even factors
    n_even = Counter(sum(1 for d in e[0] if d % 2 == 0) for e in entries)
    print(f"    #even_factors: {dict(sorted(n_even.items()))}")
    # min v2(det)
    min_v2 = min(e[1] for e in entries)
    print(f"    min v2(det) = {min_v2}")
    # v2(n/2) = 1 for n=4
    bound = 1 + max(0, kd - 1)
    print(f"    Conj bound = v2(n/2)+max(0,k-1) = 1+{max(0,kd-1)} = {bound}")
    print(f"    {'TIGHT' if min_v2 == bound else f'gap={min_v2-bound}'}")

print("\n  KEY: For n=4, kdim = #even invariant factors (always)?")
for kd in sorted(by_kdim_4):
    entries = by_kdim_4[kd]
    all_match = all(sum(1 for d in e[0] if d % 2 == 0) == kd for e in entries)
    print(f"    kdim={kd}: #even_factors==kdim? {all_match}")

sys.stdout.flush()

# ── Part 2: n=6 JM (50 samples) with SNF ────────────────────────────
print("\n>>> n=6 JM — SNF (50 samples, timeout guard)")
t0 = time.time()
by_kdim_6 = defaultdict(list)
for trial in range(50):
    L = jm_sample(6, 2000); A = get_A(L)
    det_val = bareiss(A)
    if det_val == 0: continue
    kdim = kernel_dim_F2(A)
    t1 = time.time()
    factors = snf(A)
    dt = time.time() - t1
    v2_det = v2(det_val)
    v2_profile = tuple(v2(d) for d in factors)
    by_kdim_6[kdim].append((factors, v2_det, v2_profile))
    if dt > 5:
        print(f"  WARNING: SNF took {dt:.1f}s at trial {trial}")
    if trial == 0:
        print(f"  First sample: kdim={kdim}, factors={factors}, time={dt:.2f}s")

for kd in sorted(by_kdim_6):
    entries = by_kdim_6[kd]
    profiles = Counter(e[2] for e in entries)
    print(f"\n  kdim={kd}: {len(entries)} samples")
    print(f"    v2 profiles: {dict(sorted(profiles.items()))}")
    n_even = Counter(sum(1 for d in e[0] if d % 2 == 0) for e in entries)
    print(f"    #even_factors: {dict(sorted(n_even.items()))}")
    min_v2 = min(e[1] for e in entries)
    bound = 0 + max(0, kd - 1)  # v2(6/2)=v2(3)=0
    print(f"    min v2(det) = {min_v2}, bound = {bound}")

print(f"\n  KEY: For n=6, kdim = #even invariant factors?")
for kd in sorted(by_kdim_6):
    entries = by_kdim_6[kd]
    all_match = all(sum(1 for d in e[0] if d % 2 == 0) == kd for e in entries)
    print(f"    kdim={kd}: #even_factors==kdim? {all_match}")

print(f"  Time: {time.time()-t0:.1f}s")
sys.stdout.flush()

# ── Part 3: n=8 JM (30 samples) with SNF ────────────────────────────
print("\n>>> n=8 JM — SNF (30 samples)")
t0 = time.time()
by_kdim_8 = defaultdict(list)
for trial in range(30):
    L = jm_sample(8, 5000); A = get_A(L)
    det_val = bareiss(A)
    if det_val == 0: continue
    kdim = kernel_dim_F2(A)
    t1 = time.time()
    factors = snf(A)
    dt = time.time() - t1
    v2_det = v2(det_val)
    v2_profile = tuple(v2(d) for d in factors)
    by_kdim_8[kdim].append((factors, v2_det, v2_profile))
    if (trial+1) % 10 == 0:
        print(f"  {trial+1}/30 done ({time.time()-t0:.1f}s)")

for kd in sorted(by_kdim_8):
    entries = by_kdim_8[kd]
    profiles = Counter(e[2] for e in entries)
    print(f"\n  kdim={kd}: {len(entries)} samples")
    for prof, cnt in sorted(profiles.items()):
        print(f"    v2_profile={prof}: {cnt}")
    n_even = Counter(sum(1 for d in e[0] if d % 2 == 0) for e in entries)
    print(f"    #even_factors: {dict(sorted(n_even.items()))}")
    min_v2 = min(e[1] for e in entries)
    bound = 2 + max(0, kd - 1)  # v2(8/2)=v2(4)=2
    print(f"    min v2(det) = {min_v2}, bound = {bound}")

print(f"\n  KEY: For n=8, kdim = #even invariant factors?")
for kd in sorted(by_kdim_8):
    entries = by_kdim_8[kd]
    all_match = all(sum(1 for d in e[0] if d % 2 == 0) == kd for e in entries)
    print(f"    kdim={kd}: #even_factors==kdim? {all_match}")

print(f"  Time: {time.time()-t0:.1f}s")
sys.stdout.flush()

# ── Part 4: THE KEY HYPOTHESIS ───────────────────────────────────────
print("\n" + "=" * 70)
print("  HYPOTHESIS TEST: kdim = #{even elementary divisors}")
print("=" * 70)
print("""
If kdim = dim ker_F2(B) = #{i : 2 | d_i} always, then:
  - Each of the kdim even factors has v2 >= 1
  - det(A) = prod(d_i), so v2(det) = sum of v2(d_i) for the kdim even ones
  - The MINIMUM v2(det) occurs when each even factor has v2 = 1
  - BUT we also know v2(det) >= v2(h) from the universal theorem
  - For n=4: h=2, v2(h)=1. At kdim=1: min_v2=3 (not 1!). So the last
    factor absorbs ALL 2-adic content.
  
The conjecture says: v2(det) >= v2(h) + max(0, kdim-1).
If the even factors have v2 >= v2(h), v2(h), ..., v2(h), 1, 1, ... 
then the bound would be v2(h) + (kdim-1)*1 = v2(h) + kdim - 1.
But the DIVISIBILITY CHAIN d_1|d_2|...|d_{n-1} forces the LAST factor
to have v2 >= all others. So the even factors are the LAST kdim ones,
and d_{n-kdim} | d_{n-kdim+1} | ... | d_{n-1}, all even.

The MINIMUM total v2 subject to:
  - d_{n-kdim} | d_{n-kdim+1} | ... | d_{n-1}
  - each has v2 >= 1
  - the product ... hmm no, the divisibility chain doesn't directly
    minimize the sum.
    
Actually, for the divisibility chain d_{n-kdim}|...|d_{n-1} with each 2|d_i:
  v2(d_{n-kdim}) <= v2(d_{n-kdim+1}) <= ... <= v2(d_{n-1})
  sum >= kdim * v2(d_{n-kdim}) >= kdim * 1 = kdim
  
  But we also need h | det(A), so sum >= v2(h).
  Combined: sum >= max(v2(h), kdim) ... that's weaker than the conjecture for large v2(h).
  
  Actually the conjecture says sum >= v2(h) + max(0, kdim-1) = v2(h) + kdim - 1 for kdim >= 1.
  This can't come just from counting; it requires that ONE of the even factors has v2 >= v2(h).
  
  CLAIM: d_{n-1} (the largest invariant factor) always has v2 >= v2(h).
  If so, then the kdim even factors contribute:
    v2(d_{n-1}) >= v2(h), plus (kdim-1) factors each with v2 >= 1
    sum >= v2(h) + (kdim - 1)
  THAT'S THE CONJECTURE!
""")

# Verify: does the largest invariant factor always have v2 >= v2(h)?
print("  Checking: v2(d_{n-1}) >= v2(h) for all data...")
all_data = [("n=4", by_kdim_4, 4), ("n=6", by_kdim_6, 6), ("n=8", by_kdim_8, 8)]
for label, data, nn in all_data:
    vh = v2(nn // 2)
    fails = 0; total = 0
    for kd, entries in data.items():
        for factors, vd, vp in entries:
            total += 1
            max_factor = factors[-1]
            if v2(max_factor) < vh:
                fails += 1
                print(f"  FAIL: {label} kdim={kd} factors={factors} v2(d_last)={v2(max_factor)} < {vh}")
    print(f"  {label}: {fails}/{total} violations")

print("\n  DONE")
