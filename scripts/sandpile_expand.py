"""
sandpile_expand.py — Expansion & deep verification of Sandpile Group results.

GOALS:
  1. PROVE Cor 3.7: K(A)_cyc ≅ (Z/n)^(n-2) — via explicit SNF of circulant
  2. p-Sylow conjecture: stress test on many (n,p) pairs
  3. C9: exhaustive n=5 reduced LS (56 total) — parity vs K(A)
  4. C7: characterize when |K(A)| ∤ |K(LSG)|
  5. New: SNF profile distribution for n=5 (exhaustive) and n=7
"""
import numpy as np
import random
from math import gcd
from functools import reduce
from fractions import Fraction
from collections import Counter, defaultdict
from itertools import permutations

random.seed(42)

# ═══════════════════════════════════════════════════════════════════════
# UTILITIES (same as sandpile_verify.py)
# ═══════════════════════════════════════════════════════════════════════

def bareiss_det(M_np):
    n = M_np.shape[0]
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

def snf_factors(A_np):
    M = [[int(A_np[i][j]) for j in range(A_np.shape[1])] for i in range(A_np.shape[0])]
    rows, cols = len(M), len(M[0]); m = min(rows, cols)
    for k in range(m):
        for _ in range(3000):
            best = None; best_val = None
            for i in range(k, rows):
                for j in range(k, cols):
                    if M[i][j] != 0:
                        if best_val is None or abs(M[i][j]) < best_val:
                            best = (i, j); best_val = abs(M[i][j])
            if best is None: break
            bi, bj = best
            if bi != k: M[k], M[bi] = M[bi], M[k]
            if bj != k:
                for i in range(rows): M[i][k], M[i][bj] = M[i][bj], M[i][k]
            if M[k][k] < 0:
                for j in range(cols): M[k][j] = -M[k][j]
            changed = False
            for i in range(k+1, rows):
                if M[i][k] != 0:
                    q = M[i][k] // M[k][k]
                    for j in range(cols): M[i][j] -= q * M[k][j]
                    if M[i][k] != 0: changed = True
            for j in range(k+1, cols):
                if M[k][j] != 0:
                    q = M[k][j] // M[k][k]
                    for i in range(rows): M[i][j] -= q * M[i][k]
                    if M[k][j] != 0: changed = True
            if not changed:
                ok = True
                for i in range(k+1, rows):
                    for j in range(k+1, cols):
                        if M[i][j] != 0 and M[i][j] % M[k][k] != 0:
                            for jj in range(cols): M[k][jj] += M[i][jj]
                            ok = False; break
                    if not ok: break
                if ok: break
    return [abs(M[i][i]) for i in range(m)]

def snf_nontrivial(A_np):
    return [d for d in snf_factors(A_np) if d > 1]

def vp(x, p):
    if x == 0: return float('inf')
    x = abs(int(x)); v = 0
    while x % p == 0: x //= p; v += 1
    return v

def kernel_dim_Fp(A_np, p):
    n_rows, n_cols = A_np.shape
    B = ((A_np % p).astype(int) + p) % p
    pivot_row = 0
    for col in range(n_cols):
        found = False
        for row in range(pivot_row, n_rows):
            if B[row, col] % p != 0:
                B[[pivot_row, row]] = B[[row, pivot_row]]; found = True; break
        if not found: continue
        inv_piv = pow(int(B[pivot_row, col]), p-2, p)
        for row in range(n_rows):
            if row != pivot_row and B[row, col] % p != 0:
                factor = (B[row, col] * inv_piv) % p
                B[row] = (B[row] - factor * B[pivot_row]) % p
        pivot_row += 1
    return n_cols - pivot_row

def get_A(L):
    n = L.shape[0]
    return np.array([[int(L[i,j]) - int(L[i,n-1]) for j in range(n-1)] for i in range(n-1)], dtype=int)

def cyclic_ls(n):
    return np.array([[(i+j) % n + 1 for j in range(n)] for i in range(n)], dtype=int)

def jm_sample(n, steps):
    cube = [[[0]*n for _ in range(n)] for _ in range(n)]
    for i in range(n):
        for j in range(n): cube[i][j][(i+j) % n] = 1
    proper = True; ir = ic = ik = 0; count = 0
    while count < steps or not proper:
        if proper:
            i = random.randint(0, n-1); j = random.randint(0, n-1)
            k_cur = 0
            while cube[i][j][k_cur] != 1: k_cur += 1
            k = random.randint(0, n-2)
            if k >= k_cur: k += 1
            count += 1
        else:
            i, j, k = ir, ic, ik
            cands = [x for x in range(n) if cube[i][j][x] == 1]
            k_cur = cands[random.randint(0, len(cands)-1)]
        jc = [jj for jj in range(n) if jj != j and cube[i][jj][k] == 1]
        j_p = jc[random.randint(0, len(jc)-1)]
        ic_ = [ii for ii in range(n) if ii != i and cube[ii][j][k] == 1]
        i_p = ic_[random.randint(0, len(ic_)-1)]
        cube[i][j][k] += 1; cube[i][j][k_cur] -= 1
        cube[i][j_p][k] -= 1; cube[i][j_p][k_cur] += 1
        cube[i_p][j][k] -= 1; cube[i_p][j][k_cur] += 1
        cube[i_p][j_p][k] += 1; cube[i_p][j_p][k_cur] -= 1
        if cube[i_p][j_p][k_cur] == -1:
            proper = False; ir, ic, ik = i_p, j_p, k_cur
        else: proper = True
    L = np.zeros((n,n), dtype=int)
    for i in range(n):
        for j in range(n):
            for kk in range(n):
                if cube[i][j][kk] == 1: L[i,j] = kk+1; break
    return L

def lsg_laplacian(L):
    n = L.shape[0]; N = n*n
    adj = np.zeros((N, N), dtype=int)
    for i1 in range(n):
        for j1 in range(n):
            v1 = i1*n + j1
            for i2 in range(n):
                for j2 in range(n):
                    if (i1, j1) == (i2, j2): continue
                    v2 = i2*n + j2
                    if i1 == i2 or j1 == j2 or L[i1, j1] == L[i2, j2]:
                        adj[v1, v2] = 1
    return np.diag(adj.sum(axis=1)) - adj


# ═══════════════════════════════════════════════════════════════════════
# Enumerate ALL reduced Latin squares of order n (small n only)
# ═══════════════════════════════════════════════════════════════════════

def enumerate_reduced_ls(n):
    """Enumerate all reduced LS of order n (first row = 1..n, first col = 1..n).
    Only feasible for n <= 6."""
    results = []
    L = np.zeros((n, n), dtype=int)
    # Fix first row
    for j in range(n): L[0, j] = j + 1
    # Fix first column
    for i in range(n): L[i, 0] = i + 1

    def valid(r, c, v):
        for j2 in range(c):
            if L[r][j2] == v: return False
        for i2 in range(r):
            if L[i2][c] == v: return False
        return True

    cells = [(i, j) for i in range(1, n) for j in range(1, n)]

    def fill(idx):
        if idx == len(cells):
            results.append(L.copy())
            return
        r, c = cells[idx]
        for v in range(1, n+1):
            if valid(r, c, v):
                L[r][c] = v
                fill(idx + 1)
                L[r][c] = 0
    fill(0)
    return results


# ═══════════════════════════════════════════════════════════════════════
print("=" * 70)
print("  SANDPILE GROUP — EXPANSION & DEEP VERIFICATION")
print("=" * 70)


# ═══════════════════════════════════════════════════════════════════════
# BLOCK 1: PROVE Corollary 3.7 — circulant SNF structure
# ═══════════════════════════════════════════════════════════════════════
print("\n" + "=" * 70)
print("BLOCK 1: Corollary 3.7 — WHY is K(A)_cyc = (Z/n)^(n-2)?")
print("=" * 70)

print("""
PROOF SKETCH (circulant theory):
  For cyclic LS: L_ij = (i+j) mod n + 1, so A_ij = L_ij - L_{i,n-1}
  = ((i+j) mod n) - ((i+n-1) mod n).

  A is a circulant matrix: A = circ(a_0, a_1, ..., a_{n-2})
  where a_j = (j mod n) - ((n-1) mod n) for j = 0..n-2.

  For a circulant C = circ(c_0,...,c_{m-1}) over Z, the SNF has
  d_k = gcd of all k×k minors. For the SPECIFIC circulant A from
  the cyclic LS, we need to verify that ALL invariant factors = n.

  This means: gcd of all 1×1 minors (= entries) = 1,
              gcd of all 2×2 minors / d_1 = ?,
              ...
              and the product = n^(n-2) = det(A).

  If all d_i = n, then d_1 = n means gcd(entries) = n.
  But entries are a_j = (j mod n) - (n-1 mod n), which for j=0
  gives 0 - (n-1) = -(n-1) ≡ 1 mod n. So gcd(entries) = ?
  Let's check explicitly.
""")

# Verify: what are the entries of A_cyc?
for n in [5, 7, 9]:
    L = cyclic_ls(n)
    A = get_A(L)
    entries = sorted(set(A.flatten()))
    g = reduce(gcd, [abs(x) for x in A.flatten() if x != 0])
    print(f"  n={n}: A entries = {entries}, gcd(entries) = {g}")
    # Show first row
    print(f"         first row = {list(A[0])}")

print("\n  → gcd(entries) = 1 for all n, so d_1 = 1??")
print("  → But SNF shows d_i = n. Resolution: SNF factors ≠ entry gcd.")
print("  → d_1 | d_2 | ... | d_{n-2}, product = n^{n-2}.")
print("  → If d_1 = n, then ALL minors of size 1 divisible by n. FALSE above.")
print("  → So d_1 < n. Let me re-check...\n")

# Actually compute d_1 = gcd of all entries
# d_2 = gcd of all 2x2 minors / d_1
# etc.
for n in [3, 5, 7]:
    L = cyclic_ls(n)
    A = get_A(L)
    factors = snf_factors(A)
    nontrivial = [d for d in factors if d > 1]
    print(f"  n={n}: full SNF = {factors}")
    print(f"         nontrivial = {nontrivial}")
    # d_1 = first SNF factor = gcd of all entries
    g = reduce(gcd, [abs(int(x)) for x in A.flatten() if x != 0])
    print(f"         gcd(all entries) = {g}, SNF d_1 = {factors[0]}")

print("""
  KEY INSIGHT: For n=5, SNF = [1, 1, 5, 5, ... wait no:
  nontrivial = [5, 5, 5] means SNF = [1, 5, 5, 5] (4 factors for 4x4 A).
  So d_1=1, d_2=5, d_3=5, d_4=... wait.
""")

# Let me be precise
for n in [3, 4, 5, 6, 7]:
    L = cyclic_ls(n)
    A = get_A(L)
    factors = snf_factors(A)
    print(f"  n={n}: A is {n-1}x{n-1}, SNF = {factors}")
    print(f"         → K(A) = {' × '.join(f'Z/{d}' for d in factors if d > 1)}")
    # Check: how many 1s at the start?
    n_ones = sum(1 for d in factors if d == 1)
    n_n = sum(1 for d in factors if d == n)
    print(f"         → {n_ones} ones, {n_n} copies of n={n}")
    if n_ones + n_n == n - 1:
        print(f"         → Pattern: [{n_ones} × 1, {n_n} × {n}] ✓")
    print()


# ═══════════════════════════════════════════════════════════════════════
# BLOCK 1B: ALGEBRAIC PROOF of Cor 3.7
# ═══════════════════════════════════════════════════════════════════════
print("─" * 70)
print("BLOCK 1B: Algebraic proof verification")
print("─" * 70)

print("""
For cyclic LS of order n, A is (n-1)×(n-1) with:
  A_ij = ((i+j) mod n) - ((i+n-1) mod n)    (0-indexed, i,j = 0..n-2)

Eigenvalues of A (on C^{n-1}): the DFT diagonalizes any circulant.
A has eigenvalue λ_k = Σ_{j=0}^{n-2} a_j · ω^{jk}  for k=1..n-1
where ω = e^{2πi/n}.

For the cyclic case, det(A) = ±n^{n-2} (Proposition 3.6 of paper).

KEY THEOREM (circulant SNF):
If C is an integer circulant of size m, and ALL eigenvalues λ_k
satisfy |λ_k| = same value q, then SNF has all nontrivial factors = q
(when q is a prime power times units).

For our A: eigenvalues are λ_k = n/(ω^k - 1) scaled.
|λ_k| varies with k in general! So the "equal modulus" argument
doesn't directly apply.

Alternative: use the EXPLICIT eigenvalue computation.
λ_k = n/(ω^k - 1) for the E matrix eigenvalues on V_std.
Since E = G·A where G has eigenvalue n on V_std,
A has eigenvalue λ_k(A) = λ_k(E)/n = 1/(ω^k - 1) × ... hmm.

Let me just verify the SNF pattern numerically for larger n.
""")

# Extend verification to n=12..15
for n in range(12, 16):
    L = cyclic_ls(n)
    A = get_A(L)
    factors = snf_factors(A)
    nontrivial = [d for d in factors if d > 1]
    all_n = all(d == n for d in nontrivial)
    n_ones = sum(1 for d in factors if d == 1)
    tag = "✓" if all_n and len(nontrivial) == n - 2 else "✗"
    print(f"  n={n:2d}: SNF nontrivial = {nontrivial[:5]}{'...' if len(nontrivial)>5 else ''}"
          f"  → (Z/{n})^{len(nontrivial)}  {tag}")


# ═══════════════════════════════════════════════════════════════════════
# BLOCK 2: p-Sylow STRESS TEST
# ═══════════════════════════════════════════════════════════════════════
print("\n" + "=" * 70)
print("BLOCK 2: p-Sylow generalization — stress test")
print("    Conjecture: vp(det A) >= vp(n/p) + max(0, kp-1)")
print("=" * 70)

sylow_fails = 0
sylow_total = 0

test_cases = [
    # (n, p, samples)
    (6, 2, 80), (6, 3, 80),
    (8, 2, 80),
    (9, 3, 80),
    (10, 2, 40), (10, 5, 40),
    (12, 2, 30), (12, 3, 30),
    (14, 2, 20), (14, 7, 20),
    (15, 3, 20), (15, 5, 20),
]

for n, p, samples in test_cases:
    local_fails = 0
    local_total = 0
    kp_dist = Counter()
    
    for _ in range(samples):
        L = jm_sample(n, 5 * n * n)
        A = get_A(L)
        det_val = bareiss_det(A)
        if det_val == 0: continue
        
        kp = kernel_dim_Fp(A, p)
        bound = vp(n // p, p) + max(0, kp - 1)
        actual = vp(abs(det_val), p)
        
        kp_dist[kp] += 1
        local_total += 1
        sylow_total += 1
        if actual < bound:
            local_fails += 1
            sylow_fails += 1
            print(f"  ✗ FAIL: n={n}, p={p}, vp(det)={actual}, bound={bound}, kp={kp}")
    
    status = "✓" if local_fails == 0 else f"✗ {local_fails}"
    print(f"  n={n:2d}, p={p}: {status} ({local_total} tests)  kp_dist={dict(sorted(kp_dist.items()))}")

print(f"\n  OVERALL: {'✓ ALL PASS' if sylow_fails == 0 else f'✗ {sylow_fails} FAILURES'}"
      f" ({sylow_total} tests)")


# ═══════════════════════════════════════════════════════════════════════
# BLOCK 3: EXHAUSTIVE n=5 — all 56 reduced LS
# ═══════════════════════════════════════════════════════════════════════
print("\n" + "=" * 70)
print("BLOCK 3: Exhaustive n=5 — all 56 reduced LS")
print("=" * 70)

ls5 = enumerate_reduced_ls(5)
print(f"  Enumerated {len(ls5)} reduced LS of order 5 (expected: 56)")

# Full SNF profile for each
profile_dist_5 = Counter()
order_dist_5 = Counter()
parity_to_profiles_5 = defaultdict(set)

for L in ls5:
    A = get_A(L)
    det_val = bareiss_det(A)
    if det_val == 0:
        profile_dist_5[("zero",)] += 1
        continue
    profile = tuple(snf_nontrivial(A))
    profile_dist_5[profile] += 1
    order_dist_5[abs(det_val)] += 1
    
    parity_key = tuple((L % 2).flatten())
    parity_to_profiles_5[parity_key].add(profile)

print(f"\n  SNF profile distribution:")
for profile, cnt in sorted(profile_dist_5.items(), key=lambda x: -x[1]):
    struct = ' × '.join(f'Z/{d}' for d in profile) if profile[0] != "zero" else "ZERO"
    order = 1
    for d in profile:
        if d != "zero": order *= d
    print(f"    [{cnt:3d}] {struct:30s}  |K|={order}")

print(f"\n  |K(A)| distribution:")
for order, cnt in sorted(order_dist_5.items()):
    print(f"    |K|={order:6d}: {cnt} LS")

# Max |K|
max_order = max(order_dist_5.keys())
max_cyc = 5**(5-2)
print(f"\n  Max |K(A)| = {max_order}, cyclic = {max_cyc} = 5^3")
print(f"  Cyclic is {'maximum' if max_order == max_cyc else 'NOT maximum'}")

# Parity pattern consistency at n=5 (n odd, so parity is less meaningful,
# but test anyway)
conflicts_5 = 0
for key, profiles in parity_to_profiles_5.items():
    if len(profiles) > 1:
        conflicts_5 += 1

print(f"\n  Parity patterns with multiple K(A) profiles: {conflicts_5}")
if conflicts_5 > 0:
    print(f"  → Parity does NOT determine K(A) for n=5")
else:
    print(f"  → Parity determines K(A) for all n=5 reduced LS")


# ═══════════════════════════════════════════════════════════════════════
# BLOCK 4: C7 — divisibility analysis at n=5
# ═══════════════════════════════════════════════════════════════════════
print("\n" + "=" * 70)
print("BLOCK 4: |K(A)| divides |K(LSG)| — detailed n=5 analysis")
print("=" * 70)

# Use a subset of the exhaustive n=5 LS
n5_divides = 0
n5_fails = 0
n5_fail_orders = []

# LSG computation is O(n^4) ≈ 625 per grid — feasible for 56 grids
for idx, L in enumerate(ls5):
    A = get_A(L)
    det_A = bareiss_det(A)
    if det_A == 0: continue
    KA_order = abs(det_A)
    
    Lap = lsg_laplacian(L)
    L0 = Lap[:-1, :-1]
    det_L0 = abs(bareiss_det(L0))
    
    if det_L0 == 0:
        print(f"  Grid {idx}: K(LSG) order = 0 (tree)")
        continue
    
    if KA_order > 0 and det_L0 % KA_order == 0:
        n5_divides += 1
    else:
        n5_fails += 1
        n5_fail_orders.append(KA_order)
        if n5_fails <= 5:
            print(f"  Grid {idx}: |K(A)|={KA_order}, |K(LSG)|={det_L0}, "
                  f"remainder={det_L0 % KA_order if KA_order > 0 else 'N/A'}")

print(f"\n  n=5 exhaustive: divides {n5_divides}/{n5_divides+n5_fails}"
      f" ({100*n5_divides/(n5_divides+n5_fails):.1f}%)")
if n5_fails > 0:
    print(f"  Failed |K(A)| values: {Counter(n5_fail_orders)}")
    print(f"  → Divisibility FAILS for {n5_fails} reduced LS at n=5")
else:
    print(f"  → Divisibility holds for ALL n=5 reduced LS!")


# ═══════════════════════════════════════════════════════════════════════
# BLOCK 5: Exhaustive n=4 for parity/K(A)
# ═══════════════════════════════════════════════════════════════════════
print("\n" + "=" * 70)
print("BLOCK 5: Exhaustive n=4 — all 4 reduced LS")
print("=" * 70)

ls4 = enumerate_reduced_ls(4)
print(f"  Enumerated {len(ls4)} reduced LS of order 4 (expected: 4)")

for idx, L in enumerate(ls4):
    A = get_A(L)
    det_val = bareiss_det(A)
    profile = tuple(snf_nontrivial(A))
    struct = ' × '.join(f'Z/{d}' for d in profile) if profile else "trivial"
    parity = tuple((L % 2).flatten())
    kdim = kernel_dim_Fp(A, 2)
    print(f"  LS {idx}: |K|={abs(det_val)}, K(A)={struct}, "
          f"kdim_F2={kdim}")
    print(f"    L = {L.tolist()}")

    # LSG test
    Lap = lsg_laplacian(L)
    L0 = Lap[:-1, :-1]
    det_L0 = abs(bareiss_det(L0))
    KA = abs(det_val)
    divides = (det_L0 % KA == 0) if KA > 0 else True
    print(f"    |K(LSG)|={det_L0}, divides: {divides}")


# ═══════════════════════════════════════════════════════════════════════
# BLOCK 6: NEW — what determines K(A) profile?
# ═══════════════════════════════════════════════════════════════════════
print("\n" + "=" * 70)
print("BLOCK 6: What invariant determines the K(A) profile?")
print("    Testing: cycle structure of row permutations")
print("=" * 70)

def cycle_type(perm):
    """Cycle type of a permutation (as sorted tuple)."""
    n = len(perm)
    visited = [False] * n
    cycles = []
    for i in range(n):
        if visited[i]: continue
        length = 0; j = i
        while not visited[j]:
            visited[j] = True; j = perm[j]; length += 1
        cycles.append(length)
    return tuple(sorted(cycles, reverse=True))

# For n=5 exhaustive
print("\n  n=5: cycle structure of L (as product of row permutations)")
cycle_to_profiles = defaultdict(set)
cycle_to_orders = defaultdict(list)

for L in ls5:
    A = get_A(L)
    det_val = bareiss_det(A)
    if det_val == 0: continue
    profile = tuple(snf_nontrivial(A))
    
    # Compute cycle types of each row (as permutation of {0,...,n-1})
    row_cycles = []
    for i in range(5):
        perm = [L[i, j] - 1 for j in range(5)]
        row_cycles.append(cycle_type(perm))
    cycle_key = tuple(sorted(row_cycles))
    
    cycle_to_profiles[cycle_key].add(profile)
    cycle_to_orders[cycle_key].append(abs(det_val))

for ckey in sorted(cycle_to_profiles.keys()):
    profiles = cycle_to_profiles[ckey]
    orders = cycle_to_orders[ckey]
    print(f"  Cycles {ckey}: {len(orders)} LS, profiles={profiles}, "
          f"orders={sorted(set(orders))}")


# ═══════════════════════════════════════════════════════════════════════
# BLOCK 7: NEW — is |K(A)| maximized at a specific structure?
# ═══════════════════════════════════════════════════════════════════════
print("\n" + "=" * 70)
print("BLOCK 7: Which LS structure maximizes |K(A)|?")
print("=" * 70)

# n=5: find the LS with maximum |K(A)| among the 56
max_det = 0
max_L = None
for L in ls5:
    A = get_A(L)
    d = abs(bareiss_det(A))
    if d > max_det:
        max_det = d
        max_L = L.copy()

print(f"  n=5: max |K(A)| = {max_det}")
print(f"  Maximizer:\n{max_L}")
A_max = get_A(max_L)
print(f"  SNF = {snf_nontrivial(A_max)}")
print(f"  Cyclic |K| = {5**3} = 125")

# Is the maximizer unique?
max_count = sum(1 for L in ls5 if abs(bareiss_det(get_A(L))) == max_det)
print(f"  Number of reduced LS with max |K|: {max_count}")

# n=7: sample to find larger
print(f"\n  n=7: searching for max |K(A)| (sampling)...")
best7 = 0
best7_profile = None
for _ in range(500):
    L = jm_sample(7, 5*49)
    d = abs(bareiss_det(get_A(L)))
    if d > best7:
        best7 = d
        best7_profile = tuple(snf_nontrivial(get_A(L)))

print(f"  n=7: max observed |K(A)| = {best7}")
print(f"  SNF profile = {best7_profile}")
print(f"  Cyclic |K| = {7**5} = {7**5}")
print(f"  Ratio max/cyclic = {best7/7**5:.4f}")


# ═══════════════════════════════════════════════════════════════════════
print("\n" + "=" * 70)
print("  FINAL SUMMARY")
print("=" * 70)
