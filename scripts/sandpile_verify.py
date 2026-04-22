"""
sandpile_verify.py — Independent verification of ALL claims in "analisi agente.md"

Claims to verify:
  C1: |K(A)| = |det(A)| (via SNF)
  C2: n/gcd(n,2) | |K(A)| (Theorem 3.2 / thm:unified)
  C3: v2(|K(A)|) >= v2(n/2) + max(0, k-1) for even n (Theorem 5.7 / thm:sharp-v2)
  C4: n | d_{n-1} for odd prime p|n (Theorem 5.5 / thm:odd-prime-kernel)
  C5: K(A)_cyc ~ (Z/n)^(n-2) for ALL n >= 2 (Corollary 3.7 candidate)
  C6: Cyclic LS does NOT maximize |K(A)| (Observation)
  C7: |K(A)| divides |K(LSG)| (to test — claimed to fail at n=5)
  C8: p-Sylow generalization for odd primes (Conjectured)
  C9: Parity pattern does NOT determine K(A) for n=6 (Observation)
"""
import sys as _sys
try:
    _sys.stdout.reconfigure(encoding="utf-8")
    _sys.stderr.reconfigure(encoding="utf-8")
except Exception:
    pass
import numpy as np
import random
from fractions import Fraction
from math import gcd
from functools import reduce
from collections import Counter

random.seed(42)

# ═══════════════════════════════════════════════════════════════════════
# CORE UTILITIES — all exact arithmetic
# ═══════════════════════════════════════════════════════════════════════

def bareiss_det(M_np):
    """Exact determinant via Bareiss algorithm."""
    n = M_np.shape[0]
    A = [[int(M_np[i][j]) for j in range(n)] for i in range(n)]
    sign = 1
    prev = 1
    for i in range(n):
        if A[i][i] == 0:
            found = False
            for k in range(i+1, n):
                if A[k][i] != 0:
                    A[i], A[k] = A[k], A[i]
                    sign *= -1
                    found = True
                    break
            if not found:
                return 0
        for j in range(i+1, n):
            for k in range(i+1, n):
                A[j][k] = (A[j][k] * A[i][i] - A[j][i] * A[i][k]) // prev
            A[j][i] = 0
        prev = A[i][i]
    return sign * A[n-1][n-1]


def snf_factors(A_np):
    """Smith Normal Form — returns list of diagonal entries (all >= 0)."""
    M = [[int(A_np[i][j]) for j in range(A_np.shape[1])]
         for i in range(A_np.shape[0])]
    rows, cols = len(M), len(M[0])
    m = min(rows, cols)
    
    for k in range(m):
        for iteration in range(3000):
            # Find smallest nonzero below-right of (k,k)
            best = None
            best_val = None
            for i in range(k, rows):
                for j in range(k, cols):
                    if M[i][j] != 0:
                        if best_val is None or abs(M[i][j]) < best_val:
                            best = (i, j)
                            best_val = abs(M[i][j])
            if best is None:
                break
            bi, bj = best
            if bi != k:
                M[k], M[bi] = M[bi], M[k]
            if bj != k:
                for i in range(rows):
                    M[i][k], M[i][bj] = M[i][bj], M[i][k]
            if M[k][k] < 0:
                for j in range(cols):
                    M[k][j] = -M[k][j]
            
            changed = False
            for i in range(k+1, rows):
                if M[i][k] != 0:
                    q = M[i][k] // M[k][k]
                    for j in range(cols):
                        M[i][j] -= q * M[k][j]
                    if M[i][k] != 0:
                        changed = True
            for j in range(k+1, cols):
                if M[k][j] != 0:
                    q = M[k][j] // M[k][k]
                    for i in range(rows):
                        M[i][j] -= q * M[i][k]
                    if M[k][j] != 0:
                        changed = True
            if not changed:
                ok = True
                for i in range(k+1, rows):
                    for j in range(k+1, cols):
                        if M[i][j] != 0 and M[i][j] % M[k][k] != 0:
                            for jj in range(cols):
                                M[k][jj] += M[i][jj]
                            changed = True
                            ok = False
                            break
                    if not ok:
                        break
                if ok:
                    break
    
    return [abs(M[i][i]) if i < min(rows, cols) else 0 for i in range(m)]


def snf_nontrivial(A_np):
    """SNF factors > 1."""
    return [d for d in snf_factors(A_np) if d > 1]


def vp(x, p):
    """p-adic valuation of x."""
    if x == 0:
        return float('inf')
    x = abs(int(x))
    v = 0
    while x % p == 0:
        x //= p
        v += 1
    return v


def kernel_dim_Fp(A_np, p):
    """Dimension of ker(A mod p) over F_p."""
    n_rows, n_cols = A_np.shape
    B = (A_np % p).astype(int)
    # Make positive
    B = ((B % p) + p) % p
    pivot_row = 0
    for col in range(n_cols):
        found = False
        for row in range(pivot_row, n_rows):
            if B[row, col] % p != 0:
                B[[pivot_row, row]] = B[[row, pivot_row]]
                found = True
                break
        if not found:
            continue
        inv_piv = pow(int(B[pivot_row, col]), p - 2, p)  # Fermat's little theorem
        for row in range(n_rows):
            if row != pivot_row and B[row, col] % p != 0:
                factor = (B[row, col] * inv_piv) % p
                B[row] = (B[row] - factor * B[pivot_row]) % p
        pivot_row += 1
    return n_cols - pivot_row


def get_A(L):
    """Standard representation matrix: A_ij = L_ij - L_{i,n-1}."""
    n = L.shape[0]
    return np.array([[int(L[i, j]) - int(L[i, n-1]) for j in range(n-1)]
                     for i in range(n-1)], dtype=int)


def cyclic_ls(n):
    """Cyclic Latin square: L_ij = (i+j) mod n + 1."""
    return np.array([[(i + j) % n + 1 for j in range(n)] for i in range(n)], dtype=int)


def jm_sample(n, steps):
    """Jacobson-Matthews random Latin square sampler."""
    cube = [[[0]*n for _ in range(n)] for _ in range(n)]
    for i in range(n):
        for j in range(n):
            cube[i][j][(i+j) % n] = 1
    proper = True
    ir = ic = ik = 0
    count = 0
    while count < steps or not proper:
        if proper:
            i = random.randint(0, n-1)
            j = random.randint(0, n-1)
            row_ij = cube[i][j]
            k_cur = 0
            while row_ij[k_cur] != 1:
                k_cur += 1
            k = random.randint(0, n-2)
            if k >= k_cur:
                k += 1
            count += 1
        else:
            i, j, k = ir, ic, ik
            cands = [x for x in range(n) if cube[i][j][x] == 1]
            k_cur = cands[random.randint(0, len(cands)-1)]
        
        jc = [jj for jj in range(n) if jj != j and cube[i][jj][k] == 1]
        j_p = jc[random.randint(0, len(jc)-1)]
        ic_ = [ii for ii in range(n) if ii != i and cube[ii][j][k] == 1]
        i_p = ic_[random.randint(0, len(ic_)-1)]
        
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
    
    L = np.zeros((n, n), dtype=int)
    for i in range(n):
        for j in range(n):
            for kk in range(n):
                if cube[i][j][kk] == 1:
                    L[i, j] = kk + 1
                    break
    return L


def lsg_laplacian(L):
    """Laplacian of the Latin square graph LSG(L).
    Vertices: n² cells. Edge if same row, same col, or same symbol."""
    n = L.shape[0]
    N = n * n
    adj = np.zeros((N, N), dtype=int)
    for i1 in range(n):
        for j1 in range(n):
            v1 = i1 * n + j1
            for i2 in range(n):
                for j2 in range(n):
                    if (i1, j1) == (i2, j2):
                        continue
                    v2 = i2 * n + j2
                    if i1 == i2 or j1 == j2 or L[i1, j1] == L[i2, j2]:
                        adj[v1, v2] = 1
    return np.diag(adj.sum(axis=1)) - adj


# ═══════════════════════════════════════════════════════════════════════
print("=" * 70)
print("  INDEPENDENT VERIFICATION — analisi agente.md claims")
print("=" * 70)

# ═══════════════════════════════════════════════════════════════════════
# C1 + C2 + C3 + C4: Core theorems
# ═══════════════════════════════════════════════════════════════════════
print("\n>>> C1-C4: Core theorems verification")
print("    C1: |K(A)| = |det(A)|")
print("    C2: n/gcd(n,2) | |K(A)|")
print("    C3: v2(|K(A)|) >= v2(n/2)+max(0,k-1) [even n]")
print("    C4: For odd prime p|n: p | det(A)")

fails = {"C1": 0, "C2": 0, "C3": 0, "C4": 0}
totals = {"C1": 0, "C2": 0, "C3": 0, "C4": 0}

for n in [3, 4, 5, 6, 7, 8, 9]:
    samples = 30
    for _ in range(samples):
        L = jm_sample(n, 5 * n * n)
        A = get_A(L)
        det_val = bareiss_det(A)
        if det_val == 0:
            continue
        
        factors = snf_factors(A)
        K_order = 1
        for d in factors:
            if d > 0:
                K_order *= d
        
        # C1: |K(A)| = |det(A)|
        totals["C1"] += 1
        if K_order != abs(det_val):
            fails["C1"] += 1
            print(f"  C1 FAIL: n={n}, K_order={K_order}, |det|={abs(det_val)}")
        
        # C2: n/gcd(n,2) | |K(A)|
        totals["C2"] += 1
        h_div = n // gcd(n, 2)
        if K_order % h_div != 0:
            fails["C2"] += 1
            print(f"  C2 FAIL: n={n}, K_order={K_order}, n/gcd(n,2)={h_div}")
        
        # C3: even n only — v2 bound
        if n % 2 == 0:
            totals["C3"] += 1
            k_dim = kernel_dim_Fp(A, 2)
            bound = vp(n // 2, 2) + max(0, k_dim - 1)
            if vp(K_order, 2) < bound:
                fails["C3"] += 1
                print(f"  C3 FAIL: n={n}, v2={vp(K_order,2)}, bound={bound}, kdim={k_dim}")
        
        # C4: odd primes dividing n
        for p in [3, 5, 7]:
            if p > n:
                break
            if n % p == 0:
                totals["C4"] += 1
                if abs(det_val) % p != 0:
                    fails["C4"] += 1
                    print(f"  C4 FAIL: n={n}, p={p}, det={det_val}")

for c in ["C1", "C2", "C3", "C4"]:
    status = "✓ PASS" if fails[c] == 0 else f"✗ {fails[c]} FAILURES"
    print(f"  {c}: {status} ({totals[c]} tests)")


# ═══════════════════════════════════════════════════════════════════════
# C5: K(A)_cyc ~ (Z/n)^(n-2) for ALL n >= 2
# ═══════════════════════════════════════════════════════════════════════
print("\n>>> C5: Corollary 3.7 — K(A)_cyc structure")
print("    Claim: K(A)_cyc ≅ (Z/n)^(n-2) for ALL n >= 2")

c5_pass = True
for n in range(2, 12):
    L = cyclic_ls(n)
    A = get_A(L)
    det_val = bareiss_det(A)
    factors = snf_nontrivial(A)
    
    expected = [n] * (n - 2)
    match = (factors == expected)
    
    tag = "✓" if match else "✗"
    if not match:
        c5_pass = False
    
    struct = " × ".join(f"Z/{d}" for d in factors) if factors else "trivial"
    print(f"  n={n:2d}: K(A) = {struct}, |K|={abs(det_val)}"
          f"  expected (Z/{n})^{n-2}  {tag}")

print(f"  {'✓ CONFIRMED' if c5_pass else '✗ FAILED'} for n=2..11")


# ═══════════════════════════════════════════════════════════════════════
# C6: Cyclic LS does NOT maximize |K(A)|
# ═══════════════════════════════════════════════════════════════════════
print("\n>>> C6: Does cyclic LS maximize |K(A)|?")

for n in [5, 7, 9]:
    L_cyc = cyclic_ls(n)
    det_cyc = abs(bareiss_det(get_A(L_cyc)))
    
    orders = []
    for _ in range(100):
        L = jm_sample(n, 5 * n * n)
        d = bareiss_det(get_A(L))
        if d != 0:
            orders.append(abs(d))
    
    above = sum(1 for x in orders if x > det_cyc)
    max_obs = max(orders) if orders else 0
    print(f"  n={n}: cyclic={det_cyc}=n^{n-2}, max_observed={max_obs}, "
          f"exceed_cyclic={above}/{len(orders)} ({100*above/len(orders):.1f}%)")


# ═══════════════════════════════════════════════════════════════════════
# C7: |K(A)| divides |K(LSG)| 
# ═══════════════════════════════════════════════════════════════════════
print("\n>>> C7: |K(A)| divides |K(LSG)|?")
print("    Claim: divisibility holds for n=3,4 but FAILS at n=5")

for n in [3, 4, 5]:
    n_tests = 8 if n <= 4 else 5
    results = []
    
    for trial in range(n_tests):
        if trial == 0:
            L = cyclic_ls(n)
            label = "cyclic"
        else:
            L = jm_sample(n, 5 * n * n)
            label = f"random_{trial}"
        
        A = get_A(L)
        det_A = bareiss_det(A)
        if det_A == 0:
            continue
        KA_order = abs(det_A)
        
        Lap = lsg_laplacian(L)
        L0 = Lap[:-1, :-1]
        det_L0 = abs(bareiss_det(L0))
        
        divides = (det_L0 % KA_order == 0) if KA_order > 0 else True
        results.append(divides)
        
        if not divides:
            print(f"  n={n} {label}: |K(A)|={KA_order}, |K(LSG)|={det_L0} — DOES NOT DIVIDE")
    
    n_div = sum(results)
    print(f"  n={n}: divides {n_div}/{len(results)}")


# ═══════════════════════════════════════════════════════════════════════
# C8: p-Sylow generalization
# ═══════════════════════════════════════════════════════════════════════
print("\n>>> C8: p-Sylow generalization of Theorem 5.7")
print("    Conjectured: vp(det A) >= vp(n/p) + max(0, kp-1)")
print("    where kp = dim ker(A mod p)")

c8_fails = 0
c8_total = 0

for n, p in [(6, 3), (9, 3), (10, 5), (8, 2), (6, 2)]:
    if n > 10:
        continue
    samples = 40 if n <= 8 else 20
    local_fails = 0
    local_total = 0
    
    for _ in range(samples):
        L = jm_sample(n, 5 * n * n)
        A = get_A(L)
        det_val = bareiss_det(A)
        if det_val == 0:
            continue
        
        kp = kernel_dim_Fp(A, p)
        bound = vp(n // p, p) + max(0, kp - 1)
        actual = vp(abs(det_val), p)
        
        local_total += 1
        c8_total += 1
        if actual < bound:
            local_fails += 1
            c8_fails += 1
            print(f"  FAIL: n={n}, p={p}, vp(det)={actual}, bound={bound}, kp={kp}")
    
    status = "✓" if local_fails == 0 else f"✗ {local_fails}"
    print(f"  n={n}, p={p}: {status} ({local_total} tests)")

print(f"  Overall: {'✓ PASS' if c8_fails == 0 else f'✗ {c8_fails} FAILURES'} ({c8_total} tests)")


# ═══════════════════════════════════════════════════════════════════════
# C9: Parity pattern → K(A) for n=6
# ═══════════════════════════════════════════════════════════════════════
print("\n>>> C9: Does parity pattern determine K(A) for n=6?")

from collections import defaultdict

profile_by_parity = defaultdict(list)
for _ in range(300):
    L = jm_sample(6, 200)
    A = get_A(L)
    d = bareiss_det(A)
    if d == 0:
        continue
    parity_key = tuple((L % 2).flatten())
    profile = tuple(snf_nontrivial(A))
    profile_by_parity[parity_key].append(profile)

conflicts = 0
conflict_examples = []
for key, profiles in profile_by_parity.items():
    if len(profiles) > 1 and len(set(profiles)) > 1:
        conflicts += 1
        if len(conflict_examples) < 3:
            conflict_examples.append((key, set(profiles)))

n_shared = sum(1 for profiles in profile_by_parity.values() if len(profiles) > 1)
print(f"  {len(profile_by_parity)} distinct parity patterns from 300 samples")
print(f"  {n_shared} patterns shared by >=2 LS")
print(f"  Conflicts (same parity, different K(A)): {conflicts}")

if conflicts > 0:
    print("  → Parity pattern does NOT determine K(A) structure")
    for key, profs in conflict_examples[:2]:
        print(f"    Example: profiles = {profs}")
else:
    print("  → Parity pattern DETERMINES K(A) in all tested cases")
    print("    (Note: doc claims conflicts exist — may need more samples)")


# ═══════════════════════════════════════════════════════════════════════
# C10: NEW — 3-Sylow minimum = 9 for n=9 (Theorem 3.2 tightness)
# ═══════════════════════════════════════════════════════════════════════
print("\n>>> C10: 3-Sylow tightness for n=9")

sylow3 = []
for _ in range(100):
    L = jm_sample(9, 10 * 81)
    A = get_A(L)
    d = bareiss_det(A)
    if d == 0:
        continue
    order = abs(d)
    v3 = 0
    while order % 3 == 0:
        order //= 3
        v3 += 1
    sylow3.append(3**v3)

cnt = Counter(sylow3)
print(f"  Distribution: {dict(sorted(cnt.items()))}")
print(f"  Min 3-Sylow = {min(sylow3)}, expected lower bound = 9 = 3^2")
print(f"  Bound tight: {min(sylow3) == 9}")


# ═══════════════════════════════════════════════════════════════════════
print("\n" + "=" * 70)
print("  SUMMARY")
print("=" * 70)
