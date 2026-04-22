"""
Explore the REOPENED question: for n ≡ 0 (mod 4), what characterizes n² | det(E_std)?

Conjectures tested:
  C1: dim ker_F2(A) = 1 (only **1**) ⟺ v₂(det E) = v₂(n²) - 1  (i.e. n² ∤ det E)
  C2: dim ker_F2(A) ≥ 2  ⟺ n² | det(E_std)
  C3: ker(A mod 4) non-trivial ⟺ n² | det(E_std)
  C4: Smith normal form analysis over Z — look at elementary divisors mod 2^a
"""
import numpy as np
from collections import Counter, defaultdict

# ── Core helpers (from run_T3.py) ──────────────────────────────────────
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

# ── F₂ linear algebra ─────────────────────────────────────────────────
def kernel_dim_F2(M):
    """Compute dim ker(M mod 2) over F₂ using Gaussian elimination."""
    n_rows, n_cols = M.shape
    B = M.copy() % 2
    pivot_cols = []
    row = 0
    for col in range(n_cols):
        found = None
        for r in range(row, n_rows):
            if B[r, col] % 2 == 1:
                found = r; break
        if found is None:
            continue
        B[[row, found]] = B[[found, row]]
        pivot_cols.append(col)
        for r in range(n_rows):
            if r != row and B[r, col] % 2 == 1:
                B[r] = (B[r] + B[row]) % 2
        row += 1
    rank = len(pivot_cols)
    return n_cols - rank  # nullity = n_cols - rank

def has_one_in_kernel_F2(M):
    """Check if 1-vector is in ker(M mod 2)."""
    n_rows, n_cols = M.shape
    v = M.sum(axis=1) % 2
    return all(x == 0 for x in v)

# ── mod 4 kernel analysis ─────────────────────────────────────────────
def kernel_dim_mod4(M):
    """Check if M*x ≡ 0 (mod 4) has solutions beyond trivial multiples of mod-2 kernel.
    Returns number of independent vectors x with M*x ≡ 0 mod 4."""
    n_rows, n_cols = M.shape
    # Enumerate: for each candidate vector x in {0,1,2,3}^n_cols, check Mx ≡ 0 mod 4
    # Too expensive for n=7 (4^7 = 16384) — use Gaussian elimination mod 4 instead
    # Actually: we just check specific vectors
    
    # Check 1-vector
    v1 = M.sum(axis=1) % 4
    one_in_ker4 = all(x == 0 for x in v1)
    
    # Check 2*1-vector
    v2 = (2 * M.sum(axis=1)) % 4
    two_one_in_ker4 = all(x == 0 for x in v2)
    
    return one_in_ker4, two_one_in_ker4

# ── Smith Normal Form (elementary divisors) ────────────────────────────
def smith_nf_diag(M):
    """Compute Smith Normal Form diagonal of integer matrix M.
    Returns sorted list of diagonal entries (elementary divisors)."""
    from math import gcd
    m, n = M.shape
    B = M.copy().tolist()
    B = [[int(x) for x in row] for row in B]
    
    size = min(m, n)
    for s in range(size):
        # Find nonzero pivot
        found = False
        for i in range(s, m):
            for j in range(s, n):
                if B[i][j] != 0:
                    # Swap rows
                    B[s], B[i] = B[i], B[s]
                    # Swap cols
                    for r in range(m):
                        B[r][s], B[r][j] = B[r][j], B[r][s]
                    found = True
                    break
            if found: break
        if not found:
            break
        
        # Reduce until pivot divides all entries in its row and column
        changed = True
        while changed:
            changed = False
            # Column operations
            for i in range(s+1, m):
                if B[i][s] != 0:
                    if B[s][s] != 0 and B[i][s] % B[s][s] == 0:
                        q = B[i][s] // B[s][s]
                        for j in range(n):
                            B[i][j] -= q * B[s][j]
                    else:
                        g = gcd(abs(B[s][s]), abs(B[i][s]))
                        a_, b_ = B[s][s] // g, B[i][s] // g
                        new_s = [a_ * B[i][j] - b_ * B[s][j] for j in range(n)]
                        new_i = [B[s][j] + B[i][j] for j in range(n)]  # not quite right
                        # Simple: just swap if |B[i][s]| < |B[s][s]|
                        if abs(B[i][s]) < abs(B[s][s]) and B[i][s] != 0:
                            B[s], B[i] = B[i], B[s]
                            changed = True
                        elif B[i][s] != 0:
                            q = B[i][s] // B[s][s]
                            for j in range(n):
                                B[i][j] -= q * B[s][j]
                            if B[i][s] != 0:
                                changed = True
            # Row operations
            for j in range(s+1, n):
                if B[s][j] != 0:
                    if B[s][s] != 0 and B[s][j] % B[s][s] == 0:
                        q = B[s][j] // B[s][s]
                        for i in range(m):
                            B[i][j] -= q * B[i][s]
                    else:
                        if abs(B[s][j]) < abs(B[s][s]) and B[s][j] != 0:
                            for i in range(m):
                                B[i][s], B[i][j] = B[i][j], B[i][s]
                            changed = True
                        elif B[s][j] != 0:
                            q = B[s][j] // B[s][s]
                            for i in range(m):
                                B[i][j] -= q * B[i][s]
                            if B[s][j] != 0:
                                changed = True
    
    diag = [abs(B[i][i]) for i in range(size)]
    return diag

def v2_of_smith(M):
    """Return list of v_2 of each elementary divisor."""
    diag = smith_nf_diag(M)
    return [v_p(d, 2) for d in diag]

# ═══════════════════════════════════════════════════════════════════════
# MAIN EXPLORATION
# ═══════════════════════════════════════════════════════════════════════
n = 8
L_cyc = np.array([[(i+j) % n + 1 for j in range(n)] for i in range(n)], dtype=np.int64)
L_nc8 = np.array([
    [1,2,3,4,5,6,7,8],[2,1,4,3,6,5,8,7],[3,4,1,2,7,8,5,6],
    [4,3,2,1,8,7,6,5],[5,6,7,8,1,2,3,4],[6,5,8,7,2,1,4,3],
    [7,8,5,6,3,4,1,2],[8,7,6,5,4,3,2,1]], dtype=np.int64)
bases = [L_cyc, L_nc8]

print("=" * 72)
print("EXPLORATION: n=8, characterizing n² | det(E_std) for n ≡ 0 (mod 4)")
print("=" * 72)
print(f"n={n}, n²={n**2}, v₂(n²)={v_p(n**2,2)}, n²/2={n**2//2}, v₂(n²/2)={v_p(n**2//2,2)}")
print()

# Collect data
results = []
total = 0

for bi, Lb in enumerate(bases):
    num = 1000 if bi == 0 else 500
    for seed in range(num):
        rng = np.random.default_rng(seed + bi * 10000)
        L = swap_chain(Lb, rng, 800)
        
        A = np.zeros((n-1, n-1), dtype=np.int64)
        for i in range(n-1):
            for j in range(n-1):
                A[i,j] = int(L[i,j]) - int(L[i,n-1])
        
        dA = bareiss(A.tolist())
        dE = n * dA
        if dE == 0:
            continue
        
        total += 1
        v2_dE = v_p(dE, 2)
        v2_dA = v_p(dA, 2)
        
        # F₂ analysis
        kdim = kernel_dim_F2(A)
        one_in_ker = has_one_in_kernel_F2(A)
        
        # mod 4 analysis
        one_ker4, two_one_ker4 = kernel_dim_mod4(A)
        
        results.append({
            'v2_dE': v2_dE,
            'v2_dA': v2_dA,
            'kdim_F2': kdim,
            'one_in_ker': one_in_ker,
            'one_ker4': one_ker4,
            'n2_div': (dE % (n**2) == 0),
            'bi': bi, 'seed': seed,
        })

print(f"Total non-degenerate LS sampled: {total}")
print()

# ── ANALYSIS 1: v₂ distribution ───────────────────────────────────────
print("─" * 50)
print("ANALYSIS 1: v₂(det E_std) distribution")
print("─" * 50)
v2_dist = Counter(r['v2_dE'] for r in results)
for k in sorted(v2_dist):
    pct = 100 * v2_dist[k] / total
    marker = " *** < v₂(n²)=6" if k < 6 else ""
    print(f"  v₂ = {k:3d}: {v2_dist[k]:5d}  ({pct:5.1f}%){marker}")

n2_pass = sum(1 for r in results if r['n2_div'])
n2_fail = total - n2_pass
print(f"\nn² | det E: {n2_pass}/{total} ({100*n2_pass/total:.1f}%) PASS")
print(f"n² ∤ det E: {n2_fail}/{total} ({100*n2_fail/total:.1f}%) FAIL")

# ── ANALYSIS 2: F₂ kernel dimension ──────────────────────────────────
print()
print("─" * 50)
print("ANALYSIS 2: dim ker_F₂(A mod 2)")
print("─" * 50)
kdim_dist = Counter(r['kdim_F2'] for r in results)
for k in sorted(kdim_dist):
    pct = 100 * kdim_dist[k] / total
    print(f"  dim ker = {k}: {kdim_dist[k]:5d}  ({pct:5.1f}%)")

print(f"\n1 ∈ ker(A mod 2): {sum(1 for r in results if r['one_in_ker'])}/{total}")

# ── ANALYSIS 3: CROSS-TABULATION — the key test ──────────────────────
print()
print("─" * 50)
print("ANALYSIS 3: CROSS-TABULATION — ker dim vs n²-divisibility")
print("─" * 50)
print("  Conjecture C1: dim ker_F₂ = 1 ⟺ n² ∤ det E")
print("  Conjecture C2: dim ker_F₂ ≥ 2 ⟺ n² | det E")
print()

# Cross-tab
tab = defaultdict(lambda: {'pass': 0, 'fail': 0})
for r in results:
    key = r['kdim_F2']
    if r['n2_div']:
        tab[key]['pass'] += 1
    else:
        tab[key]['fail'] += 1

print(f"  {'dim ker':>8s}  {'n²|det (PASS)':>14s}  {'n²∤det (FAIL)':>14s}  {'Total':>8s}")
print(f"  {'─'*8}  {'─'*14}  {'─'*14}  {'─'*8}")
for k in sorted(tab):
    p, f = tab[k]['pass'], tab[k]['fail']
    t = p + f
    print(f"  {k:>8d}  {p:>14d}  {f:>14d}  {t:>8d}")

# Check perfect correlation
perfect_C1 = all(
    (r['kdim_F2'] == 1) == (not r['n2_div'])
    for r in results
)
print(f"\n  C1 perfect? {'YES ✓' if perfect_C1 else 'NO ✗'}")

# ── ANALYSIS 4: mod 4 kernel ─────────────────────────────────────────
print()
print("─" * 50)
print("ANALYSIS 4: 1 ∈ ker(A mod 4)?")
print("─" * 50)

tab4 = defaultdict(lambda: {'pass': 0, 'fail': 0})
for r in results:
    key = r['one_ker4']
    if r['n2_div']:
        tab4[key]['pass'] += 1
    else:
        tab4[key]['fail'] += 1

print(f"  {'1∈ker mod4':>12s}  {'n²|det (PASS)':>14s}  {'n²∤det (FAIL)':>14s}")
print(f"  {'─'*12}  {'─'*14}  {'─'*14}")
for k in sorted(tab4):
    p, f = tab4[k]['pass'], tab4[k]['fail']
    print(f"  {str(k):>12s}  {p:>14d}  {f:>14d}")

perfect_C3 = all(
    r['one_ker4'] == r['n2_div']
    for r in results
)
print(f"\n  C3 (1∈ker mod 4 ⟺ n²|det)? {'YES ✓' if perfect_C3 else 'NO ✗'}")

# ── ANALYSIS 5: v₂(det A) by kernel dimension ────────────────────────
print()
print("─" * 50)
print("ANALYSIS 5: v₂(det A) distribution by ker_F₂ dimension")
print("─" * 50)
for kdim in sorted(set(r['kdim_F2'] for r in results)):
    sub = [r for r in results if r['kdim_F2'] == kdim]
    v2s = Counter(r['v2_dA'] for r in sub)
    print(f"\n  dim ker_F₂ = {kdim} ({len(sub)} samples):")
    for v in sorted(v2s):
        print(f"    v₂(det A) = {v}: {v2s[v]}  ({100*v2s[v]/len(sub):.1f}%)")

print()
print("=" * 72)
print("SUMMARY")
print("=" * 72)
