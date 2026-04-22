"""
Deep investigation of dim ker = 1 at n=8: what separates PASS from FAIL?

Key facts:
  - dim ker = 1 means ker_F₂(A mod 2) = span{1}
  - ALL fail: v₂(det A) = 2 exactly (|det A| ≡ 4 mod 8)
  - ALL pass: v₂(det A) ≥ 3 (|det A| ≡ 0 mod 8)
  - 77% pass, 23% fail

Hypotheses to test:
  H1: The SECOND elementary divisor of A is odd for FAIL, even for PASS
  H2: Some row-structure invariant (e.g., parity pattern of row pairs)
  H3: det(A mod 4) analysis — can we predict from A mod 4?
  H4: Rank of A² mod 2 (higher-order structure)
  H5: Number of even entries in A (or some combinatorial invariant of L)
  H6: Kernel dimension of (A^T A) mod 2 or other derived matrices
  H7: Structural property of the Latin square itself (e.g., transversal count)
"""
import numpy as np
from collections import Counter, defaultdict
from math import gcd

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

def smith_diagonal(M):
    """Compute Smith Normal Form diagonal entries via iterative Euclidean reduction."""
    m, n = M.shape
    B = [[int(M[i,j]) for j in range(n)] for i in range(m)]
    size = min(m, n)
    
    for s in range(size):
        # Step 1: find non-zero entry
        found = False
        for i in range(s, m):
            for j in range(s, n):
                if B[i][j] != 0:
                    B[s], B[i] = B[i], B[s]
                    for r in range(m):
                        B[r][s], B[r][j] = B[r][j], B[r][s]
                    found = True
                    break
            if found: break
        if not found: break
        
        # Step 2: iteratively reduce
        for _ in range(200):  # safety limit
            changed = False
            # Column reduction
            for i in range(s+1, m):
                if B[i][s] != 0:
                    q = B[i][s] // B[s][s]
                    for j in range(n):
                        B[i][j] -= q * B[s][j]
                    if B[i][s] != 0:
                        # |B[i][s]| < |B[s][s]| after reduction, swap
                        B[s], B[i] = B[i], B[s]
                        changed = True
            # Row reduction  
            for j in range(s+1, n):
                if B[s][j] != 0:
                    q = B[s][j] // B[s][s]
                    for i in range(m):
                        B[i][j] -= q * B[i][s]
                    if B[s][j] != 0:
                        for i in range(m):
                            B[i][s], B[i][j] = B[i][j], B[i][s]
                        changed = True
            # Check divisibility of remaining submatrix
            if not changed:
                all_div = True
                for i in range(s+1, m):
                    for j in range(s+1, n):
                        if B[i][j] % B[s][s] != 0:
                            # Add row i to row s, continue
                            for jj in range(n):
                                B[s][jj] += B[i][jj]
                            all_div = False
                            break
                    if not all_div: break
                if all_div:
                    break
    
    return [abs(B[i][i]) if i < m and i < n and B[i][i] != 0 else 0 for i in range(size)]

# ═══════════════════════════════════════════════════════════════════════
n = 8
print("=" * 72)
print(f"DEEP INVESTIGATION: dim ker_F₂ = 1 at n={n}")
print("=" * 72)

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

# Collect dim-ker-1 samples with full detail
dim1_pass = []
dim1_fail = []
other = 0
total = 0

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
        
        kdim = kernel_dim_F2(A)
        if kdim != 1:
            other += 1
            continue
        
        n2_div = (dE % (n**2) == 0)
        v2_dA = v_p(dA, 2)
        
        # Compute features
        A_mod2 = A % 2
        A_mod4 = A % 4
        
        # H3: det(A mod 4) — but det over Z/4Z isn't well-defined for non-prime
        # Instead: count entries divisible by 2 in A
        even_count = int(np.sum(A % 2 == 0))
        
        # H4: rank of A² mod 2
        A2 = A @ A
        kdim_A2 = kernel_dim_F2(A2)
        
        # H5: trace of A mod 2
        trace_mod2 = int(np.trace(A)) % 2
        
        # H6: (A^T A) mod 2 kernel
        ATA = A.T @ A
        kdim_ATA = kernel_dim_F2(ATA)
        
        # H7: Parity of number of even entries per row
        even_per_row = tuple(int(np.sum(A[i,:] % 2 == 0)) for i in range(n-1))
        total_even_per_row = sum(even_per_row)
        
        # H8: A mod 4 — compute (A/2 mod 2) where A entries are even
        # i.e., for even entries of A, how many are ≡ 2 mod 4 vs ≡ 0 mod 4?
        entries_0mod4 = int(np.sum((A % 4 == 0) & (A % 2 == 0)))
        entries_2mod4 = int(np.sum(np.abs(A) % 4 == 2))
        
        # H9: det(A/gcd) reduced — look at A after factoring out row GCDs
        row_gcds = []
        for i in range(n-1):
            g = 0
            for j in range(n-1):
                g = gcd(g, abs(int(A[i,j])))
            row_gcds.append(g)
        
        # H10: sum of entries mod 4
        sum_mod4 = int(np.sum(A)) % 4
        
        # H11: permanent-like invariant — parity of column sums
        col_sum_mod2 = tuple(int(np.sum(A[:,j])) % 2 for j in range(n-1))
        
        record = {
            'v2_dA': v2_dA,
            'n2_div': n2_div,
            'even_count': even_count,
            'kdim_A2': kdim_A2,
            'trace_mod2': trace_mod2,
            'kdim_ATA': kdim_ATA,
            'even_per_row': even_per_row,
            'total_even': total_even_per_row,
            'entries_0mod4': entries_0mod4,
            'entries_2mod4': entries_2mod4,
            'row_gcds': tuple(row_gcds),
            'sum_mod4': sum_mod4,
            'col_sum_mod2': col_sum_mod2,
            'A': A.copy(),
            'L': L.copy(),
        }
        
        if n2_div:
            dim1_pass.append(record)
        else:
            dim1_fail.append(record)

print(f"Total: {total}, dim ker=1: {len(dim1_pass)+len(dim1_fail)}, dim ker≥2: {other}")
print(f"dim ker=1: {len(dim1_pass)} PASS, {len(dim1_fail)} FAIL ({100*len(dim1_fail)/(len(dim1_pass)+len(dim1_fail)):.1f}% fail)")
print()

# ── Test each hypothesis ──────────────────────────────────────────────

def test_binary_hypothesis(name, key, pass_data, fail_data):
    """Test if a feature perfectly separates PASS from FAIL."""
    pass_vals = Counter(r[key] for r in pass_data)
    fail_vals = Counter(r[key] for r in fail_data)
    all_vals = sorted(set(list(pass_vals.keys()) + list(fail_vals.keys())))
    
    print(f"\n{'─'*50}")
    print(f"  {name}")
    print(f"{'─'*50}")
    print(f"  {'Value':>20s}  {'PASS':>8s}  {'FAIL':>8s}")
    print(f"  {'─'*20}  {'─'*8}  {'─'*8}")
    for v in all_vals:
        p, f = pass_vals.get(v, 0), fail_vals.get(v, 0)
        marker = " ←" if (p == 0 or f == 0) and (p + f > 5) else ""
        print(f"  {str(v):>20s}  {p:>8d}  {f:>8d}{marker}")
    
    # Check if perfectly separable by threshold
    overlap = sum(min(pass_vals.get(v, 0), fail_vals.get(v, 0)) for v in all_vals)
    total_samples = sum(pass_vals.values()) + sum(fail_vals.values())
    print(f"  Overlap: {overlap}/{total_samples} ({100*overlap/total_samples:.1f}%)")

test_binary_hypothesis("H3: even entries in A", 'even_count', dim1_pass, dim1_fail)
test_binary_hypothesis("H4: dim ker_F₂(A²)", 'kdim_A2', dim1_pass, dim1_fail)
test_binary_hypothesis("H5: trace(A) mod 2", 'trace_mod2', dim1_pass, dim1_fail)
test_binary_hypothesis("H6: dim ker_F₂(AᵀA)", 'kdim_ATA', dim1_pass, dim1_fail)
test_binary_hypothesis("H7: total even entries per row", 'total_even', dim1_pass, dim1_fail)
test_binary_hypothesis("H8a: entries ≡ 0 mod 4", 'entries_0mod4', dim1_pass, dim1_fail)
test_binary_hypothesis("H8b: entries with |a| ≡ 2 mod 4", 'entries_2mod4', dim1_pass, dim1_fail)
test_binary_hypothesis("H10: sum(A) mod 4", 'sum_mod4', dim1_pass, dim1_fail)
test_binary_hypothesis("H11: col sum parity pattern", 'col_sum_mod2', dim1_pass, dim1_fail)

# ── SMITH NORMAL FORM analysis ────────────────────────────────────────
print()
print("=" * 72)
print("SMITH NORMAL FORM — 2-adic structure of elementary divisors")
print("=" * 72)

# Sample some PASS and FAIL for SNF
snf_pass = []
snf_fail = []
for r in dim1_pass[:100]:
    d = smith_diagonal(r['A'])
    v2s = tuple(v_p(di, 2) for di in d if di > 0)
    snf_pass.append(v2s)
for r in dim1_fail[:100]:
    d = smith_diagonal(r['A'])
    v2s = tuple(v_p(di, 2) for di in d if di > 0)
    snf_fail.append(v2s)

print(f"\nSample PASS SNF v₂ profiles (first 10):")
for s in snf_pass[:10]:
    print(f"  {s}  (sum={sum(s)})")
print(f"\nSample FAIL SNF v₂ profiles (first 10):")
for s in snf_fail[:10]:
    print(f"  {s}  (sum={sum(s)})")

# Aggregate: what's different about the 2-adic elementary divisors?
print(f"\nMin entry of v₂ profile (2nd smallest if first always 0):")
pass_min_v2_1 = Counter(s[0] if len(s) > 0 else -1 for s in snf_pass)
fail_min_v2_1 = Counter(s[0] if len(s) > 0 else -1 for s in snf_fail)
print(f"  PASS — smallest v₂(d_i): {dict(sorted(pass_min_v2_1.items()))}")
print(f"  FAIL — smallest v₂(d_i): {dict(sorted(fail_min_v2_1.items()))}")

# Count how many elementary divisors are odd
pass_odd_count = Counter(sum(1 for v in s if v == 0) for s in snf_pass)
fail_odd_count = Counter(sum(1 for v in s if v == 0) for s in snf_fail)
print(f"\n  # odd elementary divisors — PASS: {dict(sorted(pass_odd_count.items()))}")
print(f"  # odd elementary divisors — FAIL: {dict(sorted(fail_odd_count.items()))}")

# Count how many elem divisors have v₂ ≥ 2
pass_v2ge2 = Counter(sum(1 for v in s if v >= 2) for s in snf_pass)
fail_v2ge2 = Counter(sum(1 for v in s if v >= 2) for s in snf_fail)
print(f"\n  # elem divisors with v₂ ≥ 2 — PASS: {dict(sorted(pass_v2ge2.items()))}")
print(f"  # elem divisors with v₂ ≥ 2 — FAIL: {dict(sorted(fail_v2ge2.items()))}")

# Sum of v₂ profile
pass_sum = Counter(sum(s) for s in snf_pass)
fail_sum = Counter(sum(s) for s in snf_fail)
print(f"\n  sum(v₂ profile) — PASS min: {min(pass_sum.keys())}, max: {max(pass_sum.keys())}")
print(f"  sum(v₂ profile) — FAIL min: {min(fail_sum.keys())}, max: {max(fail_sum.keys())}")

# ── SPECIAL: Look at the half-matrix A/2 mod 2 ───────────────────────
print()
print("=" * 72)
print("HALF-MATRIX ANALYSIS: use (A//2) mod 2 where entries are even")
print("=" * 72)
print("(For dim ker=1 at n=8, every entry of A is even due to row structure)")

# Check: at n=8, A_{ij} = L_{ij} - L_{i,8}. Range is -7..7.
# Row sums = (n/2)(odd) = 4·odd. So row sums are ≡ 0 mod 4.
# Are individual entries always even? No — they can be odd.
# But A·1 ≡ 0 mod 2, and if dim ker = 1, then 1 is the ONLY kernel vector.

# Let's look at A/2 where possible
print("\nEntry parity analysis:")
pass_all_even = 0
fail_all_even = 0
for r in dim1_pass[:200]:
    if np.all(r['A'] % 2 == 0):
        pass_all_even += 1
for r in dim1_fail[:200]:
    if np.all(r['A'] % 2 == 0):
        fail_all_even += 1
n_pass_sample = min(200, len(dim1_pass))
n_fail_sample = min(200, len(dim1_fail))
print(f"  PASS: {pass_all_even}/{n_pass_sample} have ALL entries even")
print(f"  FAIL: {fail_all_even}/{n_fail_sample} have ALL entries even")
print("  (A entries need NOT all be even — they're L_{ij} - L_{i,n})")

# ── HALF-DETERMINANT: for even A, compute det(A/2) mod 2 ─────────────
# If A has some odd entries, we can't simply halve. 
# Instead, let B = A mod 4 and look at the structure of B/2 for even entries.

# Actually, let's try a different approach:
# The KEY QUESTION: among dim-ker-1 LS, what makes det(A) ≡ 0 mod 8 vs ≡ 4 mod 8?
# Equivalently: is v₂(det A) ≥ 3 or exactly 2?

# The Hensel-lift approach: reduce A mod 4. Write A = 2C + R where R ∈ {0,1}^{7×7} 
# (the mod-2 reduction) and C captures the "half" of even entries.
# det(A) = det(2C + R) and we need this mod 8.

print()
print("=" * 72)  
print("HENSEL APPROACH: det(A) mod 8 from A mod 4")
print("=" * 72)
print("Write A = 2·Q + R where R = A mod 2, Q = (A - R)/2 = (A mod 4 - R)/2 'ish'")
print("Actually: just compute det(A mod 8) and see if it's 0 or 4.")
print()

# For a cleaner test: compute det(A) mod 8 using the cofactor expansion
# This is just bareiss mod 8... but Bareiss over Z gives exact det, which we reduce.
# Already have det(A). Just look at det(A) mod 8.

# More interesting: look at det(A/gcd(row)) properties
# For dim-ker-1 at n=8, all row sums are 4·odd. 
# The rows themselves have sum = 4·(n+1-2L_{i,n-1})/2... wait, let me recompute.
# (A·1)_i = sum_j (L_{ij} - L_{i,n}) = sum_j L_{ij} - (n-1)·L_{i,n}
# = (1+2+...+n) - n·L_{i,n} = n(n+1)/2 - n·L_{i,n} = n·((n+1)/2 - L_{i,n})
# Wait, sum_j goes j=0..n-2 (n-1 terms), and we subtract L_{i,n-1} from each.
# Actually A_{ij} = L_{i,j} - L_{i,n-1} for j=0..n-2. So:
# (A·1)_i = sum_{j=0}^{n-2} (L_{i,j} - L_{i,n-1}) = sum_{j=0}^{n-2} L_{i,j} - (n-1)·L_{i,n-1}
# = (sum all - L_{i,n-1}) - (n-1)·L_{i,n-1} = n(n+1)/2 - n·L_{i,n-1}
# = n·( (n+1)/2 - L_{i,n-1} )
# For n=8: = 8·(4.5 - L_{i,7}). But L_{i,7} is integer, so this is 8·(4.5 - k) = 4·(9-2k).
# 9-2k is always odd. So (A·1)_i = 4·(odd). And v₂((A·1)_i) = 2 exactly.

# So all row sums have EXACTLY v₂ = 2. The "extra" 2-adic power in det(A) comes from
# the LINEAR DEPENDENCE structure, not from individual row sums.

# NEW IDEA: Look at A/4 · (1...1)^T — since A·1 = 4·u (u odd), 
# the matrix B = A (mod 4) might capture the relevant structure.

# Let's try: compute rank_F₂(A // 2) — the matrix where we integer-divide by 2
# This is NOT well-defined for odd entries... but A's entries can be odd.
# Skip this.

# Instead: look at the ADJUGATE matrix approach.
# det(A) = sum of products of entries. We need those products mod 8.

# SIMPLEST WORKING TEST: 
# At dim ker_F₂ = 1, the rank mod 2 is n-2 = 6. 
# One SNF elementary divisor is even (the one from 1 ∈ ker).
# QUESTION: Is this even divisor ≡ 2 mod 4 (FAIL) or ≡ 0 mod 4 (PASS)?

print("QUESTION: In SNF, which elementary divisor is 'the even one' from ker?")
print("And is it ≡ 2 mod 4 (→ v₂=1 contribution → total v₂=2 → FAIL)")
print("or ≡ 0 mod 4 (→ v₂≥2 contribution → total v₂≥3 → PASS)?")
print()

# Look at sorted SNF profiles
print("SNF v₂ profiles (sorted ascending):")
print("\nFAIL examples:")
for s in snf_fail[:15]:
    ss = sorted(s)
    print(f"  {ss}  product v₂ = {sum(ss)}")
print("\nPASS examples:")
for s in snf_pass[:15]:
    ss = sorted(s)
    print(f"  {ss}  product v₂ = {sum(ss)}")

# Count: for FAIL, exactly how many SNF entries have v₂=1?
print("\nNumber of SNF entries with v₂ = 1:")
fail_v2eq1 = Counter(sum(1 for v in s if v == 1) for s in snf_fail)
pass_v2eq1 = Counter(sum(1 for v in s if v == 1) for s in snf_pass)
print(f"  FAIL: {dict(sorted(fail_v2eq1.items()))}")
print(f"  PASS: {dict(sorted(pass_v2eq1.items()))}")

# More refined: for FAIL, the SNF should have exactly ONE entry with v₂=1 
# and the rest with v₂=0, giving total v₂ = 1 + base.
# For PASS, the "even" entry has v₂ ≥ 2.

# DEFINITIVE: look at the second-largest v₂ in the SNF
print("\nSecond-largest v₂ in SNF:")
def second_max(lst):
    s = sorted(lst, reverse=True)
    return s[1] if len(s) > 1 else 0
pass_2nd = Counter(second_max(s) for s in snf_pass)
fail_2nd = Counter(second_max(s) for s in snf_fail)
print(f"  PASS: {dict(sorted(pass_2nd.items()))}")
print(f"  FAIL: {dict(sorted(fail_2nd.items()))}")
