"""
red_team_verify.py — Verify the reviewer's concerns:

1. Is "rank_p(A) = 1" actually proven or just observed?
   → Follows from A ≡ R (mod n) and R having identical rows.
   
2. For non-squarefree n (e.g. n=4,8,9): does "p | d_k" imply "n | d_k"?
   → Need rank over Z/p^a, not just Z/p. 
   → A ≡ R (mod p^a) for p^a | n, and R has rank 1 over Z/p^a.
   
3. Is the p-Sylow conjecture trivially true / reducible?
"""
import numpy as np
from math import gcd
from functools import reduce

def cyclic_ls(n):
    return np.array([[(i+j) % n + 1 for j in range(n)] for i in range(n)], dtype=int)

def get_A(L):
    n = L.shape[0]
    return np.array([[int(L[i,j]) - int(L[i,n-1]) for j in range(n-1)] for i in range(n-1)], dtype=int)

# ════════════════════════════════════════════════════════════════════
print("=" * 70)
print("  RED TEAM VERIFICATION")  
print("=" * 70)

# ─── POINT 1: Is A ≡ R (mod p^a) for p^a || n? ──────────────────
print("\n▶ POINT 1: A ≡ R (mod p^a) for every prime power p^a dividing n")
print("  This is the KEY claim. If true, rank(A mod p^a) = rank(R mod p^a) = 1.")

for n in [4, 8, 9, 12, 16, 25, 27, 36]:
    L = cyclic_ls(n)
    A = get_A(L)
    m = n - 1
    R = np.array([[j+1 for j in range(m)] for _ in range(m)], dtype=int)
    
    diff = A - R
    # Check: is every entry of (A - R) divisible by n?
    all_div_n = all(int(diff[i,j]) % n == 0 for i in range(m) for j in range(m))
    
    # Factor n and check divisibility by each p^a
    temp = n; factors = {}
    for p in range(2, n+1):
        while temp % p == 0:
            factors[p] = factors.get(p, 0) + 1
            temp //= p
    
    results = {}
    for p, a in factors.items():
        pa = p ** a
        all_div_pa = all(int(diff[i,j]) % pa == 0 for i in range(m) for j in range(m))
        results[f"p^a={pa}"] = all_div_pa
    
    tag = "✓" if all_div_n else "✗"
    print(f"  n={n:3d}: A ≡ R (mod {n}): {all_div_n} {tag}  |  per p^a: {results}")

print("""
  RESULT: A ≡ R (mod n) for ALL tested n. This is not coincidence —
  it follows algebraically from:
    A_ij = (i+j) mod n - (i+n-1) mod n
    ≡ (j+1) mod n = R_ij  (mod n)
  The congruence is mod n, not just mod p. So for every p^a | n,
  A ≡ R (mod p^a). QED.
""")

# ─── POINT 2: Does rank_1 over Z/p^a imply the right SNF bound? ──
print("▶ POINT 2: Rank of R over Z/p^a")
print("  R has all rows = [1,2,...,n-1]. Any 2×2 minor of R is zero")
print("  (identical rows!). Any 1×1 minor containing 1 is a unit mod p^a.")

for n in [4, 8, 9, 25, 27]:
    m = n - 1
    R = np.array([[j+1 for j in range(m)] for _ in range(m)])
    
    temp = n; factors = {}
    for p in range(2, n+1):
        while temp % p == 0:
            factors[p] = factors.get(p, 0) + 1
            temp //= p
    
    for p, a in factors.items():
        pa = p ** a
        # Check: gcd of first row entries with p^a
        row_gcd = reduce(gcd, [int(R[0,j]) % pa for j in range(m) if R[0,j] % pa != 0])
        # The entry R[0,0] = 1, which is a unit mod p^a for any p^a
        has_unit = (1 % pa != 0)  # always True since 1 is a unit
        
        # Any 2×2 minor of R
        det_22 = int(R[0,0])*int(R[1,1]) - int(R[0,1])*int(R[1,0])
        
        print(f"  n={n}, p^a={pa}: R[0,0]={R[0,0]} (unit mod {pa}? YES)  "
              f"2×2 minor = {det_22} (zero? {det_22==0})")

print("""
  RESULT: R mod p^a has rank exactly 1:
    - All rows identical → any 2×2 minor = 0 → rank ≤ 1  
    - Entry R[0,0] = 1 is a unit mod p^a → rank ≥ 1
  So rank(A mod p^a) = rank(R mod p^a) = 1.
  
  SNF IMPLICATION: #{k : v_p(d_k) < a} = rank over Z/p^a = 1.
  Since d_1 = 1 has v_p(d_1) = 0 < a, this is the only one.
  So v_p(d_k) ≥ a = v_p(n) for all k ≥ 2, i.e., n | d_k.  ✓
""")

# ─── POINT 3: The proof text in main.tex — is it sloppy? ─────────
print("▶ POINT 3: Is Step 3 in main.tex precisely correct?")
print("""
  Current text says: "A mod p has rank 1 for every prime p|n"
  Reviewer concern: this gives d_k divisible by p, not necessarily by n.
  
  For SQUAREFREE n: p | d_k for all p|n → lcm = n → n | d_k. FINE.
  For n = p^a (a ≥ 2): knowing "p | d_k" only gives d_k ≥ p, not d_k ≥ p^a.
  
  THE FIX: Replace "A mod p" with "A mod p^{v_p(n)}" in the proof text.
  The argument A ≡ R (mod p^{v_p(n)}) is equally valid (since A ≡ R mod n
  and p^{v_p(n)} | n). This gives v_p(d_k) ≥ v_p(n), hence n | d_k.
  
  VERDICT: The RESULT is correct. The PROOF TEXT has a minor imprecision
  for non-squarefree n. Easy 1-line fix.
""")

# ─── POINT 4: Is p-Sylow trivially implied? ──────────────────────
print("▶ POINT 4: Is the p-Sylow conjecture trivially true?")
print("  Claim: v_p(det A) ≥ v_p(n/p) + max(0, k_p - 1)")
print("  For p=2, this IS Theorem 5.5.")
print()

# Check: for odd p, is v_p(n/p) + max(0, k_p -1) ever STRONGER than
# just v_p(n) - 1 (which follows from Theorem 3.2)?
# Theorem 3.2 gives: n/gcd(n,2) | det(E_std), so n | det(E_std) for odd n.
# But det(E_std) = n * det(A), so n | n*det(A), giving... nothing for det(A).
# Actually Theorem 3.2 gives n/gcd(n,2) | n*det(A), i.e., 1/gcd(n,2) | det(A).
# For n odd: 1 | det(A). Not useful.
# For n even: 1/2 | det(A), i.e., det(A) can be half-integer? No, det(A) is integer.
# 
# Actually: the paper's Theorem 3.2 is about det(E_std), not det(A).
# det(E_std) = n * det(A). 
# Theorem 3.2: n^2/gcd(n,2) | det(E_std) = n*det(A)
# → n/gcd(n,2) | det(A)
#
# For odd p|n with v_p(n) = a:
# From Thm 3.2: v_p(det(A)) ≥ v_p(n/gcd(n,2)) = a (since gcd(n,2) is 1 or 2, coprime to odd p)
# p-Sylow claim: v_p(det(A)) ≥ v_p(n/p) + max(0, k_p-1) = (a-1) + max(0, k_p-1)
# 
# If k_p = 0: bound = a-1. But Thm 3.2 gives a. So Thm 3.2 is STRONGER.
# If k_p = 1: bound = a-1+0 = a-1. Thm 3.2 gives a. Still stronger.
# If k_p = 2: bound = a-1+1 = a. Same as Thm 3.2.
# If k_p = 3: bound = a-1+2 = a+1. NOW p-Sylow is stronger!
#
# So for odd primes: p-Sylow is only non-trivial when k_p ≥ 3.

print("  Comparing p-Sylow with Theorem 3.2's base bound:")
print("  Theorem 3.2 gives: v_p(det(A)) ≥ v_p(n) for odd p|n")
print("  p-Sylow gives:     v_p(det(A)) ≥ v_p(n/p) + max(0, k_p - 1) = v_p(n)-1 + max(0, k_p-1)")
print()
print("  k_p = 0: p-Sylow gives v_p(n)-1 < v_p(n) from Thm 3.2. WEAKER!")
print("  k_p = 1: p-Sylow gives v_p(n)-1 < v_p(n). WEAKER!")  
print("  k_p = 2: p-Sylow gives v_p(n). SAME as Thm 3.2.")
print("  k_p ≥ 3: p-Sylow gives v_p(n)-1+(k_p-1) > v_p(n). STRONGER!")
print()
print("  → For odd primes, the p-Sylow conjecture is WEAKER than Theorem 3.2")
print("    whenever k_p ≤ 2, which covers the vast majority of cases!")
print()

# Verify: how often is k_p ≥ 3 for odd primes?
import random
random.seed(42)

def jm_sample(n, steps=2000):
    L = list(range(n))
    grid = [L[i:] + L[:i] for i in range(n)]
    grid = [list(row) for row in grid]
    for _ in range(steps):
        i = random.randint(0, n-1)
        j = random.randint(0, n-1)
        s = random.randint(0, n-1)
        if grid[i][j] == s: continue
        i2 = None; j2 = None
        for ii in range(n):
            if ii != i and grid[ii][j] == s:
                i2 = ii; break
        for jj in range(n):
            if jj != j and grid[i][jj] == s:
                j2 = jj; break
        if i2 is None or j2 is None: continue
        t = grid[i2][j2]
        if t == grid[i][j]:
            grid[i][j] = s
            grid[i2][j] = grid[i][j2] = t
            grid[i2][j2] = grid[i][j]
            grid[i][j] = s; grid[i2][j2] = t
            grid[i2][j] = grid[i][j]; grid[i][j2] = grid[i2][j2]
            # Proper JM
            grid[i][j] = s
            grid[i][j2] = t  
            grid[i2][j] = t
            grid[i2][j2] = s
    return np.array([[v+1 for v in row] for row in grid])

def rank_Fp(M, p):
    n_rows, n_cols = M.shape
    B = (np.array(M, dtype=int) % p + p) % p
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
                factor = (int(B[row, col]) * inv_piv) % p
                for k in range(n_cols):
                    B[row, k] = (int(B[row, k]) - factor * int(B[pivot_row, k])) % p
        pivot_row += 1
    return pivot_row

print("  Checking k_p distribution for odd primes:")
for n, p in [(9, 3), (15, 3), (15, 5), (25, 5)]:
    kp_counts = {}
    for trial in range(100):
        L = jm_sample(n, steps=5000)
        A = get_A(L)
        r = rank_Fp(A, p)
        kp = (n-1) - r
        kp_counts[kp] = kp_counts.get(kp, 0) + 1
    print(f"  n={n}, p={p}: k_p distribution = {dict(sorted(kp_counts.items()))}")
    high = sum(v for k,v in kp_counts.items() if k >= 3)
    print(f"    k_p ≥ 3: {high}% → p-Sylow nontrivial in {high}% of cases")

# ─── POINT 5: Does Prop 3.6 cover even n? ────────────────────────
print("\n▶ POINT 5: Does Prop 3.6 cover even n?")
print("  Prop 3.6 is stated for ODD n. Step 2 says 'extended to even n by direct computation'.")
print("  Verify |det(A)| = n^{n-2} for even n:")

for n in [2, 4, 6, 8, 10]:
    L = cyclic_ls(n)
    A = get_A(L)
    # Use numpy for larger n
    det_val = int(round(np.linalg.det(A.astype(float))))
    expected = n ** (n-2)
    # For small n, verify exactly
    if n <= 8:
        from functools import reduce as fred
        # Bareiss for exact det
        m = n - 1
        M = [[int(A[i][j]) for j in range(m)] for i in range(m)]
        sign = 1; prev = 1
        for i in range(m):
            if M[i][i] == 0:
                for k in range(i+1, m):
                    if M[k][i] != 0:
                        M[i], M[k] = M[k], M[i]; sign *= -1; break
            for j in range(i+1, m):
                for k in range(i+1, m):
                    M[j][k] = (M[j][k]*M[i][i] - M[j][i]*M[i][k]) // prev
                M[j][i] = 0
            prev = M[i][i]
        det_exact = sign * M[m-1][m-1]
    else:
        det_exact = det_val
    
    match = (abs(det_exact) == expected)
    print(f"  n={n}: |det(A)| = {abs(det_exact)}, n^(n-2) = {expected}  {'✓' if match else '✗'}")


print("\n" + "=" * 70)
print("  REVIEWER ANALYSIS SUMMARY")
print("=" * 70)
print("""
  ┌────────────────────────────────────────────────────────────────────┐
  │ REVIEWER CLAIM                              │ VERDICT             │
  ├────────────────────────────────────────────────────────────────────┤
  │ Paper core is solid                         │ ✓ AGREE             │
  │ Cor 3.7 "rank mod p" is bottleneck          │ ⚠ PARTIALLY RIGHT  │
  │   → rank_p(A)=1 IS proven (A≡R mod n)       │   but proof TEXT    │
  │   → BUT text says "mod p" not "mod p^a"     │   needs 1-line fix  │
  │ C7 is strong science                        │ ✓ AGREE             │
  │ p-Sylow is premature (517 test insufficient)│ ✓ AGREE + WORSE:   │
  │   → for odd p, it's WEAKER than Thm 3.2!   │   partially trivial │
  │ n=6 enumeration is strong                   │ ✓ AGREE             │
  │ Missing deep "why"                          │ ✓ FAIR CRITICISM    │
  │ Overclaim on K(A)                           │ ⚠ PARTIALLY — we   │
  │                                             │   said "partial"    │
  └────────────────────────────────────────────────────────────────────┘
""")
