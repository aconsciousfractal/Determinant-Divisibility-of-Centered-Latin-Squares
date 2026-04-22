"""
explore_structural.py — Focused structural analysis of kernel vectors.

CONFIRMED: min ker=3 at n=12 is stable across mixing levels (500-10000 swaps).
Now test: WHY? Are specific vectors PROVABLY in ker(B) for ALL LS?

Key insight from explore_min_ker.py (weak mixing):
  For n=12: 1-vector, alternating, and all 3 period-3 vectors were ALWAYS in ker.
  But that was with weak mixing. Need to retest with proper mixing.
"""

import numpy as np
from collections import Counter
import random

def kernel_dim_F2(M_int):
    B = np.array(M_int % 2, dtype=int)
    rows, cols = B.shape
    pivot_row = 0
    for col in range(cols):
        found = False
        for row in range(pivot_row, rows):
            if B[row, col] == 1:
                B[[pivot_row, row]] = B[[row, pivot_row]]
                found = True
                break
        if not found:
            continue
        for row in range(rows):
            if row != pivot_row and B[row, col] == 1:
                B[row] = (B[row] + B[pivot_row]) % 2
        pivot_row += 1
    return cols - pivot_row

def random_ls(n, target_swaps=2000):
    L = np.zeros((n, n), dtype=int)
    for i in range(n):
        for j in range(n):
            L[i, j] = (i + j) % n + 1
    done = 0
    for _ in range(target_swaps * 20):
        if done >= target_swaps:
            break
        r1, r2 = random.sample(range(n), 2)
        c1, c2 = random.sample(range(n), 2)
        if L[r1, c1] == L[r2, c2] and L[r1, c2] == L[r2, c1]:
            L[r1, c1], L[r1, c2] = L[r1, c2], L[r1, c1]
            L[r2, c1], L[r2, c2] = L[r2, c2], L[r2, c1]
            done += 1
    return L, done

def get_A(L):
    n = L.shape[0]
    A = np.zeros((n - 1, n - 1), dtype=int)
    for i in range(n - 1):
        for j in range(n - 1):
            A[i, j] = L[i, j] - L[i, n - 1]
    return A

# =================================================================
print("=" * 72)
print("TEST 1: n=12, vector-in-kernel test with strong mixing")
print("=" * 72)

n = 12
m = n - 1
N_SAMPLES = 500

vectors = {}
vectors["1-vec"] = np.ones(m, dtype=int)
vectors["alt(0,1)"] = np.array([1 if j % 2 == 0 else 0 for j in range(m)], dtype=int)
for ph in range(3):
    vectors[f"p3-{ph}"] = np.array([1 if j % 3 == ph else 0 for j in range(m)], dtype=int)

results = {name: 0 for name in vectors}
kdims = []

for trial in range(N_SAMPLES):
    L, actual = random_ls(n, target_swaps=3000)
    A = get_A(L)
    B = A % 2
    kdim = kernel_dim_F2(A)
    kdims.append(kdim)
    for name, v in vectors.items():
        if all((B @ v) % 2 == 0):
            results[name] += 1

c = Counter(kdims)
print(f"  Samples: {N_SAMPLES}, dim ker dist: {dict(sorted(c.items()))}")
print(f"  min ker = {min(kdims)}")
for name, count in results.items():
    pct = 100 * count / N_SAMPLES
    status = "ALWAYS" if count == N_SAMPLES else f"{pct:.1f}%"
    print(f"  {name}: in ker {count}/{N_SAMPLES} → {status}")

# =================================================================
print()
print("=" * 72)
print("TEST 2: n=8 (no odd factor > 1), same vectors")
print("=" * 72)

n = 8
m = n - 1
vectors8 = {}
vectors8["1-vec"] = np.ones(m, dtype=int)
vectors8["alt(0,1)"] = np.array([1 if j % 2 == 0 else 0 for j in range(m)], dtype=int)
# No period-3 makes sense for n=8 (n/gcd not aligned), but test anyway:
for ph in range(3):
    vectors8[f"p3-{ph}"] = np.array([1 if j % 3 == ph else 0 for j in range(m)], dtype=int)

results8 = {name: 0 for name in vectors8}
kdims8 = []

for trial in range(N_SAMPLES):
    L, _ = random_ls(n, target_swaps=2000)
    A = get_A(L)
    B = A % 2
    kdim = kernel_dim_F2(A)
    kdims8.append(kdim)
    for name, v in vectors8.items():
        if all((B @ v) % 2 == 0):
            results8[name] += 1

c = Counter(kdims8)
print(f"  Samples: {N_SAMPLES}, dim ker dist: {dict(sorted(c.items()))}")
print(f"  min ker = {min(kdims8)}")
for name, count in results8.items():
    pct = 100 * count / N_SAMPLES
    status = "ALWAYS" if count == N_SAMPLES else f"{pct:.1f}%"
    print(f"  {name}: in ker {count}/{N_SAMPLES} → {status}")

# =================================================================
print()
print("=" * 72)
print("TEST 3: n=20 (odd factor 5), period-5 vectors")
print("=" * 72)

n = 20
m = n - 1
vectors20 = {}
vectors20["1-vec"] = np.ones(m, dtype=int)
vectors20["alt(0,1)"] = np.array([1 if j % 2 == 0 else 0 for j in range(m)], dtype=int)
for ph in range(5):
    vectors20[f"p5-{ph}"] = np.array([1 if j % 5 == ph else 0 for j in range(m)], dtype=int)
# Also test period-3 (shouldn't work since 3∤20)
for ph in range(3):
    vectors20[f"p3-{ph}"] = np.array([1 if j % 3 == ph else 0 for j in range(m)], dtype=int)

results20 = {name: 0 for name in vectors20}
kdims20 = []

for trial in range(200):
    L, _ = random_ls(n, target_swaps=3000)
    A = get_A(L)
    B = A % 2
    kdim = kernel_dim_F2(A)
    kdims20.append(kdim)
    for name, v in vectors20.items():
        if all((B @ v) % 2 == 0):
            results20[name] += 1

c = Counter(kdims20)
print(f"  Samples: 200, dim ker dist: {dict(sorted(c.items()))}")
print(f"  min ker = {min(kdims20)}")
for name, count in results20.items():
    pct = 100 * count / 200
    status = "ALWAYS" if count == 200 else f"{pct:.1f}%"
    print(f"  {name}: in ker {count}/200 → {status}")

# =================================================================
print()
print("=" * 72)
print("TEST 4: Theoretical proof check for 1 and alt vectors")
print("=" * 72)

print("""
For B = A mod 2 where A_{ij} = L_{ij} - L_{i,n}, i,j ∈ {0,..,n-2}:

1-vector in ker(B):
  (B·1)_i = Σ_j B_{ij} = Σ_j A_{ij} mod 2 = [Σ_j L_{ij} - (n-1)·L_{i,n}] mod 2
  Σ_j=0^{n-2} L_{ij} = n(n+1)/2 - L_{i,n-1}  
  Wait: L uses 1..n, and j=0..n-2 covers columns 0..n-2 (n-1 columns).
  Row i of L: a permutation of 1..n. Sum = n(n+1)/2.
  Σ_{j=0}^{n-2} L_{ij} = n(n+1)/2 - L_{i,n-1}
  So: Σ_j A_{ij} = n(n+1)/2 - L_{i,n-1} - (n-1)·L_{i,n-1}
                  = n(n+1)/2 - n·L_{i,n-1}
                  = n·[(n+1)/2 - L_{i,n-1}]
  For even n: n is even, so n·(...) ≡ 0 mod 2. ✓
  
  Actually wait, L_{i,n-1} is the LAST column (column index n-1).
  Our A uses last column as base: A_{ij} = L_{ij} - L_{i,n-1}.
  The j index in A runs 0..n-2 (first n-1 columns of L).
  Row sum of A: Σ_{j=0}^{n-2} (L_{ij} - L_{i,n-1}) = [Σ L_{ij}] - (n-1)·L_{i,n-1}
                = [n(n+1)/2 - L_{i,n-1}] - (n-1)·L_{i,n-1}
                = n(n+1)/2 - n·L_{i,n-1}
  For n even: n(n+1)/2 = n/2·(n+1). Since n is even, n/2 is an integer.
  n(n+1)/2 mod 2: if n ≡ 0 mod 4, then n/2 is even, so n(n+1)/2 ≡ 0 mod 2.
                   if n ≡ 2 mod 4, then n/2 is odd, so n(n+1)/2 ≡ (n+1) mod 2 ≡ 1 mod 2.
  n·L_{i,n-1} mod 2: n even → always 0 mod 2.
  So row sum mod 2 = n(n+1)/2 mod 2:
    n ≡ 0 mod 4: row sum ≡ 0 mod 2 → 1 ∈ ker(B) ✓
    n ≡ 2 mod 4: row sum ≡ 1 mod 2 → 1 ∉ ker(B) ✗

This is why the n ≡ 2 mod 4 case is different!
""")

# =================================================================
print()
print("=" * 72)
print("TEST 5: Theoretical check — alternating vector")
print("=" * 72)

print("""
alt = (1,0,1,0,...) so (B·alt)_i = Σ_{j even, 0≤j≤n-2} A_{ij} mod 2.

For n ≡ 0 mod 4, n-1 is odd, so j runs 0,1,...,n-2.
Even j in {0,2,4,...,n-2}: there are n/2 such values.

(B·alt)_i = Σ_{j even} [L_{ij} - L_{i,n-1}] mod 2
          = [Σ_{j even} L_{ij}] - (n/2)·L_{i,n-1} mod 2

Σ_{j even} L_{ij}: sum of n/2 entries from row i of L (a permutation of 1..n).
We pick every other column. There's no obvious reason this is always even.

Unless there's a deeper structure from the Latin square property...
""")

# Let's verify numerically: gather the actual even-column sums
print("Numerical check at n=8:")
n_test = 8
for trial in range(10):
    L, _ = random_ls(n_test, target_swaps=2000)
    A = get_A(L)
    B = A % 2
    alt = np.array([1 if j % 2 == 0 else 0 for j in range(n_test - 1)], dtype=int)
    Balt = (B @ alt) % 2
    # Even-column sums
    even_sums = []
    for i in range(n_test - 1):
        s = sum(A[i, j] for j in range(n_test - 1) if j % 2 == 0)
        even_sums.append(s)
    print(f"  even col sums of A: {even_sums}, mod 2: {[s%2 for s in even_sums]}, "
          f"B·alt mod 2: {list(Balt)}")

print()
print("Numerical check at n=12:")
n_test = 12
for trial in range(10):
    L, _ = random_ls(n_test, target_swaps=3000)
    A = get_A(L)
    even_sums = []
    for i in range(n_test - 1):
        s = sum(A[i, j] for j in range(n_test - 1) if j % 2 == 0)
        even_sums.append(s)
    mods = [s % 2 for s in even_sums]
    print(f"  even col sums mod 2: {mods} (all 0? {all(m == 0 for m in mods)})")

# =================================================================
print()
print("=" * 72)
print("TEST 6: Period-p vector theoretical analysis for p|n")
print("=" * 72)

print("""
For p|n, consider the period-p vector v_p^{(r)} with v_j = 1 iff j ≡ r mod p.

(B·v)_i = Σ_{j ≡ r mod p, 0≤j≤n-2} A_{ij} mod 2
        = Σ_{j ≡ r mod p} [L_{ij} - L_{i,n-1}] mod 2

Number of such j in {0,...,n-2}: 
  If r < n-1 mod p or equal, the count depends on n-1 and p.
  For n = 12, p = 3: j ∈ {0,3,6,9} for r=0 (4 terms), {1,4,7,10} for r=1 (4 terms), 
                      {2,5,8} for r=2 (3 terms, since n-2=10 and 11 is next ≡2 mod 3).
  Wait: n-2 = 10. j runs 0..10.
  r=0: j ∈ {0,3,6,9} → 4 terms
  r=1: j ∈ {1,4,7,10} → 4 terms  
  r=2: j ∈ {2,5,8} → 3 terms

  So (B·v^{(r)})_i = [Σ L_{ij} for j≡r mod p] - count·L_{i,n-1} mod 2

For this to be 0 mod 2 for ALL L, we need:
  Σ_{j≡r mod p} L_{ij} ≡ count·L_{i,n-1} mod 2 for ALL LS.

This is a constraint on the sub-sums of each row of the LS.
""")

# Numerical: check sub-sums for n=12
print("n=12, p=3: sub-sums modulo 2 check:")
n_test = 12
p = 3
for trial in range(10):
    L, _ = random_ls(n_test, target_swaps=3000)  
    for r in range(p):
        js = [j for j in range(n_test - 1) if j % p == r]
        for i in range(min(3, n_test - 1)):
            s = sum(L[i, j] for j in js)
            base = L[i, n_test - 1]
            # (B·v)_i = [s - len(js)*base] mod 2
            val = (s - len(js) * base) % 2
            if trial == 0 and i == 0:
                print(f"  phase {r}: L entries at js={js}: {[L[i,j] for j in js]}, "
                      f"sum={s}, base={base}, count={len(js)}, "
                      f"s-count*base={s-len(js)*base}, mod 2 = {val}")

# Better: check statistics
print()
for r in range(p):
    js = [j for j in range(n_test - 1) if j % p == r]
    count_j = len(js)
    all_zero = 0
    total_checks = 0
    for trial in range(200):
        L, _ = random_ls(n_test, target_swaps=3000)
        all_rows_zero = True
        for i in range(n_test - 1):
            s = sum(int(L[i, j]) for j in js) - count_j * int(L[i, n_test - 1])
            if s % 2 != 0:
                all_rows_zero = False
                break
        total_checks += 1
        if all_rows_zero:
            all_zero += 1
    print(f"  p3-phase-{r}: {count_j} cols, all rows zero mod 2: {all_zero}/{total_checks} "
          f"({100*all_zero/total_checks:.0f}%)")

# =================================================================
print()
print("=" * 72)
print("TEST 7: n=20, period-5 sub-sums")
print("=" * 72)

n_test = 20
p = 5
for r in range(p):
    js = [j for j in range(n_test - 1) if j % p == r]
    count_j = len(js)
    all_zero = 0
    total_checks = 0
    for trial in range(100):
        L, _ = random_ls(n_test, target_swaps=3000)
        all_rows_zero = True
        for i in range(n_test - 1):
            s = sum(int(L[i, j]) for j in js) - count_j * int(L[i, n_test - 1])
            if s % 2 != 0:
                all_rows_zero = False
                break
        total_checks += 1
        if all_rows_zero:
            all_zero += 1
    print(f"  p5-phase-{r}: {count_j} cols, all rows zero mod 2: {all_zero}/{total_checks} "
          f"({100*all_zero/total_checks:.0f}%)")
