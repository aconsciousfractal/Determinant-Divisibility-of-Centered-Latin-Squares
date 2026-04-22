"""
explore_mixing.py — Test with aggressive mixing to find TRUE minimum ker dim.

Previous results:
  - explore_unify.py (unknown swaps): min ker n=12 = 3, n=20 = 5
  - explore_min_ker.py (500 swaps):   min ker n=12 = 6, n=20 = 16

The discrepancy means 500 swaps is insufficient mixing.
Need to dramatically increase swaps to find the true minimum.

Also: test whether specific vectors (1, period-3, alternating) are 
truly ALWAYS in ker(B) or just artifacts of weak mixing.
"""

import numpy as np
from collections import Counter
import random

def kernel_dim_F2(M_int):
    """dim ker over F_2 of an integer matrix (reduced mod 2)."""
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
    rank = pivot_row
    return cols - rank

def random_latin_square(n, swaps):
    """Generate random LS via cyclic base + many intercalate swaps."""
    L = np.zeros((n, n), dtype=int)
    for i in range(n):
        for j in range(n):
            L[i, j] = (i + j) % n + 1
    
    successful = 0
    attempts = 0
    max_attempts = swaps * 20  # allow many attempts to find valid swaps
    
    while successful < swaps and attempts < max_attempts:
        attempts += 1
        r1, r2 = random.sample(range(n), 2)
        c1, c2 = random.sample(range(n), 2)
        v1, v2 = L[r1, c1], L[r1, c2]
        v3, v4 = L[r2, c1], L[r2, c2]
        if v1 == v4 and v2 == v3:
            L[r1, c1], L[r1, c2] = v2, v1
            L[r2, c1], L[r2, c2] = v4, v3
            successful += 1
    
    return L, successful

def get_A_matrix(L):
    n = L.shape[0]
    A = np.zeros((n - 1, n - 1), dtype=int)
    for i in range(n - 1):
        for j in range(n - 1):
            A[i, j] = L[i, j] - L[i, n - 1]
    return A

def test_vector_in_ker(B, v):
    """Returns True if B·v ≡ 0 mod 2."""
    return all((B @ v) % 2 == 0)

# =========================================================================
print("=" * 72)
print("PART 1: Mixing quality test at n=12")
print("=" * 72)

n = 12
for num_swaps in [500, 2000, 10000, 50000]:
    ker_dims = []
    actual_swaps_list = []
    for trial in range(100):
        L, actual = random_latin_square(n, num_swaps)
        A = get_A_matrix(L)
        kdim = kernel_dim_F2(A)
        ker_dims.append(kdim)
        actual_swaps_list.append(actual)
    c = Counter(ker_dims)
    avg_actual = sum(actual_swaps_list) / len(actual_swaps_list)
    print(f"  swaps={num_swaps:>6d}, actual_avg={avg_actual:>8.0f}: "
          f"min ker={min(ker_dims)}, max={max(ker_dims)}, "
          f"dist={dict(sorted(c.items()))}")

# =========================================================================
print()
print("=" * 72)
print("PART 2: Mixing quality test at n=8")
print("=" * 72)

n = 8
for num_swaps in [200, 1000, 5000, 20000]:
    ker_dims = []
    actual_swaps_list = []
    for trial in range(200):
        L, actual = random_latin_square(n, num_swaps)
        A = get_A_matrix(L)
        kdim = kernel_dim_F2(A)
        ker_dims.append(kdim)
        actual_swaps_list.append(actual)
    c = Counter(ker_dims)
    avg_actual = sum(actual_swaps_list) / len(actual_swaps_list)
    print(f"  swaps={num_swaps:>6d}, actual_avg={avg_actual:>8.0f}: "
          f"min ker={min(ker_dims)}, max={max(ker_dims)}, "
          f"dist={dict(sorted(c.items()))}")

# =========================================================================
print()
print("=" * 72)
print("PART 3: With maximum mixing at n=12, test specific vectors")
print("=" * 72)

n = 12
m = n - 1  # = 11
ones = np.ones(m, dtype=int)
alt = np.array([1 if j % 2 == 0 else 0 for j in range(m)], dtype=int)
p3_0 = np.array([1 if j % 3 == 0 else 0 for j in range(m)], dtype=int)
p3_1 = np.array([1 if j % 3 == 1 else 0 for j in range(m)], dtype=int)
p3_2 = np.array([1 if j % 3 == 2 else 0 for j in range(m)], dtype=int)

vectors = {
    "1-vec": ones,
    "alt": alt,
    "p3-0": p3_0,
    "p3-1": p3_1,
    "p3-2": p3_2,
}

print("With 50000 swaps, 200 samples:")
results = {name: 0 for name in vectors}
for trial in range(200):
    L, _ = random_latin_square(n, 50000)
    A = get_A_matrix(L)
    B = A % 2
    for name, v in vectors.items():
        if test_vector_in_ker(B, v):
            results[name] += 1

for name, count in results.items():
    print(f"  {name}: in ker {count}/200 ({100*count/200:.0f}%)")

# =========================================================================
print()
print("=" * 72)
print("PART 4: n=20 mixing test")
print("=" * 72)

n = 20
for num_swaps in [1000, 5000, 20000, 100000]:
    ker_dims = []
    actual_swaps_list = []
    for trial in range(50):
        L, actual = random_latin_square(n, num_swaps)
        A = get_A_matrix(L)
        kdim = kernel_dim_F2(A)
        ker_dims.append(kdim)
        actual_swaps_list.append(actual)
    c = Counter(ker_dims)
    avg_actual = sum(actual_swaps_list) / len(actual_swaps_list)
    print(f"  swaps={num_swaps:>6d}, actual_avg={avg_actual:>8.0f}: "
          f"min ker={min(ker_dims)}, max={max(ker_dims)}, "
          f"dist={dict(sorted(c.items()))}")

# =========================================================================
print()
print("=" * 72)
print("PART 5: n=16 (pure power of 2) mixing test")
print("=" * 72)

n = 16
for num_swaps in [500, 5000, 50000]:
    ker_dims = []
    for trial in range(100):
        L, _ = random_latin_square(n, num_swaps)
        A = get_A_matrix(L)
        kdim = kernel_dim_F2(A)
        ker_dims.append(kdim)
    c = Counter(ker_dims)
    print(f"  swaps={num_swaps:>6d}: min ker={min(ker_dims)}, "
          f"max={max(ker_dims)}, dist={dict(sorted(c.items()))}")

# =========================================================================
print()
print("=" * 72)
print("PART 6: Theoretical analysis — 1 always in ker for n even")
print("=" * 72)

print("For ANY n×n Latin square L with even n:")
print("  A_{ij} = L_{ij} - L_{i,n}")
print("  row sum of A = sum_{j=0}^{n-2} A_{ij}")
print("       = sum_{j=0}^{n-2} L_{ij} - (n-1)·L_{i,n}")
print("       = [n(n+1)/2 - L_{i,n}] - (n-1)·L_{i,n}")
print("       = n(n+1)/2 - n·L_{i,n}")
print("       = n·[(n+1)/2 - L_{i,n}]")
print("  For n even: n·(...) is always even → row sums ≡ 0 mod 2")
print("  Therefore 1 ∈ ker(B) for ALL even n ✓")
print()
print("For n=12: need to check if period-3 vectors are PROVABLY in ker.")
print("  A_{ij} = L_{ij} - L_{i,12}")
print("  For p3-phase-0 vector v = (1,0,0,1,0,0,1,0,0,1,0):")
print("  (B·v)_i = sum of A_{i,j} mod 2 for j ∈ {0,3,6,9}")
print("          = [L_{i,0}+L_{i,3}+L_{i,6}+L_{i,9} - 4·L_{i,12}] mod 2")
print("  s = L_{i,0}+L_{i,3}+L_{i,6}+L_{i,9}: sum of 4 of 12 symbols in row i")
print("  This is NOT obviously always even. It depends on the LS.")
