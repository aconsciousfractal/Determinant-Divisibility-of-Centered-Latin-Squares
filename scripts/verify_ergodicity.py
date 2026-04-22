"""
verify_ergodicity.py — Critical test: is the swap-chain from cyclic LS biased?

HYPOTHESIS: Intercalate swaps from cyclic LS preserve the parity matrix
P_{ij} = L_{ij} mod 2, because all reachable intercalates swap same-parity values.
This would make ALL previous "always in kernel" results ARTIFACTS.

TEST: Generate LS using symbol permutations (isotopy) instead of swap chain.
If period-p vectors FAIL to be in kernel, the previous results were misleading.
"""

import numpy as np
from collections import Counter
import random, itertools

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

def bareiss(M):
    n = M.shape[0]
    A = M.astype(object)
    sign = 1
    prev = 1
    for i in range(n):
        if A[i, i] == 0:
            found = False
            for k in range(i + 1, n):
                if A[k, i] != 0:
                    A[[i, k]] = A[[k, i]]
                    sign *= -1
                    found = True
                    break
            if not found:
                return 0
        for j in range(i + 1, n):
            for k in range(i + 1, n):
                A[j, k] = (A[j, k] * A[i, i] - A[j, i] * A[i, k]) // prev
            A[j, i] = 0
        prev = A[i, i]
    return sign * A[n - 1, n - 1]

def get_A(L):
    n = L.shape[0]
    A = np.zeros((n - 1, n - 1), dtype=int)
    for i in range(n - 1):
        for j in range(n - 1):
            A[i, j] = L[i, j] - L[i, n - 1]
    return A

def cyclic_ls(n):
    L = np.zeros((n, n), dtype=int)
    for i in range(n):
        for j in range(n):
            L[i, j] = (i + j) % n + 1
    return L

def isotope_ls(L):
    """Random isotope: permute symbols, rows, and columns."""
    n = L.shape[0]
    # Symbol permutation
    sigma = list(range(1, n + 1))
    random.shuffle(sigma)
    L2 = np.zeros_like(L)
    for i in range(n):
        for j in range(n):
            L2[i, j] = sigma[L[i, j] - 1]
    # Row permutation
    row_perm = list(range(n))
    random.shuffle(row_perm)
    L3 = L2[row_perm, :]
    # Column permutation
    col_perm = list(range(n))
    random.shuffle(col_perm)
    L4 = L3[:, col_perm]
    return L4

def swap_chain_ls(n, swaps=2000):
    L = cyclic_ls(n)
    done = 0
    for _ in range(swaps * 20):
        if done >= swaps:
            break
        r1, r2 = random.sample(range(n), 2)
        c1, c2 = random.sample(range(n), 2)
        if L[r1, c1] == L[r2, c2] and L[r1, c2] == L[r2, c1]:
            L[r1, c1], L[r1, c2] = L[r1, c2], L[r1, c1]
            L[r2, c1], L[r2, c2] = L[r2, c2], L[r2, c1]
            done += 1
    return L

def check_parity_matrix(L):
    """Return parity matrix P_{ij} = L_{ij} mod 2."""
    return L % 2

# =====================================================================
print("=" * 72)
print("TEST 0: Verify parity invariance of swap chain from cyclic")
print("=" * 72)

n = 12
L_cyc = cyclic_ls(n)
P_cyc = check_parity_matrix(L_cyc)

print("Cyclic LS parity matrix (rows 0-3, cols 0-5):")
for i in range(4):
    print(f"  {list(P_cyc[i, :6])}")

# Check if parity changes after many swaps
parity_changed = 0
for trial in range(100):
    L = swap_chain_ls(n, swaps=5000)
    P = check_parity_matrix(L)
    if not np.array_equal(P, P_cyc):
        parity_changed += 1

print(f"\nParity matrix changed in {parity_changed}/100 swap-chain trials")
if parity_changed == 0:
    print("*** CONFIRMED: swap chain PRESERVES parity matrix! ***")
    print("*** ALL previous 'always in kernel' results are ARTIFACTS! ***")

# =====================================================================
print()
print("=" * 72)
print("TEST 1: Isotope LS of order 12 — period-3 vectors")
print("=" * 72)

n = 12
m = n - 1
N = 500

p3_vecs = []
for ph in range(3):
    p3_vecs.append(np.array([1 if j % 3 == ph else 0 for j in range(m)], dtype=int))
ones = np.ones(m, dtype=int)

p3_in_ker = [0, 0, 0]
ones_in_ker = 0
ker_dims = []

for trial in range(N):
    L = isotope_ls(cyclic_ls(n))
    A = get_A(L)
    B = A % 2
    kdim = kernel_dim_F2(A)
    ker_dims.append(kdim)
    
    if all((B @ ones) % 2 == 0):
        ones_in_ker += 1
    for ph in range(3):
        if all((B @ p3_vecs[ph]) % 2 == 0):
            p3_in_ker[ph] += 1

c = Counter(ker_dims)
print(f"Isotope LS (n=12), {N} samples:")
print(f"  ker dim dist: {dict(sorted(c.items()))}")
print(f"  min ker = {min(ker_dims)}")
print(f"  1-vec in ker: {ones_in_ker}/{N} ({100*ones_in_ker/N:.1f}%)")
for ph in range(3):
    print(f"  p3-{ph} in ker: {p3_in_ker[ph]}/{N} ({100*p3_in_ker[ph]/N:.1f}%)")

# =====================================================================
print()
print("=" * 72)
print("TEST 2: Isotope LS of order 8 — ker dim distribution")
print("=" * 72)

n = 8
m = n - 1
N = 500
ker_dims_8 = []
ones_8 = np.ones(m, dtype=int)
ones_in_ker_8 = 0

for trial in range(N):
    L = isotope_ls(cyclic_ls(n))
    A = get_A(L)
    kdim = kernel_dim_F2(A)
    ker_dims_8.append(kdim)
    B = A % 2
    if all((B @ ones_8) % 2 == 0):
        ones_in_ker_8 += 1

c = Counter(ker_dims_8)
print(f"Isotope LS (n=8), {N} samples:")
print(f"  ker dim dist: {dict(sorted(c.items()))}")
print(f"  min ker = {min(ker_dims_8)}")
print(f"  1-vec in ker: {ones_in_ker_8}/{N} ({100*ones_in_ker_8/N:.1f}%)")

# =====================================================================
print()
print("=" * 72)
print("TEST 3: Isotope LS of order 20 — period-5 vectors")
print("=" * 72)

n = 20
m = n - 1
N = 200

p5_vecs = [np.array([1 if j % 5 == ph else 0 for j in range(m)], dtype=int) 
           for ph in range(5)]
ones_20 = np.ones(m, dtype=int)
p5_in_ker = [0] * 5
ones_in_ker_20 = 0
ker_dims_20 = []

for trial in range(N):
    L = isotope_ls(cyclic_ls(n))
    A = get_A(L)
    B = A % 2
    kdim = kernel_dim_F2(A)
    ker_dims_20.append(kdim)
    if all((B @ ones_20) % 2 == 0):
        ones_in_ker_20 += 1
    for ph in range(5):
        if all((B @ p5_vecs[ph]) % 2 == 0):
            p5_in_ker[ph] += 1

c = Counter(ker_dims_20)
print(f"Isotope LS (n=20), {N} samples:")
print(f"  ker dim dist: {dict(sorted(c.items()))}")
print(f"  min ker = {min(ker_dims_20)}")
print(f"  1-vec in ker: {ones_in_ker_20}/{N} ({100*ones_in_ker_20/N:.1f}%)")
for ph in range(5):
    print(f"  p5-{ph} in ker: {p5_in_ker[ph]}/{N} ({100*p5_in_ker[ph]/N:.1f}%)")

# =====================================================================
print()
print("=" * 72)
print("TEST 4: Isotope LS of order 16 — ker dims")
print("=" * 72)

n = 16
m = n - 1
N = 200
ker_dims_16 = []
ones_16 = np.ones(m, dtype=int)
ones_in_ker_16 = 0

for trial in range(N):
    L = isotope_ls(cyclic_ls(n))
    A = get_A(L)
    kdim = kernel_dim_F2(A)
    ker_dims_16.append(kdim)
    B = A % 2
    if all((B @ ones_16) % 2 == 0):
        ones_in_ker_16 += 1

c = Counter(ker_dims_16)
print(f"Isotope LS (n=16), {N} samples:")
print(f"  ker dim dist: {dict(sorted(c.items()))}")
print(f"  min ker = {min(ker_dims_16)}")
print(f"  1-vec in ker: {ones_in_ker_16}/{N} ({100*ones_in_ker_16/N:.1f}%)")

# =====================================================================
print()
print("=" * 72)
print("TEST 5: Unified criterion check with isotope LS")
print("=" * 72)

print("n=8: checking dim ker(A²) ≥ 2 ↔ n²|det")
n = 8
m = n - 1
N = 300
match = 0
total = 0
fail_count = 0

for trial in range(N):
    L = isotope_ls(cyclic_ls(n))
    A = get_A(L)
    det_val = bareiss(A)
    B = A % 2
    B2 = (B @ B) % 2
    kdim2 = kernel_dim_F2(B2)
    
    criterion = (kdim2 >= 2)
    actual = (det_val % (n * n) == 0) if det_val != 0 else True
    
    if criterion == actual:
        match += 1
    total += 1
    if not actual:
        fail_count += 1

print(f"  Samples: {N}, matches: {match}/{total}, fail_rate: {fail_count}/{total}")
if match == total:
    print(f"  ★★★ UNIFIED CRITERION STILL PERFECT WITH ISOTOPE LS ★★★")
else:
    print(f"  *** MISMATCH FOUND! Criterion may be wrong! ***")

print()
print("n=12: checking dim ker(A²) ≥ 2 ↔ n²|det")
n = 12
m = n - 1
N = 200
match = 0
total = 0
fail_count = 0

for trial in range(N):
    L = isotope_ls(cyclic_ls(n))
    A = get_A(L)
    det_val = bareiss(A)
    B = A % 2
    B2 = (B @ B) % 2
    kdim2 = kernel_dim_F2(B2)
    
    criterion = (kdim2 >= 2)
    actual = (det_val % (n * n) == 0) if det_val != 0 else True
    
    if criterion == actual:
        match += 1
    total += 1
    if not actual:
        fail_count += 1

print(f"  Samples: {N}, matches: {match}/{total}, fail_rate: {fail_count}/{total}")
if match == total:
    print(f"  ★★★ UNIFIED CRITERION STILL PERFECT WITH ISOTOPE LS ★★★")
else:
    print(f"  *** MISMATCH FOUND! Criterion may be wrong! ***")

# =====================================================================
print()
print("=" * 72)
print("TEST 6: Multiple independent cyclic bases")
print("=" * 72)

def back_circulant_ls(n):
    """Back-circulant (different base LS)."""
    L = np.zeros((n, n), dtype=int)
    for i in range(n):
        for j in range(n):
            L[i, j] = (i * 2 + j) % n + 1 if n % 2 == 0 else (i + j) % n + 1
    # verify
    for i in range(n):
        if len(set(L[i])) != n:
            return None
    for j in range(n):
        if len(set(L[:, j])) != n:
            return None
    return L

def general_affine_ls(n, a, b):
    """L_{ij} = (a*i + j + b) mod n + 1, valid if gcd(a,n)=1."""
    import math
    if math.gcd(a, n) != 1:
        return None
    L = np.zeros((n, n), dtype=int)
    for i in range(n):
        for j in range(n):
            L[i, j] = (a * i + j + b) % n + 1
    return L

print("Testing different affine bases for n=12:")
n = 12
import math
for a in range(1, n):
    if math.gcd(a, n) != 1:
        continue
    L0 = general_affine_ls(n, a, 0)
    if L0 is None:
        continue
    # test with 50 isotopes
    ker_dims = []
    for _ in range(50):
        L = isotope_ls(L0)
        A = get_A(L)
        kdim = kernel_dim_F2(A)
        ker_dims.append(kdim)
    print(f"  a={a}: min ker={min(ker_dims)}, dist={dict(sorted(Counter(ker_dims).items()))}")
