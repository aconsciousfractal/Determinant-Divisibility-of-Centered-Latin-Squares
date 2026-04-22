"""
verify_base_column.py — Diagnose whether criterion failure is due to base column change.

Key question: does n²|det(A) depend on which column is chosen as base?
If yes, the isotopy test with column permutation was unfair.

Test:
1. For a single LS, compute A with ALL possible base columns
2. Check if n²|det varies across bases
3. Test criterion with isotopy WITHOUT column permutation
4. Test criterion with isotopy WITH column permutation
"""
import sys as _sys
try:
    _sys.stdout.reconfigure(encoding="utf-8")
    _sys.stderr.reconfigure(encoding="utf-8")
except Exception:
    pass

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

def get_A_base(L, base_col):
    """A matrix using specified column as base, dropping last row."""
    n = L.shape[0]
    cols = [j for j in range(n) if j != base_col]
    A = np.zeros((n - 1, n - 1), dtype=int)
    for i in range(n - 1):
        for idx, j in enumerate(cols):
            A[i, idx] = L[i, j] - L[i, base_col]
    return A

def cyclic_ls(n):
    L = np.zeros((n, n), dtype=int)
    for i in range(n):
        for j in range(n):
            L[i, j] = (i + j) % n + 1
    return L

def swap_chain_ls(n, swaps=3000):
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

# =====================================================================
print("=" * 72)
print("TEST 1: Base column dependence for a single LS")
print("=" * 72)

n = 8
# Generate a swap chain LS
L = swap_chain_ls(n, swaps=2000)
print(f"LS (n={n}):")
for i in range(n):
    print(f"  {list(L[i])}")

print(f"\nBase column scan (n²={n*n}):")
for base in range(n):
    A = get_A_base(L, base)
    det_val = bareiss(A)
    div = (det_val % (n * n) == 0) if det_val != 0 else True
    B = A % 2
    B2 = (B @ B) % 2
    kdim = kernel_dim_F2(A)
    kdim2 = kernel_dim_F2(B2)
    print(f"  base={base}: det={det_val:>12d}, {n}²|det={div}, "
          f"kdim={kdim}, kdim(B²)={kdim2}")

# =====================================================================
print()
print("=" * 72)
print("TEST 2: Same test, different LS")
print("=" * 72)

for trial in range(3):
    L = swap_chain_ls(n, swaps=2000)
    all_div = []
    for base in range(n):
        A = get_A_base(L, base)
        det_val = bareiss(A)
        div = (det_val % (n * n) == 0) if det_val != 0 else True
        all_div.append(div)
    n_pass = sum(all_div)
    print(f"  Trial {trial}: {n_pass}/{n} bases pass n²|det → {all_div}")

# =====================================================================
print()
print("=" * 72)
print("TEST 3: Isotopy WITHOUT column permutation (symbol + row only)")
print("=" * 72)

def isotope_no_col(L):
    """Symbol permutation + row permutation only. Base column stays n-1."""
    n = L.shape[0]
    sigma = list(range(1, n + 1))
    random.shuffle(sigma)
    L2 = np.zeros_like(L)
    for i in range(n):
        for j in range(n):
            L2[i, j] = sigma[L[i, j] - 1]
    row_perm = list(range(n))
    random.shuffle(row_perm)
    L3 = L2[row_perm, :]
    return L3

n = 8
N = 500
match = 0
total = 0
mismatches = []

for trial in range(N):
    L = isotope_no_col(cyclic_ls(n))
    A = get_A_base(L, n - 1)  # always use last column as base
    det_val = bareiss(A)
    B = A % 2
    B2 = (B @ B) % 2
    kdim2 = kernel_dim_F2(B2)
    
    criterion = (kdim2 >= 2)
    actual = (det_val % (n * n) == 0) if det_val != 0 else True
    
    if criterion == actual:
        match += 1
    else:
        if len(mismatches) < 5:
            mismatches.append((det_val, kdim2, criterion, actual))
    total += 1

print(f"n={n}: isotope (no col perm), {N} samples: {match}/{total} match")
if match < total:
    print("  Mismatches:")
    for det_val, k2, crit, act in mismatches:
        dir_str = "crit=PASS,actual=FAIL" if crit and not act else "crit=FAIL,actual=PASS"
        print(f"    det={det_val}, kdim(B²)={k2}, {dir_str}")

# =====================================================================
print()
print("Same for n=12:")
n = 12
N = 300
match = 0
total = 0
mismatches = []

for trial in range(N):
    L = isotope_no_col(cyclic_ls(n))
    A = get_A_base(L, n - 1)
    det_val = bareiss(A)
    B = A % 2
    B2 = (B @ B) % 2
    kdim2 = kernel_dim_F2(B2)
    
    criterion = (kdim2 >= 2)
    actual = (det_val % (n * n) == 0) if det_val != 0 else True
    
    if criterion == actual:
        match += 1
    else:
        if len(mismatches) < 5:
            mismatches.append((det_val, kdim2, criterion, actual))
    total += 1

print(f"n={n}: isotope (no col perm), {N} samples: {match}/{total} match")
if match < total:
    print("  Mismatches:")
    for det_val, k2, crit, act in mismatches:
        dir_str = "crit=PASS,actual=FAIL" if crit and not act else "crit=FAIL,actual=PASS"
        print(f"    det={det_val}, kdim(B²)={k2}, {dir_str}")

# =====================================================================
print()
print("=" * 72)
print("TEST 4: Isotopy WITH column permutation — same")
print("=" * 72)

def isotope_full(L):
    n = L.shape[0]
    sigma = list(range(1, n + 1))
    random.shuffle(sigma)
    L2 = np.zeros_like(L)
    for i in range(n):
        for j in range(n):
            L2[i, j] = sigma[L[i, j] - 1]
    row_perm = list(range(n))
    random.shuffle(row_perm)
    L3 = L2[row_perm, :]
    col_perm = list(range(n))
    random.shuffle(col_perm)
    L4 = L3[:, col_perm]
    return L4

n = 8
N = 500
match_full = 0
mismatches_full = []

for trial in range(N):
    L = isotope_full(cyclic_ls(n))
    A = get_A_base(L, n - 1)
    det_val = bareiss(A)
    B = A % 2
    B2 = (B @ B) % 2
    kdim2 = kernel_dim_F2(B2)
    
    criterion = (kdim2 >= 2)
    actual = (det_val % (n * n) == 0) if det_val != 0 else True
    
    if criterion == actual:
        match_full += 1
    else:
        if len(mismatches_full) < 5:
            mismatches_full.append((det_val, kdim2, criterion, actual))

print(f"n={n}: isotope (FULL), {N} samples: {match_full}/{N} match")
if match_full < N:
    print("  Mismatches:")
    for det_val, k2, crit, act in mismatches_full:
        dir_str = "crit=PASS,actual=FAIL" if crit and not act else "crit=FAIL,actual=PASS"
        print(f"    det={det_val}, kdim(B²)={k2}, {dir_str}")

# =====================================================================
print()
print("=" * 72)
print("TEST 5: Swap chain — reconfirm criterion")
print("=" * 72)

for n in [8, 12]:
    N = 300
    match_sc = 0
    for trial in range(N):
        L = swap_chain_ls(n, swaps=3000)
        A = get_A_base(L, n - 1)
        det_val = bareiss(A)
        B = A % 2
        B2 = (B @ B) % 2
        kdim2 = kernel_dim_F2(B2)
        criterion = (kdim2 >= 2)
        actual = (det_val % (n * n) == 0) if det_val != 0 else True
        if criterion == actual:
            match_sc += 1
    print(f"n={n}: swap chain, {N} samples: {match_sc}/{N} match")

# =====================================================================
print()
print("=" * 72)
print("TEST 6: ker dim distribution comparison")
print("=" * 72)

for n in [8, 12]:
    print(f"\nn={n}:")
    # Swap chain
    dims_sc = []
    for _ in range(300):
        L = swap_chain_ls(n, swaps=3000)
        A = get_A_base(L, n - 1)
        dims_sc.append(kernel_dim_F2(A))
    c_sc = Counter(dims_sc)
    print(f"  Swap chain: min={min(dims_sc)}, dist={dict(sorted(c_sc.items()))}")
    
    # Isotope no col
    dims_iso = []
    for _ in range(300):
        L = isotope_no_col(cyclic_ls(n))
        A = get_A_base(L, n - 1)
        dims_iso.append(kernel_dim_F2(A))
    c_iso = Counter(dims_iso)
    print(f"  Iso(no col): min={min(dims_iso)}, dist={dict(sorted(c_iso.items()))}")
    
    # Isotope full
    dims_full = []
    for _ in range(300):
        L = isotope_full(cyclic_ls(n))
        A = get_A_base(L, n - 1)
        dims_full.append(kernel_dim_F2(A))
    c_full = Counter(dims_full)
    print(f"  Iso(full):   min={min(dims_full)}, dist={dict(sorted(c_full.items()))}")
