"""
verify_definitive.py — Final definitive test of the criterion.

The "unified criterion" dim ker(B²) ≥ 2 ↔ n²|det appeared perfect with
swap chain from cyclic (few swaps). But isotopy breaks it.

This script:
1. Re-verify n=4 exhaustive (should still work — it's all LS)
2. Test swap chain with INCREASING number of swaps to see where it breaks
3. Show that the criterion was an artifact of proximity to cyclic LS
4. Determine what IS true universally
"""
import sys as _sys
try:
    _sys.stdout.reconfigure(encoding="utf-8")
    _sys.stderr.reconfigure(encoding="utf-8")
except Exception:
    pass

import numpy as np
from collections import Counter
from itertools import permutations
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

def get_A(L):
    n = L.shape[0]
    A = np.zeros((n - 1, n - 1), dtype=int)
    for i in range(n - 1):
        for j in range(n - 1):
            A[i, j] = L[i, j] - L[i, n - 1]
    return A

def v2(x):
    """2-adic valuation."""
    if x == 0:
        return 999
    x = abs(int(x))
    v = 0
    while x % 2 == 0:
        x //= 2
        v += 1
    return v

def cyclic_ls(n):
    L = np.zeros((n, n), dtype=int)
    for i in range(n):
        for j in range(n):
            L[i, j] = (i + j) % n + 1
    return L

def swap_chain_ls(n, target_successes):
    L = cyclic_ls(n)
    done = 0
    for _ in range(target_successes * 50):
        if done >= target_successes:
            break
        r1, r2 = random.sample(range(n), 2)
        c1, c2 = random.sample(range(n), 2)
        if L[r1, c1] == L[r2, c2] and L[r1, c2] == L[r2, c1]:
            L[r1, c1], L[r1, c2] = L[r1, c2], L[r1, c1]
            L[r2, c1], L[r2, c2] = L[r2, c2], L[r2, c1]
            done += 1
    return L

def isotope_ls(L):
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

# =================================================================
print("=" * 72)
print("PART 1: n=4 EXHAUSTIVE — The ground truth")
print("=" * 72)

squares_4 = []
perms = list(permutations(range(1, 5)))
for r0 in perms:
    if r0 != (1, 2, 3, 4):
        continue
    for r1 in perms:
        if len(set(zip(r0, r1))) < 4:
            continue
        cu = [set() for _ in range(4)]
        for j in range(4):
            cu[j].add(r0[j])
            cu[j].add(r1[j])
        for r2 in perms:
            ok = all(r2[j] not in cu[j] for j in range(4))
            if not ok:
                continue
            cu2 = [s.copy() for s in cu]
            for j in range(4):
                cu2[j].add(r2[j])
            for r3 in perms:
                ok = all(r3[j] not in cu2[j] for j in range(4))
                if not ok:
                    continue
                squares_4.append(np.array([r0, r1, r2, r3], dtype=int))

print(f"Total 4×4 LS (first row fixed): {len(squares_4)}")
match_4 = 0
mismatch_4 = 0
for L in squares_4:
    A = get_A(L)
    det_val = bareiss(A)
    B = A % 2
    B2 = (B @ B) % 2
    kdim2 = kernel_dim_F2(B2)
    criterion = (kdim2 >= 2)
    actual = (det_val % 16 == 0) if det_val != 0 else True
    if criterion == actual:
        match_4 += 1
    else:
        mismatch_4 += 1
        print(f"  MISMATCH: det={det_val}, v2={v2(det_val)}, kdim(B²)={kdim2}, "
              f"crit={'PASS' if criterion else 'FAIL'}, actual={'PASS' if actual else 'FAIL'}")

print(f"Result: {match_4}/{match_4+mismatch_4} match, {mismatch_4} mismatches")
if mismatch_4 == 0:
    print("★★★ n=4 EXHAUSTIVE: CRITERION HOLDS UNIVERSALLY ★★★")

# =================================================================
print()
print("=" * 72)
print("PART 2: n=8 swap chain — criterion vs number of swaps")
print("=" * 72)

for target_swaps in [5, 10, 20, 50, 100, 200, 500, 1000, 3000]:
    N = 200
    match = 0
    for _ in range(N):
        L = swap_chain_ls(8, target_swaps)
        A = get_A(L)
        det_val = bareiss(A)
        B = A % 2
        B2 = (B @ B) % 2
        kdim2 = kernel_dim_F2(B2)
        criterion = (kdim2 >= 2)
        actual = (det_val % 64 == 0) if det_val != 0 else True
        if criterion == actual:
            match += 1
    print(f"  swaps={target_swaps:>5d}: {match}/{N} match ({100*match/N:.0f}%)")

# =================================================================
print()
print("=" * 72)
print("PART 3: n=8 isotope — mismatch analysis")
print("=" * 72)

N = 300
type_counts = Counter()  # (criterion, actual) pairs

for _ in range(N):
    L = isotope_ls(cyclic_ls(8))
    A = get_A(L)
    det_val = bareiss(A)
    B = A % 2
    B2 = (B @ B) % 2
    kdim2 = kernel_dim_F2(B2)
    criterion = (kdim2 >= 2)
    actual = (det_val % 64 == 0) if det_val != 0 else True
    type_counts[(criterion, actual)] += 1

print(f"  (crit=PASS, act=PASS): {type_counts[(True, True)]}")
print(f"  (crit=PASS, act=FAIL): {type_counts[(True, False)]}")
print(f"  (crit=FAIL, act=PASS): {type_counts[(False, True)]}")
print(f"  (crit=FAIL, act=FAIL): {type_counts[(False, False)]}")

# =================================================================
print()
print("=" * 72)
print("PART 4: What IS universally true? Check 1 ∈ ker(B)")
print("=" * 72)

for n in [4, 8, 12]:
    N = 200
    ones_in_ker = 0
    for _ in range(N):
        L = isotope_ls(cyclic_ls(n))
        A = get_A(L)
        B = A % 2
        ones = np.ones(n - 1, dtype=int)
        if all((B @ ones) % 2 == 0):
            ones_in_ker += 1
    print(f"  n={n}: 1 ∈ ker(B) in {ones_in_ker}/{N} isotope LS ({100*ones_in_ker/N:.0f}%)")

# =================================================================
print()
print("=" * 72)
print("PART 5: What about dim ker(A) ≥ 1 → n|det? (weaker claim)")
print("=" * 72)

for n in [4, 8, 12]:
    N = 200
    match = 0
    for _ in range(N):
        L = isotope_ls(cyclic_ls(n))
        A = get_A(L)
        det_val = bareiss(A)
        kdim = kernel_dim_F2(A)
        criterion = (kdim >= 1)
        actual = (det_val % n == 0) if det_val != 0 else True
        if criterion == actual:
            match += 1
    print(f"  n={n}: dim ker≥1 ↔ n|det: {match}/{N} ({100*match/N:.0f}%)")

# =================================================================
print()
print("=" * 72)
print("PART 6: v₂(det) analysis for isotope LS")
print("=" * 72)

for n in [4, 8]:
    v2_vals = []
    kdims = []
    N = 300
    for _ in range(N):
        L = isotope_ls(cyclic_ls(n))
        A = get_A(L)
        det_val = bareiss(A)
        v2_vals.append(v2(det_val))
        kdims.append(kernel_dim_F2(A))
    
    print(f"n={n}:")
    print(f"  v₂(det) dist: {dict(sorted(Counter(v2_vals).items()))}")
    print(f"  kdim dist:     {dict(sorted(Counter(kdims).items()))}")
    print(f"  min v₂(det) = {min(v2_vals)}, need ≥ {2 * v2(n)} for n²|det")

    # Cross-tabulate
    cross = Counter()
    for v2v, k in zip(v2_vals, kdims):
        cross[(k, v2v >= 2 * v2(n))] += 1
    print(f"  Cross-tab (kdim, n²|det):")
    for (k, p), cnt in sorted(cross.items()):
        print(f"    kdim={k}, {'PASS' if p else 'FAIL'}: {cnt}")

# =================================================================
print()
print("=" * 72)
print("PART 7: Is the PAPER's Theorem 4.3 (n≡2 mod 4) still valid?")
print("=" * 72)

for n in [6, 10]:
    N = 200
    all_pass = True
    for _ in range(N):
        L = isotope_ls(cyclic_ls(n))
        A = get_A(L)
        det_val = bareiss(A)
        if det_val != 0 and det_val % (n * n) != 0:
            all_pass = False
            break
    if all_pass:
        print(f"  n={n}: n²|det for ALL {N} isotope LS ✓")
    else:
        print(f"  n={n}: FOUND FAILURE! n²∤det for some isotope LS ✗")
