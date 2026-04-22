"""
verify_final_clean.py — FINAL definitive test. NO isotopy with symbol permutation.

Key insight from previous tests:
- Symbol permutation changes A values and can break divisibility.
- The paper's theorem is about STANDARD symbols {1,...,n}.
- Swap chain preserves standard symbols.
- Previous "perfect" criterion was with too-few swaps (close to cyclic).

This test:
1. Swap chain with various swap counts — check dim ker(B) ≥ 2 → n²|det
2. Check BOTH directions of biconditional separately
3. Re-verify paper's Theorem 4.3 (n≡2 mod 4)
4. n=4 exhaustive with no-symbol-permutation isotopy
"""

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

def v2(x):
    if x == 0:
        return 999
    x = abs(int(x))
    v = 0
    while x % 2 == 0:
        x //= 2
        v += 1
    return v

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

def swap_chain_ls(n, target):
    L = cyclic_ls(n)
    done = 0
    for _ in range(target * 50):
        if done >= target:
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
print("PART 1: n=4 EXHAUSTIVE — dim ker(B) ≥ 2 vs n²|det")
print("=" * 72)

squares_4 = []
perms_list = list(permutations(range(1, 5)))
for r0 in perms_list:
    if r0 != (1, 2, 3, 4):
        continue
    for r1 in perms_list:
        if len(set(zip(r0, r1))) < 4:
            continue
        cu = [set() for _ in range(4)]
        for j in range(4):
            cu[j].add(r0[j])
            cu[j].add(r1[j])
        for r2 in perms_list:
            ok = all(r2[j] not in cu[j] for j in range(4))
            if not ok:
                continue
            cu2 = [s.copy() for s in cu]
            for j in range(4):
                cu2[j].add(r2[j])
            for r3 in perms_list:
                ok = all(r3[j] not in cu2[j] for j in range(4))
                if not ok:
                    continue
                squares_4.append(np.array([r0, r1, r2, r3], dtype=int))

print(f"All 4×4 LS (row 0 fixed): {len(squares_4)}")

# Cross-tabulate
cross_4 = Counter()
for L in squares_4:
    A = get_A(L)
    det_val = bareiss(A)
    kdim = kernel_dim_F2(A)
    passes = (det_val % 16 == 0) if det_val != 0 else True
    cross_4[(kdim, passes)] += 1

print("Cross-tab (kdim(B), 16|det):")
for (k, p), cnt in sorted(cross_4.items()):
    print(f"  kdim={k}, {'PASS' if p else 'FAIL'}: {cnt}")

# Also check without det=0
print("\nExcluding det=0:")
cross_4b = Counter()
for L in squares_4:
    A = get_A(L)
    det_val = bareiss(A)
    if det_val == 0:
        continue
    kdim = kernel_dim_F2(A)
    passes = (det_val % 16 == 0)
    cross_4b[(kdim, passes)] += 1

for (k, p), cnt in sorted(cross_4b.items()):
    print(f"  kdim={k}, {'PASS' if p else 'FAIL'}: {cnt}")

# =====================================================================
print()
print("=" * 72)
print("PART 2: n=8, swap chain with increasing swaps")
print("       Checking dim ker(B) ≥ 2 and dim ker(B²) ≥ 2 separately")
print("=" * 72)

for swaps in [5, 10, 20, 50, 100, 500, 2000, 5000]:
    N = 300
    # Counts: (kdim_criterion, actual)
    kB_cross = Counter()
    kB2_cross = Counter()
    
    for _ in range(N):
        L = swap_chain_ls(8, swaps)
        A = get_A(L)
        det_val = bareiss(A)
        B = A % 2
        B2 = (B @ B) % 2
        kdim = kernel_dim_F2(A)
        kdim2 = kernel_dim_F2(B2)
        actual = (det_val % 64 == 0) if det_val != 0 else True
        
        kB_cross[(kdim >= 2, actual)] += 1
        kB2_cross[(kdim2 >= 2, actual)] += 1
    
    # dim ker(B) biconditional
    kB_match = kB_cross[(True, True)] + kB_cross[(False, False)]
    kB_suff = kB_cross[(True, False)]  # crit PASS, actual FAIL (sufficiency fails)
    kB_nec = kB_cross[(False, True)]   # crit FAIL, actual PASS (necessity fails)
    
    # dim ker(B²) biconditional
    kB2_match = kB2_cross[(True, True)] + kB2_cross[(False, False)]
    kB2_suff = kB2_cross[(True, False)]
    kB2_nec = kB2_cross[(False, True)]
    
    print(f"  swaps={swaps:>5d}: ker(B)≥2 match={kB_match}/{N} (suff_fail={kB_suff}, nec_fail={kB_nec})"
          f" | ker(B²)≥2 match={kB2_match}/{N} (suff_fail={kB2_suff}, nec_fail={kB2_nec})")

# =====================================================================
print()
print("=" * 72)
print("PART 3: n=8, swap chain 5000 — full cross-tabulation")
print("=" * 72)

N = 500
cross_8 = Counter()

for _ in range(N):
    L = swap_chain_ls(8, 5000)
    A = get_A(L)
    det_val = bareiss(A)
    kdim = kernel_dim_F2(A)
    v2_det = v2(det_val)
    actual = (det_val % 64 == 0) if det_val != 0 else True
    cross_8[(kdim, actual)] += 1

print("Cross-tab (kdim(B), 64|det) for n=8 swap chain (5000 swaps):")
for (k, p), cnt in sorted(cross_8.items()):
    print(f"  kdim={k}, {'PASS' if p else 'FAIL'}: {cnt}")

# =====================================================================
print()
print("=" * 72)
print("PART 4: n=12, swap chain 5000")
print("=" * 72)

N = 300
cross_12 = Counter()

for _ in range(N):
    L = swap_chain_ls(12, 5000)
    A = get_A(L)
    det_val = bareiss(A)
    kdim = kernel_dim_F2(A)
    actual = (det_val % 144 == 0) if det_val != 0 else True
    cross_12[(kdim, actual)] += 1

print("Cross-tab (kdim(B), 144|det) for n=12 swap chain (5000 swaps):")
for (k, p), cnt in sorted(cross_12.items()):
    print(f"  kdim={k}, {'PASS' if p else 'FAIL'}: {cnt}")

# =====================================================================
print()
print("=" * 72)
print("PART 5: Paper's Theorem 4.3 (n≡2 mod 4) with swap chain only")
print("=" * 72)

for n in [6, 10]:
    N = 300
    all_pass = True
    fail_example = None
    for trial in range(N):
        L = swap_chain_ls(n, 3000)
        A = get_A(L)
        det_val = bareiss(A)
        if det_val != 0 and det_val % (n * n) != 0:
            all_pass = False
            fail_example = (det_val, v2(det_val))
            break
    if all_pass:
        print(f"  n={n}: n²|det for ALL {N} swap-chain LS ✓")
    else:
        print(f"  n={n}: FAILURE at trial {trial}: det={fail_example[0]}, v₂={fail_example[1]}")

# =====================================================================
print()
print("=" * 72)
print("PART 6: n≡2 mod 4 — exhaustive at n=6")
print("=" * 72)

# Generate all 6×6 LS with first row = (1,2,3,4,5,6)?
# There are 9408 such LS. Too many to enumerate by brute force with permutations.
# Let's test with heavy swap chain instead.
n = 6
N = 1000
n2 = n * n
fails = 0
v2_dets = []
for _ in range(N):
    L = swap_chain_ls(n, 5000)
    A = get_A(L)
    det_val = bareiss(A)
    v2_dets.append(v2(det_val))
    if det_val != 0 and det_val % n2 != 0:
        fails += 1

print(f"n={n}: {fails}/{N} fail n²|det with swap chain (5000 swaps)")
print(f"  v₂(det) dist: {dict(sorted(Counter(v2_dets).items()))}")
print(f"  Need v₂ ≥ {2*v2(n)} for n²|det")

# =====================================================================
print()
print("=" * 72)
print("PART 7: Verify that base column choice preserves |det|")
print("=" * 72)

n = 8
for trial in range(5):
    L = swap_chain_ls(n, 3000)
    dets = []
    for base in range(n):
        cols = [j for j in range(n) if j != base]
        A = np.zeros((n - 1, n - 1), dtype=int)
        for i in range(n - 1):
            for idx, j in enumerate(cols):
                A[i, idx] = L[i, j] - L[i, base]
        dets.append(abs(int(bareiss(A))))
    unique = set(dets)
    print(f"  Trial {trial}: |det| values = {unique} ({'ALL EQUAL' if len(unique)==1 else 'DIFFER!'})")
