"""Run T3 from test_T2_T6.py to check if v_2 = 5 is achievable at n=8."""
import numpy as np
from collections import Counter

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

n = 8
L_cyc = np.array([[(i+j) % n + 1 for j in range(n)] for i in range(n)], dtype=np.int64)

# Non-cyclic base: Cayley table of Z4 x Z2
L_nc8 = np.array([
    [1,2,3,4,5,6,7,8],[2,1,4,3,6,5,8,7],[3,4,1,2,7,8,5,6],
    [4,3,2,1,8,7,6,5],[5,6,7,8,1,2,3,4],[6,5,8,7,2,1,4,3],
    [7,8,5,6,3,4,1,2],[8,7,6,5,4,3,2,1]], dtype=np.int64)

bases = [L_cyc]
if all(len(set(L_nc8[i,:])) == 8 for i in range(8)) and all(len(set(L_nc8[:,j])) == 8 for j in range(8)):
    bases.append(L_nc8)
    print(f"Using {len(bases)} base LS")

v2_dist = Counter()
min_v2 = 999
min_ex = None
total = 0
n2_fails = 0

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
        v = v_p(dE, 2)
        v2_dist[v] += 1
        if v < min_v2:
            min_v2 = v
            min_ex = (bi, seed, dA, dE, L.copy())
        if dE % (n**2) != 0:
            n2_fails += 1

print(f"\nTotal non-zero: {total}")
print(f"\nv_2(det E_std) distribution:")
for k in sorted(v2_dist):
    pct = 100 * v2_dist[k] / total
    print(f"  v_2 = {k:3d}: {v2_dist[k]:5d}  ({pct:.1f}%)")
print(f"\nMin v_2(det E_std) = {min_v2}")
print(f"v_2(n^2) = v_2({n**2}) = {v_p(n**2, 2)}")
print(f"v_2(n^2/2) = v_2({n**2//2}) = {v_p(n**2//2, 2)}")
print(f"\nn^2 divisibility failures (n^2 ∤ det E): {n2_fails}/{total} ({100*n2_fails/total:.1f}%)")

if min_v2 == 5:
    print("\n*** CONFIRMED: v_2 = 5 IS achievable. n^2 does NOT always divide det E.")
    print("*** Corollary 3.3 (n ≡ 0 mod 4 case) is FALSE.")
    print("*** Corollary 4.4 is FALSE.")
    print("*** Theorem 4.3 (for a ≥ 2) is FALSE.")
elif min_v2 >= 6:
    print(f"\nAll {total} squares have v_2 ≥ 6 = v_2(n^2). n^2 ALWAYS divides det E.")
    print("Corollary 3.3 and 4.4 appear correct. Dossier's v_2=5 claim may be wrong.")

if min_ex:
    bi, s, dA, dE, Lex = min_ex
    print(f"\nMinimum example: base={bi}, seed={s}")
    print(f"  det(A) = {dA}, v_2(A) = {v_p(dA, 2)}")
    print(f"  det(E) = {dE}, v_2(E) = {v_p(dE, 2)}")
    print(f"  n^2 | det(E)? {dE % (n**2) == 0}")
    print(f"  n^2/2 | det(E)? {dE % (n**2//2) == 0}")
