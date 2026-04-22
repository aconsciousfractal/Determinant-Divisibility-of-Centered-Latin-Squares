"""Definitive verification: n^2 | det(E_std) for n=8 using Jacobson-Matthews.
This settles whether Corollary 3.3/4.4 are correct or the dossier's v_2=5 claim."""
import numpy as np


def latin_to_incidence(L):
    n = L.shape[0]
    cube = np.zeros((n, n, n), dtype=int)
    for i in range(n):
        for j in range(n):
            cube[i, j, L[i, j] - 1] = 1
    return cube


def incidence_to_latin(cube):
    n = cube.shape[0]
    L = np.zeros((n, n), dtype=int)
    for i in range(n):
        for j in range(n):
            for k in range(n):
                if cube[i, j, k] == 1:
                    L[i, j] = k + 1
                    break
    return L


def is_proper(cube):
    return not np.any(cube < 0)


def jm_step(cube, rng):
    n = cube.shape[0]
    neg = np.argwhere(cube == -1)
    if len(neg) == 0:
        zeros = np.argwhere(cube == 0)
        idx = rng.integers(len(zeros))
        i0, j0, k0 = zeros[idx]
    else:
        i0, j0, k0 = neg[0]

    k1 = np.where(cube[i0, j0, :] == 1)[0][0]
    j1 = np.where(cube[i0, :, k0] == 1)[0][0]
    i1 = np.where(cube[:, j0, k0] == 1)[0][0]

    new_cube = cube.copy()
    new_cube[i0, j0, k0] += 1
    new_cube[i0, j1, k1] += 1
    new_cube[i1, j0, k1] += 1
    new_cube[i1, j1, k0] += 1
    new_cube[i0, j0, k1] -= 1
    new_cube[i0, j1, k0] -= 1
    new_cube[i1, j0, k0] -= 1
    new_cube[i1, j1, k1] -= 1
    return new_cube, is_proper(new_cube)


def gen_random_ls(n, rng, burn_in=300):
    L = np.array([[(i + j) % n + 1 for j in range(n)] for i in range(n)])
    cube = latin_to_incidence(L)
    proper_count = 0
    steps = 0
    while proper_count < burn_in and steps < burn_in * 50:
        cube, prop = jm_step(cube, rng)
        steps += 1
        if prop:
            proper_count += 1
    if not is_proper(cube):
        for _ in range(5000):
            cube, prop = jm_step(cube, rng)
            if prop:
                break
    return incidence_to_latin(cube)


def bareiss_det(M):
    n = len(M)
    A = [list(row) for row in M]
    sign = 1
    for k in range(n):
        if A[k][k] == 0:
            found = False
            for i in range(k+1, n):
                if A[i][k] != 0:
                    A[k], A[i] = A[i], A[k]
                    sign *= -1
                    found = True
                    break
            if not found:
                return 0
        for i in range(k+1, n):
            for j in range(k+1, n):
                A[i][j] = A[k][k]*A[i][j] - A[i][k]*A[k][j]
                if k > 0:
                    A[i][j] //= A[k-1][k-1]
            A[i][k] = 0
    return sign * A[n-1][n-1]


def v2(x):
    if x == 0: return float('inf')
    x = abs(x); v = 0
    while x % 2 == 0: v += 1; x //= 2
    return v


def is_valid_ls(L):
    n = L.shape[0]
    for i in range(n):
        if sorted(L[i]) != list(range(1, n+1)):
            return False
        if sorted(L[:, i]) != list(range(1, n+1)):
            return False
    return True


rng = np.random.default_rng(42)
n = 8
print(f"=== Definitive test: n = {n}, n^2 = {n**2}, v_2(n^2) = {v2(n**2)} ===")
print(f"Generating random LS via Jacobson-Matthews (burn-in=300)")

min_v2_E = float('inf')
min_v2_A = float('inf')
fails_n2 = 0
count = 0
v2_hist = {}

for trial in range(500):
    L = gen_random_ls(n, rng)
    if not is_valid_ls(L):
        continue
    
    A = [[int(L[i][j] - L[i][n-1]) for j in range(n-1)] for i in range(n-1)]
    dA = bareiss_det(A)
    dE = n * dA
    if dE == 0:
        continue
    count += 1
    
    ve = v2(dE)
    va = v2(dA)
    v2_hist[ve] = v2_hist.get(ve, 0) + 1
    
    if ve < min_v2_E:
        min_v2_E = ve
        print(f"  Trial {trial}: NEW MIN v_2(det E) = {ve}, det(A)={dA}, v_2(A)={va}")
    
    if va < min_v2_A:
        min_v2_A = va
    
    if dE % (n**2) != 0:
        fails_n2 += 1
        if fails_n2 <= 5:
            print(f"  *** n^2 FAILS at trial {trial}: det(A)={dA}, v_2(A)={va}, v_2(E)={ve}")

print(f"\n=== Results: n = {n} ===")
print(f"Tested: {count} valid LS")
print(f"Min v_2(det E_std) = {min_v2_E}, required for n^2: {v2(n**2)}")
print(f"Min v_2(det A) = {min_v2_A}, required for n: {v2(n)}")
print(f"n^2-divisibility failures: {fails_n2} / {count} ({100*fails_n2/max(count,1):.1f}%)")
print(f"\nv_2 histogram: {dict(sorted(v2_hist.items()))}")

if fails_n2 > 0:
    print(f"\n*** THEOREM IS FALSE for n={n}: n^2 does NOT always divide det(E_std)")
    print(f"*** Paper's Corollary 3.3 and 4.4 need correction")
else:
    print(f"\nCONFIRMED: n^2={n**2} | det(E_std) for all {count} tested.")
