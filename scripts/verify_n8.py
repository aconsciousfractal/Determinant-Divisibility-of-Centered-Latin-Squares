"""Verify n^2 | det(E_std) for n=8 Latin squares.

Check whether Corollary 3.3 / Corollary 4.4 are correct:
for n ≡ 0 mod 4, does n^2 always divide det(E_std)?

Also check the Open Question 5 claim: v_2(det E_std) = 5 minimum at n=8.
"""
import sys as _sys
try:
    _sys.stdout.reconfigure(encoding="utf-8")
    _sys.stderr.reconfigure(encoding="utf-8")
except Exception:
    pass

import random
import math

def bareiss_det(M):
    """Exact integer determinant via Bareiss algorithm."""
    n = len(M)
    A = [row[:] for row in M]
    sign = 1
    for k in range(n):
        # Pivot
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
                A[i][j] = A[k][k] * A[i][j] - A[i][k] * A[k][j]
                if k > 0:
                    A[i][j] //= A[k-1][k-1]
            A[i][k] = 0
    return sign * A[n-1][n-1]


def v2(x):
    """2-adic valuation of integer x."""
    if x == 0:
        return float('inf')
    x = abs(x)
    v = 0
    while x % 2 == 0:
        v += 1
        x //= 2
    return v


def generate_random_ls(n, num_swaps=10000):
    """Generate random Latin square of order n via Jacobson-Matthews-like swaps."""
    # Start with cyclic Latin square
    L = [[(i + j) % n + 1 for j in range(n)] for i in range(n)]
    
    for _ in range(num_swaps):
        # Random row swap
        i1, i2 = random.sample(range(n), 2)
        L[i1], L[i2] = L[i2], L[i1]
        
        # Random column swap
        j1, j2 = random.sample(range(n), 2)
        for i in range(n):
            L[i][j1], L[i][j2] = L[i][j2], L[i][j1]
        
        # Random symbol swap  
        s1, s2 = random.sample(range(1, n+1), 2)
        for i in range(n):
            for j in range(n):
                if L[i][j] == s1:
                    L[i][j] = s2
                elif L[i][j] == s2:
                    L[i][j] = s1
    
    return L


def is_latin_square(L):
    n = len(L)
    for i in range(n):
        if sorted(L[i]) != list(range(1, n+1)):
            return False
        col = sorted(L[i][j] for j in range(n) for _ in [None] if False)
        col = [L[j][i] for j in range(n)]
        if sorted(col) != list(range(1, n+1)):
            return False
    return True


def compute_A(L):
    """Compute the (n-1)x(n-1) difference matrix A."""
    n = len(L)
    A = [[L[i][j] - L[i][n-1] for j in range(n-1)] for i in range(n-1)]
    return A


def compute_det_Estd(L):
    n = len(L)
    A = compute_A(L)
    det_A = bareiss_det(A)
    det_Estd = n * det_A
    return det_Estd, det_A


def main():
    random.seed(42)
    n = 8
    print(f"=== Verifying n^2 | det(E_std) for n = {n} ===")
    print(f"n^2 = {n**2}, v_2(n^2) = {v2(n**2)}")
    print()
    
    min_v2_Estd = float('inf')
    min_v2_detA = float('inf')
    num_tested = 0
    num_fail_n2 = 0
    
    for trial in range(500):
        L = generate_random_ls(n)
        if not is_latin_square(L):
            continue
        
        det_Estd, det_A = compute_det_Estd(L)
        
        if det_Estd == 0:
            continue
            
        num_tested += 1
        v2_Estd = v2(det_Estd)
        v2_A = v2(det_A)
        
        if v2_Estd < min_v2_Estd:
            min_v2_Estd = v2_Estd
            print(f"  New min v_2(det E_std) = {v2_Estd} at trial {trial}")
            print(f"    det(A) = {det_A}, v_2(det A) = {v2_A}")
            print(f"    det(E_std) = {det_Estd}")
            print(f"    n^2 | det(E_std)? {det_Estd % (n**2) == 0}")
        
        if v2_A < min_v2_detA:
            min_v2_detA = v2_A
        
        if det_Estd % (n**2) != 0:
            num_fail_n2 += 1
            print(f"  *** COUNTEREXAMPLE at trial {trial}: n^2 DOES NOT divide det(E_std)")
            print(f"    det(E_std) = {det_Estd}, v_2 = {v2_Estd}")
    
    print()
    print(f"=== Results for n = {n} ===")
    print(f"Tested: {num_tested} Latin squares")
    print(f"Min v_2(det E_std): {min_v2_Estd}")
    print(f"Min v_2(det A): {min_v2_detA}")
    print(f"n^2-divisibility failures: {num_fail_n2}")
    
    if num_fail_n2 == 0:
        print(f"\nCONFIRMED: n^2 = {n**2} divides det(E_std) for all {num_tested} tested squares.")
        print(f"The minimum v_2(det E_std) = {min_v2_Estd}, v_2(n^2) = {v2(n**2)}")
    else:
        print(f"\n*** THEOREM IS FALSE: {num_fail_n2} counterexamples found!")
    
    # Also test n=4 and n=12
    for n2 in [4, 12]:
        print(f"\n=== Quick check n = {n2} ===")
        fails = 0
        minv = float('inf')
        tested = 0
        for _ in range(200):
            L = generate_random_ls(n2)
            if not is_latin_square(L):
                continue
            det_Estd, det_A = compute_det_Estd(L)
            if det_Estd == 0:
                continue
            tested += 1
            v = v2(det_Estd)
            if v < minv:
                minv = v
            if det_Estd % (n2**2) != 0:
                fails += 1
        print(f"Tested: {tested}, Min v_2: {minv}, v_2(n^2) = {v2(n2**2)}, Failures: {fails}")
        if fails == 0:
            print(f"CONFIRMED: n^2 = {n2**2} | det(E_std) for all tested.")
        else:
            print(f"*** FAILURES at n={n2}: {fails}")


if __name__ == "__main__":
    main()
