"""
Verify the unified sign formula:
  det(A_cyc) = (-1)^{floor(n/2)} * n^{n-2}   for all n >= 2.

Proof sketch for the sign:
1. E_{ij} = (i+j) mod n - (n-1)/2 is anti-circulant, not circulant.
2. Row permutation sigma: i -> (-i) mod n converts E to a standard circulant.
   sigma has sign (-1)^{floor((n-1)/2)}.
3. The circulant has eigenvalues on V_std: lambda_k = n/(omega^k-1).
4. Cyclotomic identity: prod_{k=1}^{n-1} (omega^k-1) = (-1)^{n-1} * n.
5. Product of eigenvalues = n^{n-1} / ((-1)^{n-1} * n) = (-1)^{n-1} * n^{n-2}.
   (unsigned product of the circulant, including cyclotomic sign)
6. Total sign = permutation_sign * cyclotomic_sign:
   = (-1)^{floor((n-1)/2)} * (-1)^{n-1}
   = (-1)^{floor((n-1)/2) + n - 1}
   
   For odd n: floor((n-1)/2) + n-1 = (n-1)/2 + n-1 = (3n-3)/2 = (3(n-1))/2
     Since n is odd, n-1 is even, so (3(n-1))/2 is even * 3/2... hmm.
     Actually: (n-1)/2 + n - 1. For n=3: 1+2=3, (-1)^3 = -1. floor(3/2)=1, (-1)^1=-1. ✓
     For n=5: 2+4=6, (-1)^6=1. floor(5/2)=2, (-1)^2=1. ✓
     For n=7: 3+6=9, (-1)^9=-1. floor(7/2)=3, (-1)^3=-1. ✓ 
   
   For even n: floor((n-1)/2) + n-1 = (n-2)/2 + n-1 = (3n-4)/2
     For n=2: 0+1=1, (-1)^1=-1. floor(2/2)=1, (-1)^1=-1. ✓
     For n=4: 1+3=4, (-1)^4=1. floor(4/2)=2, (-1)^2=1. ✓
     For n=6: 2+5=7, (-1)^7=-1. floor(6/2)=3, (-1)^3=-1. ✓

   General: floor((n-1)/2) + n - 1 ≡ floor(n/2) (mod 2)?
   Let's verify:
   - n even: (n-2)/2 + n-1 = (n-2+2n-2)/2 = (3n-4)/2 = 3n/2 - 2.
     floor(n/2) = n/2. Parity of 3n/2-2 vs n/2: differ by n-2, which is even. Same parity. ✓
   - n odd: (n-1)/2 + n-1 = (n-1+2n-2)/2 = (3n-3)/2.
     floor(n/2) = (n-1)/2. Parity of (3n-3)/2 vs (n-1)/2: differ by n-1, even for odd n. Same parity. ✓

So: (-1)^{floor((n-1)/2) + n-1} = (-1)^{floor(n/2)} for all n >= 2.

UNIFIED FORMULA: det(A_cyc) = (-1)^{floor(n/2)} * n^{n-2}.
"""

import numpy as np

def cyclic_ls(n):
    return [[(i+j) % n + 1 for j in range(n)] for i in range(n)]

def bareiss(M):
    m = len(M); A = [row[:] for row in M]; sign = 1; prev = 1
    for i in range(m):
        if A[i][i] == 0:
            for k in range(i+1, m):
                if A[k][i] != 0:
                    A[i], A[k] = A[k], A[i]; sign *= -1; break
            else: return 0
        for j in range(i+1, m):
            for k in range(i+1, m):
                A[j][k] = (A[j][k]*A[i][i] - A[j][i]*A[i][k]) // prev
            A[j][i] = 0
        prev = A[i][i]
    return sign * A[m-1][m-1]

print("Unified sign formula verification: det(A) = (-1)^{floor(n/2)} * n^{n-2}")
print("=" * 70)

all_ok = True
for n in range(2, 20):
    L = cyclic_ls(n)
    A = [[L[i][j] - L[i][n-1] for j in range(n-1)] for i in range(n-1)]
    det_A = bareiss(A)
    predicted = ((-1) ** (n // 2)) * (n ** (n - 2))
    ok = det_A == predicted
    if not ok:
        all_ok = False
    print(f"  n={n:2d}: det(A)={det_A:>20d}  predicted={predicted:>20d}  {'✓' if ok else '✗ FAIL'}")

print()
if all_ok:
    print("ALL PASS: det(A_cyc) = (-1)^{floor(n/2)} * n^{n-2} for n=2..19")
else:
    print("SOME FAILURES")

# Also verify: det(E_std) = n * det(A), so det(E_std) = (-1)^{floor(n/2)} * n^{n-1}
print()
print("Verify det(E_std) = n * det(A):")
for n in range(2, 12):
    L = cyclic_ls(n)
    E = np.array(L, dtype=float) - (n+1)/2
    # E_std means "the reduced (n-1)x(n-1) matrix of the centered matrix"
    # Or does it? Let me check: det(E full n×n) = 0 since rows sum to 0.
    # The cofactor of E (delete last row and col):
    E_reduced = E[:n-1, :n-1]
    det_E_red = round(np.linalg.det(E_reduced))
    
    A = [[L[i][j] - L[i][n-1] for j in range(n-1)] for i in range(n-1)]
    det_A = bareiss(A)
    
    # In the paper, E_std is the FULL n×n centered matrix, and
    # det(E_std) is computed as n * det(A) (well-known relation).
    # Actually, det(E_std) = 0 (singular). The paper likely means
    # something else by "det(E_std)" — maybe the product of nonzero eigenvalues.
    
    # Actually re-reading the paper: "det(E_std)_cyc" in the proposition
    # likely refers to det(A) since A is the "standard" reduced matrix.
    # Let me check what E_std means in the paper context.
    
    print(f"  n={n:2d}: det(A)={det_A:>14d}  det(E[:n-1,:n-1])={det_E_red:>14d}  "
          f"ratio={det_E_red/det_A if det_A != 0 else 'N/A':>8}")
