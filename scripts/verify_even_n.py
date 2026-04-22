"""
Verify: |det(A_cyc)| = n^{n-2} for even n.
Then: prove why Prop 3.6's DFT proof works for even n.

For odd n: prod_{k=1}^{n-1} (omega^k - 1) = n  (cyclotomic identity)
For even n: prod_{k=1}^{n-1} (omega^k - 1) = ???

The cyclotomic identity: x^n - 1 = prod_{k=0}^{n-1} (x - omega^k)
Differentiating or evaluating at x=1:
n = prod_{k=1}^{n-1} (1 - omega^k)

So prod_{k=1}^{n-1} (omega^k - 1) = prod_{k=1}^{n-1} (-1)(1-omega^k) 
                                    = (-1)^{n-1} * prod_{k=1}^{n-1} (1 - omega^k)
                                    = (-1)^{n-1} * n

For odd n: (-1)^{n-1} = 1, so product = n. ✓ (matches Prop 3.6)
For even n: (-1)^{n-1} = -1, so product = -n.

det(A) = n^{n-1} / prod_{k=1}^{n-1} (omega^k - 1) = n^{n-1} / ((-1)^{n-1} * n) = (-1)^{n-1} * n^{n-2}

For odd n: det(A) has sign from conjugate pairs: (-1)^{(n-1)/2}
For even n: det(A) has sign... let me compute.

Actually, the EIGENVALUES are lambda_k = n / (omega^k - 1).
det(A) = prod_{k=1}^{n-1} lambda_k = n^{n-1} / prod_{k=1}^{n-1}(omega^k - 1)
        = n^{n-1} / ((-1)^{n-1} · n) = (-1)^{n-1} · n^{n-2}

For even n: (-1)^{n-1} = (-1)^{odd} = -1.
So det(A_cyc) = -n^{n-2} for even n ≥ 2.
|det(A_cyc)| = n^{n-2} for ALL n ≥ 2. 

Wait — but det(E_std) = n · det(A). For even n:
det(E_std) = n · (-n^{n-2}) = -n^{n-1}.

Let me verify numerically.
"""
import sys as _sys
try:
    _sys.stdout.reconfigure(encoding="utf-8")
    _sys.stderr.reconfigure(encoding="utf-8")
except Exception:
    pass
import numpy as np

def cyclic_ls(n):
    return np.array([[(i+j) % n + 1 for j in range(n)] for i in range(n)])

def get_A(L):
    n = L.shape[0]
    return np.array([[int(L[i,j]) - int(L[i,n-1]) for j in range(n-1)] for i in range(n-1)])

def bareiss(M):
    m = len(M); A = [row[:] for row in M]; sign = 1; prev = 1
    for i in range(m):
        if A[i][i] == 0:
            for k in range(i+1, m):
                if A[k][i] != 0: A[i], A[k] = A[k], A[i]; sign *= -1; break
            else: return 0
        for j in range(i+1, m):
            for k in range(i+1, m):
                A[j][k] = (A[j][k]*A[i][i] - A[j][i]*A[i][k]) // prev
            A[j][i] = 0
        prev = A[i][i]
    return sign * A[m-1][m-1]

print("Verify |det(A_cyc)| = n^{n-2} for ALL n >= 2")
print("=" * 60)

for n in range(2, 16):
    L = cyclic_ls(n)
    A = get_A(L)
    det_val = bareiss(A.tolist())
    expected_abs = n ** (n-2)
    
    # Predicted sign: (-1)^{n-1}
    predicted_sign = (-1) ** (n-1)
    predicted_det = predicted_sign * expected_abs
    
    match_abs = (abs(det_val) == expected_abs)
    match_sign = (det_val == predicted_det)
    
    print(f"  n={n:2d}: det(A)={det_val:>15d}  "
          f"|det|={abs(det_val):>12d}  n^(n-2)={expected_abs:>12d}  "
          f"|det|=n^(n-2):{match_abs}  "
          f"sign=(-1)^(n-1):{match_sign}")

print()
print("CONCLUSION:")
print("  The cyclotomic identity prod_{k=1}^{n-1}(omega^k-1) = (-1)^{n-1} · n")
print("  holds for ALL n >= 2, not just odd n.")
print("  Therefore: det(A_cyc) = (-1)^{n-1} · n^{n-2} for all n >= 2.")
print("  |det(A_cyc)| = n^{n-2} for all n >= 2.")
print()
print("  The proof in Prop 3.6 uses the formula for ODD n only because")
print("  the sign analysis referenced (n-1)/2 conjugate pairs.")
print("  For even n, the sign is simply (-1)^{n-1} = -1.")
print()
print("  The circulant structure A_{ij} = ((i+j) mod n + 1) - ((i+n-1) mod n + 1)")
print("  = (i+j) mod n - (i+n-1) mod n works identically for even n.")
print("  The DFT diagonalization and product formula are the same.")
print()

# Actually let's check: is E circulant for even n?  
# E_ij = L_ij - (n+1)/2. For even n, (n+1)/2 is not an integer!
# So the "centered matrix" E uses a different centering for even n.
# Let me re-read what centring the paper uses.

# Actually the paper defines A_ij = L_{ij} - L_{i,n-1}, which is always integer.
# The circulant structure of A for the cyclic LS:
# A_ij = ((i+j) mod n + 1) - ((i + n-1) mod n + 1) = (i+j) mod n - (i+n-1) mod n
# For n=4, i=0: A_{0,j} = j mod 4 - 3 mod 4 = j - 3 for j=0,1,2
#   = [-3, -2, -1] ✓ (first row of circulant)
# For n=4, i=1: A_{1,j} = (1+j) mod 4 - 0 = 1+j for j=0,1,2
#   = [1, 2, 3] 
# For n=4, i=2: A_{2,j} = (2+j) mod 4 - 1  
#   j=0: 2-1=1, j=1: 3-1=2, j=2: 0-1=-1... wait, that's [1, 2, -1]
# Hmm, that doesn't look circulant.

print("Is A circulant for even n?")
for n in [4, 5, 6]:
    L = cyclic_ls(n)
    A = get_A(L)
    m = n - 1
    first_row = tuple(A[0])
    is_circulant = True
    for i in range(1, m):
        expected = tuple(A[0, (j - i) % m] for j in range(m))
        actual = tuple(A[i])
        if actual != expected:
            is_circulant = False
    print(f"  n={n}: circulant={is_circulant}")
    if not is_circulant:
        print(f"    Row 0: {list(A[0])}")
        print(f"    Row 1: {list(A[1])}")
        print(f"    Row 2: {list(A[2])}")

# Even if A is not circulant, |det(A)| = n^{n-2} still holds.
# The proof path for even n might be different from the circulant eigenvalue argument.
# Let me check if A is circulant mod something.

print()
print("Is A_cyc circulant (n-1 × n-1) for general n?")
print("A_ij = (i+j) mod n - (i+n-1) mod n, for i,j = 0,...,n-2")
print()

for n in [3, 4, 5, 6, 7, 8]:
    L = cyclic_ls(n)
    A = get_A(L)
    m = n - 1
    # Check if it's a circulant of size (n-1)x(n-1)
    first_row = A[0]
    circ = True
    for i in range(1, m):
        for j in range(m):
            if A[i][j] != first_row[(j - i) % m]:
                circ = False; break
        if not circ: break
    print(f"  n={n}: A is {m}×{m} circulant: {circ}")
