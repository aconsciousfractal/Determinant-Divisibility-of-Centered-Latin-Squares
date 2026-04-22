"""
Prove det(A_cyc) for ALL n via the n×n circulant E_std.

The key insight: E_std (n×n) IS circulant for ALL n. 
E_std = L - (n+1)/2 · J/n · J... no wait.

Actually L_{ij} = (i+j) mod n + 1  is circulant (n×n).
E_std_{ij} = L_{ij} - c_i where c_i is the row average = (n+1)/2.
So E_std = L - (n+1)/2 · J where J is the all-ones matrix? No.
E_std_{ij} = L_{ij} - mean of row i.
Every row of L is a permutation of {1,...,n}, so mean = (n+1)/2.
E_std_{ij} = L_{ij} - (n+1)/2  for all i,j.

Since L is circulant, E_std is also circulant (subtracting a constant from 
all entries preserves circulant structure).

E_std is n×n with row sums = 0 (since each row is {1,...,n} - (n+1)/2).
So det(E_std) = 0 (rank ≤ n-1). The reduced determinant det(A) relates to
the product of nonzero eigenvalues.

The eigenvalues of E_std: since E_std is n×n circulant with first row
c_j = (j mod n) + 1 - (n+1)/2 = j + 1 - (n+1)/2 = j - (n-1)/2
(for j = 0,...,n-1).

Wait, E_std first row: j=0: 0+1-(n+1)/2 = (1-n)/2
                         j=1: 1+1-(n+1)/2 = (3-n)/2
                         ...
                         j=n-1: n-1+1-(n+1)/2 = (n-1)/2

For even n, these are half-integers. But E_std is the n×n matrix, 
and its eigenvalues are:
lambda_k = sum_{j=0}^{n-1} c_j * omega^{jk}  for k=0,...,n-1
where omega = e^{2pi i/n}.

c_j = j - (n-1)/2.

lambda_k = sum_{j=0}^{n-1} j omega^{jk} - (n-1)/2 * sum_{j=0}^{n-1} omega^{jk}

For k ≥ 1: sum omega^{jk} = 0, so lambda_k = sum j omega^{jk}.

sum_{j=0}^{n-1} j omega^{jk} = omega^k * d/d(omega^k) [sum omega^{jk}]
= ... this is a standard DFT identity.

Actually: sum_{j=0}^{n-1} j z^j = z * d/dz * (z^n - 1)/(z-1) = z * (nz^{n-1}(z-1) - (z^n-1)) / (z-1)^2
For z = omega^k (n-th root of unity, k ≠ 0):
z^n = 1, so nz^{n-1}(z-1) - (1-1) = nz^{n-1}(z-1)
= n * omega^{k(n-1)} * (omega^k - 1)
= n * omega^{-k} * (omega^k - 1)  [since omega^{kn} = 1]

Denominator: (z-1)^2 = (omega^k - 1)^2

So sum j z^j = z * n omega^{-k} (omega^k - 1) / (omega^k - 1)^2
= omega^k * n * omega^{-k} / (omega^k - 1)
= n / (omega^k - 1)

Therefore lambda_k = n / (omega^k - 1) for k = 1,...,n-1.

THIS IS VALID FOR ALL n, not just odd n!

det(A) = det(E_std) / n... wait no. The relationship between A and E_std.

Actually det(E_std) = 0 (singular, row sums = 0).
The standard identity: det(A) where A is the (n-1)x(n-1) "reduced" matrix,
satisfies n * det(A) = det(E_std restricted to V_std)... hmm.

Actually: E_std has eigenvalues lambda_0 = 0 (k=0) and 
lambda_k = n/(omega^k - 1) for k = 1,...,n-1.

The product of nonzero eigenvalues = prod_{k=1}^{n-1} n/(omega^k - 1)
= n^{n-1} / prod_{k=1}^{n-1} (omega^k - 1).

By the cyclotomic identity: prod_{k=1}^{n-1} (x - omega^k) = (x^n-1)/(x-1)
Evaluating at x=1: prod_{k=1}^{n-1} (1 - omega^k) = n.
So: prod_{k=1}^{n-1} (omega^k - 1) = (-1)^{n-1} * prod (1-omega^k) = (-1)^{n-1} * n.

Product of nonzero eigenvalues = n^{n-1} / ((-1)^{n-1} * n) = (-1)^{n-1} * n^{n-2}.

Now: what's the relation to det(A)?
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

# The relation between det(A) and eigenvalues of E_std:
# E_std has eigenvalue 0 with eigenvector (1,...,1).
# The other eigenvalues are lambda_1,...,lambda_{n-1}.
# 
# The relation: det(A) = (1/n) * prod_{k=1}^{n-1} lambda_k ?
# No: det(E_std) = n * det(A) (well-known for centered matrix of LS).
# But det(E_std) = 0. Hmm.
#
# Actually: E_std = G * A where G is the n×(n-1) matrix.
# No. Let me think again.
#
# The standard relationship: if L is n×n with constant row sums s,
# and A = "reduced" (n-1)×(n-1) by deleting last column and row:
# Actually A_{ij} = L_{ij} - L_{i,n-1} for i,j = 0,...,n-2.
#
# Then det(A) is related to the cofactor of E_std... not simply.
#
# For the CYCLIC LS: E_std is circulant, so its eigenvalues on the 
# (n-1)-dimensional subspace V_std = {v : sum v_i = 0} are exactly
# lambda_k = n/(omega^k - 1) for k = 1,...,n-1.
#
# The matrix A acts on V_std via conjugation. If U is the DFT matrix 
# restricted to V_std, then A = U^{-1} diag(lambda_1,...,lambda_{n-1}) U.
# Wait no, A is (n-1)×(n-1) and acts on R^{n-1}, not V_std ⊂ R^n.
#
# There's a known identity: for a matrix with row sums = s,
# the cofactor = det of (n-1)×(n-1) minor = sum of all n cofactors / n = ...
# 
# For a circulant C with eigenvalues mu_0, mu_1, ..., mu_{n-1}:
# det(C) = prod mu_k.
# The (i,j) cofactor = (1/n) * prod_{k ≠ 0} mu_k * omega^{(i-j)k} ... no.
# 
# For our case: the matrix E is n×n circulant with one zero eigenvalue.
# All cofactors of E are equal (for circulant with simple zero eigenvalue).
# Each cofactor = prod_{k=1}^{n-1} lambda_k / n ... actually:
# 
# By the matrix-tree theorem / cofactor formula:
# If E has eigenvalues 0, lambda_1, ..., lambda_{n-1}:
# Each cofactor = (1/n) * prod_{k=1}^{n-1} lambda_k.
#
# Now: which cofactor relates to det(A)?
# A_{ij} = L_{ij} - L_{i,n-1} = E_std_{ij} - E_std_{i,n-1}
# The transformation L → A is: take columns 0..n-2, subtract column n-1.
# In matrix form: A = E_std[:n-1, :n-1] - E_std[:n-1, n-1] @ ones^T
# Hmm, this isn't simply a cofactor.
#
# Actually for ANY n×n matrix M with row sums = 0:
# Adding all other columns to the last column gives a zero column.
# So rank(M) ≤ n-1. The (n-1)×(n-1) "reduced" matrix obtained by 
# deleting the last row and column has:
# det(reduced) = (-1)^{n-1} * cofactor(n-1, n-1) of M.
#
# But our A is obtained differently: A_ij = M[i][j] - M[i][n-1].
# This is: taking columns 0,...,n-2 and subtracting column n-1 from each,
# then taking the first n-1 rows.
#
# Column operations don't change the absolute value of the determinant 
# when they're elementary (add multiple of one column to another).
# Subtracting column n-1 from columns 0,...,n-2 gives:
# New col j = old col j - old col n-1, for j = 0,...,n-2.
# The new column n-1 is unchanged.
# The det of this n×n matrix = original det (elementary column ops, det unchanged).
# 
# But the last column is unchanged, and then we take the (n-1)×(n-1) upper-left block.
# If M has row sums = 0, then after the column subtraction:
# New row sum = sum_{j=0}^{n-2} (M_{ij} - M_{i,n-1}) + M_{i,n-1}
#             = sum_{j=0}^{n-2} M_{ij} - (n-1)M_{i,n-1} + M_{i,n-1}
#             = sum_j M_{ij} - (n-1)M_{i,n-1} + M_{i,n-1} - M_{i,n-1}
# Hmm, this is getting complicated. Let me just verify the relation numerically.

print("Verify: det(A) vs product of nonzero eigenvalues of E_std")
print("=" * 60)

for n in range(2, 13):
    L = cyclic_ls(n)
    A = get_A(L)
    det_A = bareiss(A.tolist())
    
    # E_std as n×n circulant
    E = np.array(L, dtype=float) - (n+1)/2
    eigs = np.linalg.eigvals(E)
    # Sort by magnitude (zero eigenvalue should be smallest)
    eigs_sorted = sorted(eigs, key=lambda x: abs(x))
    nonzero_eigs = eigs_sorted[1:]  # skip the zero eigenvalue
    prod_nonzero = np.prod(nonzero_eigs)
    
    # Test: is det(A) = (1/n) * prod_nonzero?  
    ratio = det_A / prod_nonzero.real if abs(prod_nonzero) > 1e-6 else float('nan')
    
    print(f"  n={n:2d}: det(A)={det_A:>14d}  prod_nonzero={prod_nonzero.real:>16.1f}  "
          f"ratio={ratio:>8.4f}")

# For circulant with row sums 0, the cofactor formula gives:
# cofactor(i,j) = (1/n) * prod nonzero eigenvalues * omega^{...}
# For i=j=n-1: cofactor(n-1,n-1) = (1/n) * prod nonzero eigenvalues
# And det(A) should relate to this cofactor somehow.

print()
print("The relationship: det(A) = ? × (1/n) × prod nonzero eigenvalues")
print("Let's check the sign pattern:")
for n in range(2, 13):
    L = cyclic_ls(n)
    E = np.array(L, dtype=float) - (n+1)/2
    # Exact product of nonzero eigenvalues
    # lambda_k = n / (omega^k - 1) for k=1,...,n-1
    omega = np.exp(2j * np.pi / n)
    eigenvals = [n / (omega**k - 1) for k in range(1, n)]
    prod_eig = 1
    for e in eigenvals:
        prod_eig *= e
    
    A = get_A(L)
    det_A = bareiss(A.tolist())
    
    # cofactor(n-1,n-1) of E = (-1)^{2(n-1)} det(top-left (n-1)x(n-1) block of E)
    E_reduced = E[:n-1, :n-1]
    det_E_red = np.linalg.det(E_reduced)
    
    # Our det(A) 
    print(f"  n={n}: det(A)={det_A:>14d}  "
          f"det(E[:n-1,:n-1])={det_E_red:>14.1f}  "
          f"prod_eig/n = {(prod_eig/n).real:>14.1f}  "
          f"|prod_eig/n|={abs(prod_eig/n):>14.1f}")
