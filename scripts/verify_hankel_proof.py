"""
Verify the DIRECT Hankel-symmetric proof for Prop 3.6.

Key claim: E_{ij} = (i+j) mod n - (n-1)/2 is a real symmetric matrix.
Its eigenvalues on V_std can be computed from the DFT of h(j) = j-(n-1)/2.

Since E_{ij} = h((i+j) mod n) and E is real symmetric:
  E v^(k) = ĥ(k) v^(n-k)   for DFT vectors v^(k)_j = ω^{jk}

On the 2D real subspace {Re v^(k), Im v^(k)} (for k ≠ 0, n/2):
  Eigenvalues = ±|ĥ(k)| = ±n/|ω^k - 1|

For k = n/2 (even n only):
  Eigenvalue = ĥ(n/2) = -n/2

Product of |eigenvalues| = ∏|ĥ(k)| for k=1..n-1 = n^{n-2}
Sign: ⌊n/2⌋ negative eigenvalues → (-1)^{⌊n/2⌋}
"""
import numpy as np

def cyclic_ls(n):
    return [[(i+j) % n + 1 for j in range(n)] for i in range(n)]

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

print("=" * 70)
print("TEST 1: E is real symmetric and eigenvalues of A are all real")
print("=" * 70)
for n in range(2, 10):
    L = cyclic_ls(n)
    E = np.array([[((i+j) % n) - (n-1)/2 for j in range(n)] for i in range(n)])
    # Check symmetric
    assert np.allclose(E, E.T), f"E not symmetric at n={n}"
    # eigenvalues of A (= eigenvalues of E|_Vstd)
    A = np.array([[L[i][j] - L[i][n-1] for j in range(n-1)] for i in range(n-1)], dtype=float)
    eigs_A = np.linalg.eigvals(A)
    max_imag = max(abs(e.imag) for e in eigs_A)
    print(f"  n={n}: E symmetric ✓, max|Im(eigenval of A)|={max_imag:.2e}, "
          f"eigenvals(A) = {sorted([round(e.real,4) for e in eigs_A])}")

print()
print("=" * 70)
print("TEST 2: Ev^(k) = ĥ(k) * v^(n-k) for DFT vectors")
print("=" * 70)
for n in [3, 4, 5, 6, 7, 8]:
    E = np.array([[((i+j) % n) - (n-1)/2 for j in range(n)] for i in range(n)], dtype=complex)
    omega = np.exp(2j * np.pi / n)
    ok = True
    for k in range(1, n):
        vk = np.array([omega**(j*k) for j in range(n)])
        vnk = np.array([omega**(j*(n-k)) for j in range(n)])
        hhat_k = n / (omega**k - 1)
        Evk = E @ vk
        expected = hhat_k * vnk
        if not np.allclose(Evk, expected):
            ok = False
            print(f"  FAIL n={n}, k={k}")
    print(f"  n={n}: Ev^(k) = ĥ(k)·v^(n-k) for all k=1..{n-1}: {'✓' if ok else '✗'}")

print()
print("=" * 70)
print("TEST 3: Eigenvalues on 2D subspaces are ±|ĥ(k)|")
print("=" * 70)
for n in [3, 4, 5, 6, 7, 8]:
    E = np.array([[((i+j) % n) - (n-1)/2 for j in range(n)] for i in range(n)], dtype=float)
    omega = np.exp(2j * np.pi / n)
    
    # Eigenvalues of E (real symmetric → all real)
    eigs_E = sorted(np.linalg.eigvalsh(E))  # eigvalsh for symmetric
    nonzero_eigs = [e for e in eigs_E if abs(e) > 1e-10]
    
    # Predicted eigenvalues: for each pair (k, n-k), get ±|ĥ(k)|; for k=n/2, get ĥ(n/2)
    predicted = []
    visited = set()
    for k in range(1, n):
        if k in visited:
            continue
        hhat = n / (omega**k - 1)
        if k == n - k:  # k = n/2, only for even n
            predicted.append(hhat.real)  # should be -n/2
            visited.add(k)
        else:
            predicted.extend([abs(hhat), -abs(hhat)])
            visited.add(k)
            visited.add(n - k)
    
    predicted.sort()
    nonzero_eigs.sort()
    match = np.allclose(sorted(predicted), sorted(nonzero_eigs))
    print(f"  n={n}: predicted={[round(x,4) for x in sorted(predicted)]}")
    print(f"         actual  ={[round(x,4) for x in sorted(nonzero_eigs)]}")
    print(f"         match: {'✓' if match else '✗'}")

print()
print("=" * 70)
print("TEST 4: ∏|eigenvalues| = ∏|ĥ(k)| = n^{n-2}")
print("=" * 70)
for n in range(2, 15):
    omega = np.exp(2j * np.pi / n)
    prod_hhat = 1.0
    for k in range(1, n):
        prod_hhat *= abs(n / (omega**k - 1))
    
    L = cyclic_ls(n)
    A_int = [[L[i][j] - L[i][n-1] for j in range(n-1)] for i in range(n-1)]
    det_A = bareiss(A_int)
    
    expected = n ** (n - 2)
    print(f"  n={n:2d}: ∏|ĥ(k)|={prod_hhat:>16.1f}  n^(n-2)={expected:>16d}  "
          f"|det(A)|={abs(det_A):>16d}  match={'✓' if abs(prod_hhat - expected) < 0.5 and abs(det_A) == expected else '✗'}")

print()
print("=" * 70)
print("TEST 5: Sign count — number of negative eigenvalues = ⌊n/2⌋")
print("=" * 70)
for n in range(2, 15):
    E = np.array([[((i+j) % n) - (n-1)/2 for j in range(n)] for i in range(n)], dtype=float)
    eigs = np.linalg.eigvalsh(E)
    n_neg = sum(1 for e in eigs if e < -1e-10)
    expected_neg = n // 2
    
    L = cyclic_ls(n)
    A_int = [[L[i][j] - L[i][n-1] for j in range(n-1)] for i in range(n-1)]
    det_A = bareiss(A_int)
    predicted_det = ((-1) ** (n // 2)) * (n ** (n - 2))
    
    print(f"  n={n:2d}: #neg_eigs={n_neg}, ⌊n/2⌋={expected_neg}, "
          f"sign_match={'✓' if n_neg == expected_neg else '✗'}  "
          f"det_match={'✓' if det_A == predicted_det else '✗'}")

print()
print("=" * 70)
print("TEST 6: Cyclotomic identity ∏_{k=1}^{n-1}(1-ω^k) = n for ALL n")
print("=" * 70)
for n in range(2, 20):
    omega = np.exp(2j * np.pi / n)
    prod = 1.0 + 0j
    for k in range(1, n):
        prod *= (1 - omega**k)
    print(f"  n={n:2d}: ∏(1-ω^k) = {prod.real:>12.6f} + {prod.imag:>12.6f}i  "
          f"≈ {n}? {'✓' if abs(prod - n) < 1e-6 else '✗'}")

print()
print("=" * 70)
print("TEST 7: Why row permutation σ does NOT give A'=P_σ·A")
print("=" * 70)
for n in [3, 4, 5]:
    L = cyclic_ls(n)
    A = [[L[i][j] - L[i][n-1] for j in range(n-1)] for i in range(n-1)]
    
    # L' via σ: i -> (-i) mod n
    Lp = [[((-i+j) % n) + 1 for j in range(n)] for i in range(n)]
    Ap = [[Lp[i][j] - Lp[i][n-1] for j in range(n-1)] for i in range(n-1)]
    
    # σ maps {0,1,...,n-2} to {0,n-1,n-2,...,2}
    sigma = [(-i) % n for i in range(n-1)]
    print(f"  n={n}: σ({list(range(n-1))}) = {sigma}")
    print(f"         Row 1 of A' uses L-row σ(1)={sigma[1]}, "
          f"which is {'OUTSIDE' if sigma[1] == n-1 else 'inside'} A's row range [0,{n-2}]")
    print(f"         det(A)={bareiss(A)}, det(A')={bareiss(Ap)}, ratio={bareiss(A)/bareiss(Ap) if bareiss(Ap)!=0 else 'undef'}")
    print()
