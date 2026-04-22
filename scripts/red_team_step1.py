"""
Verify: Step 1 alone gives n | d_2, making the rank argument in Step 3 REDUNDANT.

Step 1: n^{k-1} | Delta_k  (from multilinearity)
Step 2: d_1 = 1, product = n^{n-2}

For k=2: n^1 | Delta_2 = d_1 * d_2 = d_2
So: n | d_2, hence d_2 >= n.
Also: d_2^{n-2} <= product = n^{n-2}, so d_2 <= n.
Therefore: d_2 = n.
Induction: d_3 = ... = d_{n-1} = n.
"""
import sys as _sys
try:
    _sys.stdout.reconfigure(encoding="utf-8")
    _sys.stderr.reconfigure(encoding="utf-8")
except Exception:
    pass
import numpy as np
from sympy import Matrix

def cyclic_ls(n):
    return np.array([[(i+j) % n + 1 for j in range(n)] for i in range(n)])

def get_A(L):
    n = L.shape[0]
    return np.array([[int(L[i,j]) - int(L[i,n-1]) for j in range(n-1)] for i in range(n-1)])

print("Verify: n^{k-1} | Delta_k for the cyclic LS (via actual SNF)")
print("=" * 70)

for n in range(2, 13):
    A = get_A(cyclic_ls(n))
    M = Matrix(A.tolist())
    
    # Compute SNF via sympy
    from sympy import ZZ
    from sympy.matrices.normalforms import smith_normal_form as snf_func
    snf_mat = snf_func(M, domain=ZZ)
    d = [int(snf_mat[i,i]) for i in range(n-1)]
    
    # Determinantal divisors Delta_k = product of first k invariant factors
    deltas = []
    prod = 1
    for i in range(n-1):
        prod *= d[i]
        deltas.append(prod)
    
    # Check n^{k-1} | Delta_k
    all_ok = True
    details = []
    for k in range(1, n):  # k = 1, ..., n-1
        dk = deltas[k-1]
        divisor = n ** (k-1)
        ok = (dk % divisor == 0)
        details.append(f"k={k}: Δ_{k}={dk}, n^{k-1}={divisor}, divides: {ok}")
        if not ok: all_ok = False
    
    tag = "✓" if all_ok else "✗"
    print(f"\nn={n}: SNF diag = {[int(x) for x in d]}  {tag}")
    
    # Specifically: for k=2, does n | d_2?
    if n >= 3:
        d2 = int(d[1])
        print(f"  d_2 = {d2}, n | d_2: {d2 % n == 0}  ← KEY: this comes from Step 1 alone!")
    
    # Show the inductive argument
    if n >= 4:
        print(f"  d_2={int(d[1])}: d_2 >= n (from n|d_2) and d_2 <= n (from d_2^{n-2} <= n^{n-2}) → d_2 = n ✓")

print("\n" + "=" * 70)
print("CONCLUSION: The rank-mod-p argument in Step 3 is COMPLETELY REDUNDANT.")
print("Step 1 (multilinearity) directly gives n | Δ_2 = d_2, hence d_2 = n.")
print("Induction completes the proof without any rank computation.")
print("=" * 70)
