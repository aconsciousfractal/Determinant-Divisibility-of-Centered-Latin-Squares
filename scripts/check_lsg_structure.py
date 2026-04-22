"""
Check if K(LSG) structure (not just order) is really constant.
Also check if the graphs are actually isomorphic or just co-spectral.
"""
import sys as _sys
try:
    _sys.stdout.reconfigure(encoding="utf-8")
    _sys.stderr.reconfigure(encoding="utf-8")
except Exception:
    pass
import numpy as np
from math import gcd

def get_A(L):
    n = L.shape[0]
    return np.array([[int(L[i,j]) - int(L[i,n-1]) for j in range(n-1)] for i in range(n-1)])

def lsg_laplacian(L):
    n = L.shape[0]; N = n*n
    adj = np.zeros((N,N), dtype=int)
    for i1 in range(n):
        for j1 in range(n):
            v1 = i1*n+j1
            for i2 in range(n):
                for j2 in range(n):
                    if (i1,j1)==(i2,j2): continue
                    if i1==i2 or j1==j2 or L[i1,j1]==L[i2,j2]:
                        adj[v1,i2*n+j2] = 1
    deg = adj.sum(axis=1)
    return np.diag(deg) - adj

def snf_diag(M):
    """Smith Normal Form — return diagonal."""
    A = [[int(M[i,j]) for j in range(M.shape[1])] for i in range(M.shape[0])]
    rows, cols = len(A), len(A[0]); m = min(rows, cols)
    for k in range(m):
        for _ in range(2000):
            best = None; best_val = None
            for i in range(k,rows):
                for j in range(k,cols):
                    if A[i][j] != 0:
                        if best_val is None or abs(A[i][j]) < best_val:
                            best=(i,j); best_val=abs(A[i][j])
            if best is None: break
            bi,bj = best
            if bi != k: A[k],A[bi] = A[bi],A[k]
            if bj != k:
                for i in range(rows): A[i][k],A[i][bj] = A[i][bj],A[i][k]
            if A[k][k] < 0:
                for j in range(cols): A[k][j] = -A[k][j]
            changed = False
            for i in range(k+1,rows):
                if A[i][k] != 0:
                    q = A[i][k]//A[k][k]
                    for j in range(cols): A[i][j] -= q*A[k][j]
                    if A[i][k] != 0: changed = True
            for j in range(k+1,cols):
                if A[k][j] != 0:
                    q = A[k][j]//A[k][k]
                    for i in range(rows): A[i][j] -= q*A[i][k]
                    if A[k][j] != 0: changed = True
            if not changed:
                ok = True
                for i in range(k+1,rows):
                    for j in range(k+1,cols):
                        if A[i][j] != 0 and A[i][j]%A[k][k] != 0:
                            for jj in range(cols): A[k][jj] += A[i][jj]
                            changed = True; ok = False; break
                    if not ok: break
                if ok: break
    return [abs(A[i][i]) for i in range(m)]

def generate_reduced_ls(n):
    results = []
    def backtrack(grid, row, col):
        if row == n:
            results.append([r[:] for r in grid]); return
        next_row, next_col = (row, col+1) if col+1 < n else (row+1, 1)
        used_row = set(grid[row])
        used_col = {grid[r][col] for r in range(row)}
        for val in range(1, n+1):
            if val not in used_row and val not in used_col:
                grid[row][col] = val
                backtrack(grid, next_row, next_col)
                grid[row][col] = 0
    grid = [[0]*n for _ in range(n)]
    for j in range(n): grid[0][j] = j + 1
    for i in range(n): grid[i][0] = i + 1
    backtrack(grid, 1, 1)
    return results

# n=4: all 4 reduced LS  
print("n=4: K(LSG) for all reduced LS")
ls4 = generate_reduced_ls(4)
print(f"  {len(ls4)} reduced LS")

for idx, ls in enumerate(ls4):
    L = np.array(ls)
    Lap = lsg_laplacian(L)
    L0 = Lap[:-1, :-1]
    s = snf_diag(L0)
    nontrivial = [d for d in s if d > 1]
    order = 1
    for d in nontrivial: order *= d
    print(f"  LS{idx}: K(LSG) = {nontrivial}, |K| = {order}")

# n=5: sample a few
print("\nn=5: K(LSG) for first 5 reduced LS")
ls5 = generate_reduced_ls(5)
print(f"  {len(ls5)} reduced LS total")

snf_profiles = set()
for idx in range(min(5, len(ls5))):
    L = np.array(ls5[idx])
    Lap = lsg_laplacian(L)
    L0 = Lap[:-1, :-1]
    s = snf_diag(L0)
    nontrivial = tuple(d for d in s if d > 1)
    snf_profiles.add(nontrivial)
    order = 1
    for d in nontrivial: order *= d
    print(f"  LS{idx}: |K(LSG)| = {order}, nontrivial factors count = {len(nontrivial)}")

if len(snf_profiles) == 1:
    print("  → K(LSG) SNF is IDENTICAL across all tested LS ✓")
else:
    print(f"  → K(LSG) has {len(snf_profiles)} DISTINCT SNF profiles!")
    for p in snf_profiles:
        print(f"    {p}")

# Key question: are the GRAPHS isomorphic?
# A sufficient test: check if adjacency matrices have same spectrum
print("\nn=4: Checking graph spectra")
spectra = set()
for idx, ls in enumerate(ls4):
    L = np.array(ls)
    Lap = lsg_laplacian(L)
    # Adjacency = diag(deg) - Lap... actually Lap = D - Adj
    deg = Lap[0,0]
    Adj = np.diag(np.diag(Lap)) - Lap
    eigs = sorted(np.round(np.linalg.eigvalsh(Adj.astype(float)), 6))
    spec_key = tuple(eigs)
    spectra.add(spec_key)
    print(f"  LS{idx}: spectrum eigenvalues (sorted): min={eigs[0]:.1f}, max={eigs[-1]:.1f}")

if len(spectra) == 1:
    print("  → All LSG(4) are CO-SPECTRAL (same adjacency spectrum)")
else:
    print(f"  → {len(spectra)} distinct spectra found!")

# The LSG is a strongly regular graph with parameters (n^2, 3(n-1), n-2, 6)
# All SRGs with same parameters are co-spectral (same eigenvalues)
# This means |K(LSG)| = product of nonzero Laplacian eigenvalues / n^2
# is determined by the SRG parameters alone.
print(f"\n--- CONCLUSION ---")
print(f"LSG(n) is an SRG({4}^2, 3·{3}, {2}, 6) = SRG(16, 9, 4, 6) for n=4")
print(f"All SRGs with same parameters have same Laplacian eigenvalues")
print(f"→ |K(LSG)| is constant (product of eigenvalues)")
print(f"→ But the GROUP STRUCTURE K(LSG) might differ (same order, different SNF)")
print(f"→ And the graphs are NOT necessarily isomorphic!")
