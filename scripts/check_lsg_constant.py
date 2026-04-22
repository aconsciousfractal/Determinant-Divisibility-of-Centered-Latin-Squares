"""
Quick check: is K(LSG(n)) really constant across all LS of order n?
The paper claims "since all such graphs are isomorphic" — verify this.
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

def back_circ_ls(n):
    """Another LS: L_{i,j} = (i-j) mod n + 1"""
    return np.array([[(i-j) % n + 1 for j in range(n)] for i in range(n)])

def lsg_edges(L):
    """Return sorted edge set of LSG(L)."""
    n = L.shape[0]
    edges = set()
    for i1 in range(n):
        for j1 in range(n):
            for i2 in range(n):
                for j2 in range(n):
                    if (i1, j1) >= (i2, j2): continue
                    if i1 == i2 or j1 == j2 or L[i1,j1] == L[i2,j2]:
                        edges.add((i1*n+j1, i2*n+j2))
    return edges

def lsg_degree_seq(L):
    n = L.shape[0]; N = n*n
    deg = [0]*N
    for i1 in range(n):
        for j1 in range(n):
            for i2 in range(n):
                for j2 in range(n):
                    if (i1,j1) == (i2,j2): continue
                    if i1 == i2 or j1 == j2 or L[i1,j1] == L[i2,j2]:
                        deg[i1*n+j1] += 1
    return sorted(deg)

# For n=3, there's only 1 reduced LS, so trivially constant.
# For n=4, check:
print("n=4: Checking if all LS give isomorphic LSG...")
# LS 1: cyclic
L1 = cyclic_ls(4)
# LS 2: non-cyclic
L2 = np.array([
    [1, 2, 3, 4],
    [2, 1, 4, 3],
    [3, 4, 1, 2],
    [4, 3, 2, 1]
])
# LS 3: another
L3 = np.array([
    [1, 2, 3, 4],
    [2, 3, 4, 1],
    [3, 4, 1, 2],
    [4, 1, 2, 3]
])
# LS 4: yet another
L4 = np.array([
    [1, 2, 3, 4],
    [2, 4, 1, 3],
    [3, 1, 4, 2],
    [4, 3, 2, 1]
])

for i, L in enumerate([L1, L2, L3, L4], 1):
    ds = lsg_degree_seq(L)
    e = len(lsg_edges(L))
    print(f"  LS{i}: edges={e}, degree_seq={ds}")

print("\nIf degree sequences differ → graphs NOT isomorphic")
print("If degree sequences same → might be isomorphic (need more checks)")

# Actually, for Latin square graphs, each vertex (i,j) is adjacent to:
# - n-1 vertices in same row (same i)
# - n-1 vertices in same column (same j)
# - (multiplicity of L[i,j] in L) - 1 vertices with same symbol
# But actually each symbol appears exactly n times in L (once per row),
# so same-symbol gives n-1 more adjacencies.
# Total degree = (n-1) + (n-1) + (n-1) = 3(n-1) per vertex.
# But we might double-count: vertices (i,j) and (i,j') share row AND
# possibly same symbol. The adjacency is still 1 (unweighted).
# So degree = #{(i',j'): same row OR same col OR same symbol} = 
# Actually this overcounts. Let me compute properly.

print("\n--- Degree analysis ---")
for n in [4, 5]:
    L = cyclic_ls(n)
    N = n*n
    for v in range(min(3, N)):
        i, j = v // n, v % n
        neighbors = set()
        for i2 in range(n):
            for j2 in range(n):
                if (i2,j2) == (i,j): continue
                if i2 == i or j2 == j or L[i2,j2] == L[i,j]:
                    neighbors.add((i2,j2))
        print(f"  n={n}, vertex ({i},{j}): degree={len(neighbors)}")
