"""
snf_explore.py — Smith Normal Form analysis of A for Latin squares.

Goal: Understand the elementary divisor structure to prove Conjecture 28.11.B:
  v2(det A) >= v2(n/2) + max(0, kdim - 1)

Strategy:
  1. Compute SNF of A over Z (via sympy) for n=4 (exhaustive) and n=6,8 (JM samples)
  2. Extract invariant factors d_1 | d_2 | ... | d_{n-1}
  3. Analyze v2(d_i) vs kdim = dim ker_F2(B)
  4. Look for the pattern: kdim extra 0s in (d_i mod 2) <=> kdim elementary divisors even
"""
import sys as _sys
try:
    _sys.stdout.reconfigure(encoding="utf-8")
    _sys.stderr.reconfigure(encoding="utf-8")
except Exception:
    pass

import numpy as np
import random, time, sys
from collections import Counter, defaultdict

# ── Helpers (from jm_sharp_bound.py) ──────────────────────────────

def kernel_dim_F2(M_int):
    B = np.array(M_int % 2, dtype=int)
    rows, cols = B.shape; pr = 0
    for col in range(cols):
        found = False
        for row in range(pr, rows):
            if B[row, col] == 1:
                B[[pr, row]] = B[[row, pr]]; found = True; break
        if not found: continue
        for row in range(rows):
            if row != pr and B[row, col] == 1:
                B[row] = (B[row] + B[pr]) % 2
        pr += 1
    return cols - pr

def v2(x):
    if x == 0: return 999
    x = abs(int(x)); v = 0
    while x % 2 == 0: x //= 2; v += 1
    return v

def get_A(L):
    n = L.shape[0]
    return np.array([[L[i,j]-L[i,n-1] for j in range(n-1)] for i in range(n-1)], dtype=int)

def jm_sample(n, steps):
    cube = [[[0]*n for _ in range(n)] for _ in range(n)]
    for i in range(n):
        for j in range(n):
            cube[i][j][(i+j)%n] = 1
    proper = True; ir = ic = ik = 0; count = 0
    while count < steps or not proper:
        if proper:
            i = random.randint(0,n-1); j = random.randint(0,n-1)
            row_ij = cube[i][j]; k_cur = 0
            while row_ij[k_cur] != 1: k_cur += 1
            k = random.randint(0,n-2)
            if k >= k_cur: k += 1
            count += 1
        else:
            i,j,k = ir,ic,ik
            cands = [x for x in range(n) if cube[i][j][x]==1]
            k_cur = cands[random.randint(0,len(cands)-1)]
        jc = [jj for jj in range(n) if jj!=j and cube[i][jj][k]==1]
        j_p = jc[random.randint(0,len(jc)-1)]
        ic_ = [ii for ii in range(n) if ii!=i and cube[ii][j][k]==1]
        i_p = ic_[random.randint(0,len(ic_)-1)]
        cube[i][j][k]+=1; cube[i][j][k_cur]-=1
        cube[i][j_p][k]-=1; cube[i][j_p][k_cur]+=1
        cube[i_p][j][k]-=1; cube[i_p][j][k_cur]+=1
        cube[i_p][j_p][k]+=1; cube[i_p][j_p][k_cur]-=1
        if cube[i_p][j_p][k_cur]==-1:
            proper=False; ir,ic,ik = i_p,j_p,k_cur
        else:
            proper=True
    L = np.zeros((n,n),dtype=int)
    for i in range(n):
        for j in range(n):
            for k in range(n):
                if cube[i][j][k]==1: L[i,j]=k+1; break
    return L

def all_ls_4():
    results = []
    L = [[0]*4 for _ in range(4)]
    def ok(r,c,v):
        for j2 in range(c):
            if L[r][j2]==v: return False
        for i2 in range(r):
            if L[i2][c]==v: return False
        return True
    def solve(pos):
        if pos==16:
            results.append([row[:] for row in L]); return
        r,c = pos//4, pos%4
        for v in range(1,5):
            if ok(r,c,v): L[r][c]=v; solve(pos+1); L[r][c]=0
    solve(0)
    return results

def snf_invariants(A_np):
    """Compute Smith Normal Form of integer matrix A via row/col operations.
    Returns list of invariant factors d_1 | d_2 | ... | d_m (all >= 0)."""
    M = [[int(A_np[i,j]) for j in range(A_np.shape[1])] for i in range(A_np.shape[0])]
    rows, cols = len(M), len(M[0])
    m = min(rows, cols)

    for k in range(m):
        # Reduce to make M[k][k] the gcd of the submatrix M[k:,k:]
        changed = True
        while changed:
            changed = False
            # Find nonzero pivot in submatrix
            pivot = None
            for i in range(k, rows):
                for j in range(k, cols):
                    if M[i][j] != 0:
                        if pivot is None or abs(M[i][j]) < abs(M[pivot[0]][pivot[1]]):
                            pivot = (i, j)
            if pivot is None:
                break
            # Swap pivot to (k,k)
            pi, pj = pivot
            if pi != k:
                M[k], M[pi] = M[pi], M[k]
            if pj != k:
                for i in range(rows):
                    M[i][k], M[i][pj] = M[i][pj], M[i][k]

            # Make M[k][k] positive
            if M[k][k] < 0:
                for j in range(cols):
                    M[k][j] = -M[k][j]

            # Eliminate column k below row k
            for i in range(k+1, rows):
                if M[i][k] != 0:
                    q = M[i][k] // M[k][k]
                    for j in range(cols):
                        M[i][j] -= q * M[k][j]
                    if M[i][k] != 0:
                        changed = True

            # Eliminate row k right of column k
            for j in range(k+1, cols):
                if M[k][j] != 0:
                    q = M[k][j] // M[k][k]
                    for i in range(rows):
                        M[i][j] -= q * M[i][k]
                    if M[k][j] != 0:
                        changed = True

            # Check divisibility: M[k][k] must divide all entries in submatrix
            divisibility_ok = True
            for i in range(k+1, rows):
                for j in range(k+1, cols):
                    if M[i][j] % M[k][k] != 0:
                        # Add row i to row k to break divisibility
                        for jj in range(cols):
                            M[k][jj] += M[i][jj]
                        changed = True
                        divisibility_ok = False
                        break
                if not divisibility_ok:
                    break

    return [abs(M[i][i]) if i < rows and i < cols else 0 for i in range(m)]

# ══════════════════════════════════════════════════════════════════════
random.seed(271828)
print("=" * 70)
print("  SNF EXPLORATION — Elementary Divisors vs Kernel Dimension")
print("=" * 70)

# ── Part 1: n=4 exhaustive ───────────────────────────────────────────
print("\n>>> PART 1: n=4 EXHAUSTIVE — SNF of A")
all4 = all_ls_4()
snf_data_4 = defaultdict(list)  # kdim -> list of (invariants, v2_det)

for ls in all4:
    L = np.array(ls, dtype=int); A = get_A(L)
    det_A = int(np.round(np.linalg.det(A.astype(float))))
    if det_A == 0: continue  # skip singular
    kdim = kernel_dim_F2(A)
    try:
        invs = snf_invariants(A)
        snf_data_4[kdim].append((invs, v2(det_A)))
    except Exception as e:
        print(f"  SNF error: {e}")

for kd in sorted(snf_data_4):
    entries = snf_data_4[kd]
    print(f"\n  kdim={kd}: {len(entries)} Latin squares")
    # Show distribution of v2 of each invariant factor
    n_inv = len(entries[0][0])
    for j in range(n_inv):
        v2_vals = Counter(v2(e[0][j]) for e in entries)
        print(f"    d_{j+1}: v2 distribution = {dict(sorted(v2_vals.items()))}")
    # Show a few examples
    for ex in entries[:3]:
        invs, vd = ex
        print(f"    example: factors={invs}, v2(det)={vd}, v2(factors)={[v2(d) for d in invs]}")

sys.stdout.flush()

# ── Part 2: n=6 JM (100 samples) ────────────────────────────────────
print("\n>>> PART 2: n=6 JM — SNF of A (100 samples)")
t0 = time.time()
snf_data_6 = defaultdict(list)

for trial in range(100):
    L = jm_sample(6, 2000); A = get_A(L)
    det_A = int(np.round(np.linalg.det(A.astype(float))))
    if det_A == 0: continue
    kdim = kernel_dim_F2(A)
    try:
        invs = snf_invariants(A)
        snf_data_6[kdim].append((invs, v2(det_A)))
    except Exception as e:
        print(f"  SNF error at trial {trial}: {e}")

for kd in sorted(snf_data_6):
    entries = snf_data_6[kd]
    print(f"\n  kdim={kd}: {len(entries)} samples")
    n_inv = len(entries[0][0])
    for j in range(n_inv):
        v2_vals = Counter(v2(e[0][j]) for e in entries)
        print(f"    d_{j+1}: v2 distribution = {dict(sorted(v2_vals.items()))}")
    for ex in entries[:2]:
        invs, vd = ex
        print(f"    example: factors={invs}, v2(det)={vd}, v2(factors)={[v2(d) for d in invs]}")

print(f"  [{time.time()-t0:.1f}s]")
sys.stdout.flush()

# ── Part 3: n=8 JM (80 samples — SNF is slow for 7x7) ──────────────
print("\n>>> PART 3: n=8 JM — SNF of A (80 samples)")
t0 = time.time()
snf_data_8 = defaultdict(list)

for trial in range(80):
    L = jm_sample(8, 5000); A = get_A(L)
    det_A = int(np.round(np.linalg.det(A.astype(float))))
    if det_A == 0: continue
    kdim = kernel_dim_F2(A)
    try:
        invs = snf_invariants(A)
        snf_data_8[kdim].append((invs, v2(det_A)))
    except Exception as e:
        print(f"  SNF error at trial {trial}: {e}")
    if (trial+1) % 20 == 0:
        print(f"    {trial+1}/80...")
        sys.stdout.flush()

for kd in sorted(snf_data_8):
    entries = snf_data_8[kd]
    print(f"\n  kdim={kd}: {len(entries)} samples")
    n_inv = len(entries[0][0])
    for j in range(n_inv):
        v2_vals = Counter(v2(e[0][j]) for e in entries)
        print(f"    d_{j+1}: v2 distribution = {dict(sorted(v2_vals.items()))}")
    # Show a few kdim-specific examples
    for ex in entries[:2]:
        invs, vd = ex
        print(f"    example: factors={invs}, v2(det)={vd}, v2(factors)={[v2(d) for d in invs]}")

print(f"  [{time.time()-t0:.1f}s]")
sys.stdout.flush()

# ── Part 4: Analysis summary ────────────────────────────────────────
print("\n" + "=" * 70)
print("  SUMMARY — Key Patterns in SNF")
print("=" * 70)

for label, data, nn in [("n=4", snf_data_4, 4), ("n=6", snf_data_6, 6), ("n=8", snf_data_8, 8)]:
    print(f"\n  {label}: v2(n/2)={v2(nn//2)}")
    for kd in sorted(data):
        entries = data[kd]
        n_inv = len(entries[0][0])
        # Count how many invariant factors are even (per position)
        even_counts = []
        for j in range(n_inv):
            n_even = sum(1 for e in entries if e[0][j] % 2 == 0)
            even_counts.append(n_even)
        print(f"    kdim={kd}: #even_d_j = {even_counts} (out of {len(entries)})")
        # How many have exactly k even invariant factors?
        n_even_dist = Counter(sum(1 for d in e[0] if d % 2 == 0) for e in entries)
        print(f"    kdim={kd}: #even_total distribution = {dict(sorted(n_even_dist.items()))}")

print(f"\n  Total time: {time.time() - t0:.1f}s")
print("  DONE")
