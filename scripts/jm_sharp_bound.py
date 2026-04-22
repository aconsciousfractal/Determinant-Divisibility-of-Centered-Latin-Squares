"""
jm_sharp_bound.py — Test the sharp bound conjecture:
  v2(det A) >= v2(n/2) + max(0, kdim - 1)
where kdim = dim ker_F2(A mod 2).

Also: prove kdim >= 2 => n | det(A) via adj(A) argument.
"""
import sys as _sys
try:
    _sys.stdout.reconfigure(encoding="utf-8")
    _sys.stderr.reconfigure(encoding="utf-8")
except Exception:
    pass

import numpy as np
from collections import Counter
import random, time, sys

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

def bareiss(M):
    n = M.shape[0]; A = M.astype(object); sign = 1; prev = 1
    for i in range(n):
        if A[i, i] == 0:
            found = False
            for k in range(i + 1, n):
                if A[k, i] != 0:
                    A[[i, k]] = A[[k, i]]; sign *= -1; found = True; break
            if not found: return 0
        for j in range(i + 1, n):
            for k in range(i + 1, n):
                A[j, k] = (A[j, k] * A[i, i] - A[j, i] * A[i, k]) // prev
            A[j, i] = 0
        prev = A[i, i]
    return sign * A[n - 1, n - 1]

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

# ══════════════════════════════════════════════════════════════════════
random.seed(314159)
print("=" * 70)
print("  SHARP BOUND:  v2(det) >= v2(n/2) + max(0, kdim-1)")
print("=" * 70)
T0 = time.time()

def test_bound(n, data):
    """Test v2(det) >= v2(n/2) + max(0, kdim-1) and report min v2 per kdim."""
    v2h = v2(n // 2)
    min_v2_per_kd = {}
    violations = 0
    for kd, vv in data:
        if vv == 999: continue  # det=0
        bound = v2h + max(0, kd - 1)
        if vv < bound:
            violations += 1
        if kd not in min_v2_per_kd or vv < min_v2_per_kd[kd]:
            min_v2_per_kd[kd] = vv
    print(f"  v2(n/2) = v2({n//2}) = {v2h}")
    print(f"  Violations: {violations}")
    for kd in sorted(min_v2_per_kd):
        bound = v2h + max(0, kd - 1)
        mv = min_v2_per_kd[kd]
        tight = "TIGHT!" if mv == bound else f"gap={mv - bound}"
        print(f"    kdim={kd}: min_v2={mv}, bound={bound} [{tight}]")

# ── n=4 exhaustive ──────────────────────────────────────────────────
print("\n>>> n=4 EXHAUSTIVE")
all4 = all_ls_4()
data4 = []
for ls in all4:
    L = np.array(ls, dtype=int); A = get_A(L)
    data4.append((kernel_dim_F2(A), v2(bareiss(A))))
test_bound(4, data4)
sys.stdout.flush()

# ── n=6 JM ──────────────────────────────────────────────────────────
print("\n>>> n=6 JM (500 samples, 2000 steps)")
t0 = time.time(); data6 = []
for _ in range(500):
    L = jm_sample(6, 2000); A = get_A(L)
    data6.append((kernel_dim_F2(A), v2(bareiss(A))))
test_bound(6, data6)
print(f"  [{time.time()-t0:.1f}s]"); sys.stdout.flush()

# ── n=8 JM ──────────────────────────────────────────────────────────
print("\n>>> n=8 JM (500 samples, 5000 steps)")
t0 = time.time(); data8 = []
for trial in range(500):
    L = jm_sample(8, 5000); A = get_A(L)
    data8.append((kernel_dim_F2(A), v2(bareiss(A))))
    if (trial+1)%100==0: print(f"    {trial+1}/500..."); sys.stdout.flush()
test_bound(8, data8)
print(f"  [{time.time()-t0:.1f}s]"); sys.stdout.flush()

# ── n=10 JM ─────────────────────────────────────────────────────────
print("\n>>> n=10 JM (200 samples, 3000 steps)")
t0 = time.time(); data10 = []
for _ in range(200):
    L = jm_sample(10, 3000); A = get_A(L)
    data10.append((kernel_dim_F2(A), v2(bareiss(A))))
test_bound(10, data10)
print(f"  [{time.time()-t0:.1f}s]"); sys.stdout.flush()

# ── n=12 JM ─────────────────────────────────────────────────────────
print("\n>>> n=12 JM (200 samples, 5000 steps)")
t0 = time.time(); data12 = []
for trial in range(200):
    L = jm_sample(12, 5000); A = get_A(L)
    data12.append((kernel_dim_F2(A), v2(bareiss(A))))
    if (trial+1)%50==0: print(f"    {trial+1}/200..."); sys.stdout.flush()
test_bound(12, data12)
print(f"  [{time.time()-t0:.1f}s]"); sys.stdout.flush()

# ── n=16 JM (if time permits) ───────────────────────────────────────
elapsed = time.time() - T0
if elapsed < 300:
    print("\n>>> n=16 JM (50 samples, 8000 steps)")
    t0 = time.time(); data16 = []
    for trial in range(50):
        L = jm_sample(16, 8000); A = get_A(L)
        data16.append((kernel_dim_F2(A), v2(bareiss(A))))
        if (trial+1)%10==0: print(f"    {trial+1}/50..."); sys.stdout.flush()
    test_bound(16, data16)
    print(f"  [{time.time()-t0:.1f}s]"); sys.stdout.flush()

# ── PROOF CHECK: adj(A) ≡ 0 mod 2 when kdim >= 2 ────────────────────
print("\n>>> PROOF CHECK: adj(A) mod 2 = 0 when kdim >= 2")
print("    (Verifying that all (n-2)x(n-2) minors of B are 0 over F2)")
checks = 0; fails = 0
for trial in range(100):
    L = jm_sample(8, 5000); A = get_A(L)
    kd = kernel_dim_F2(A)
    if kd < 2: continue
    B = A % 2
    n1 = B.shape[0]  # 7
    adj_zero = True
    for i in range(n1):
        for j in range(n1):
            # (n-2)x(n-2) minor: remove row j, col i from B
            rows = [r for r in range(n1) if r != j]
            cols = [c for c in range(n1) if c != i]
            sub = B[np.ix_(rows, cols)]
            # det over F2
            det_f2 = kernel_dim_F2(sub)  # if det=0, ker dim >= 1
            # Actually need det directly. rank = n-2 iff det != 0.
            rank_sub = (n1-1) - det_f2  # n1-1 = 6 for n=8
            if rank_sub == n1 - 1:  # full rank = n-2, det != 0
                adj_zero = False; break
        if not adj_zero: break
    checks += 1
    if not adj_zero: fails += 1
print(f"  Checked {checks} LS with kdim>=2: adj(B)=0 for {checks-fails}/{checks}")
sys.stdout.flush()

print(f"\n{'='*70}")
print(f"  TOTAL: {time.time()-T0:.1f}s")
print(f"{'='*70}")
