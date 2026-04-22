"""
jm_crosscheck.py — Unbiased cross-tabulation via Jacobson-Matthews sampler.

The JM chain is provably ergodic on the set of ALL Latin squares of order n.
This addresses the swap-chain sampling bias discovered in §28.10.
Standard symbols {1,...,n} are preserved (no symbol permutation).
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

# ── Helpers ──────────────────────────────────────────────────────────

def kernel_dim_F2(M_int):
    B = np.array(M_int % 2, dtype=int)
    rows, cols = B.shape
    pivot_row = 0
    for col in range(cols):
        found = False
        for row in range(pivot_row, rows):
            if B[row, col] == 1:
                B[[pivot_row, row]] = B[[row, pivot_row]]
                found = True
                break
        if not found:
            continue
        for row in range(rows):
            if row != pivot_row and B[row, col] == 1:
                B[row] = (B[row] + B[pivot_row]) % 2
        pivot_row += 1
    return cols - pivot_row

def bareiss(M):
    n = M.shape[0]
    A = M.astype(object)
    sign, prev = 1, 1
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
    return np.array([[L[i, j] - L[i, n-1] for j in range(n-1)]
                      for i in range(n-1)], dtype=int)

def is_latin_square(L):
    n = L.shape[0]; s = set(range(1, n+1))
    return all(set(L[i]) == s and set(L[:, i]) == s for i in range(n))

def one_in_ker(A):
    for i in range(A.shape[0]):
        if sum(int(A[i, j]) for j in range(A.shape[1])) % 2 != 0:
            return False
    return True

# ── Jacobson-Matthews ────────────────────────────────────────────────

def jm_sample(n, steps):
    """Uniform random LS of order n via Jacobson-Matthews chain."""
    cube = [[[0]*n for _ in range(n)] for _ in range(n)]
    for i in range(n):
        for j in range(n):
            cube[i][j][(i+j) % n] = 1

    proper = True; ir = ic = ik = 0
    count = 0
    while count < steps or not proper:
        if proper:
            i = random.randint(0, n-1); j = random.randint(0, n-1)
            row_ij = cube[i][j]
            k_cur = 0
            while row_ij[k_cur] != 1: k_cur += 1
            k = random.randint(0, n-2)
            if k >= k_cur: k += 1
            count += 1
        else:
            i, j, k = ir, ic, ik
            cands = [x for x in range(n) if cube[i][j][x] == 1]
            k_cur = cands[random.randint(0, len(cands)-1)]

        # j': col in row i with symbol k, j'!=j
        jc = [jj for jj in range(n) if jj != j and cube[i][jj][k] == 1]
        j_p = jc[random.randint(0, len(jc)-1)]
        # i': row in col j with symbol k, i'!=i
        ic_ = [ii for ii in range(n) if ii != i and cube[ii][j][k] == 1]
        i_p = ic_[random.randint(0, len(ic_)-1)]

        cube[i][j][k] += 1;     cube[i][j][k_cur] -= 1
        cube[i][j_p][k] -= 1;   cube[i][j_p][k_cur] += 1
        cube[i_p][j][k] -= 1;   cube[i_p][j][k_cur] += 1
        cube[i_p][j_p][k] += 1; cube[i_p][j_p][k_cur] -= 1

        if cube[i_p][j_p][k_cur] == -1:
            proper = False; ir, ic, ik = i_p, j_p, k_cur
        else:
            proper = True

    L = np.zeros((n, n), dtype=int)
    for i in range(n):
        for j in range(n):
            for k in range(n):
                if cube[i][j][k] == 1:
                    L[i, j] = k + 1; break
    return L

# ── n=4 exhaustive ──────────────────────────────────────────────────

def all_ls_4():
    results = []
    L = [[0]*4 for _ in range(4)]
    def ok(r, c, v):
        for j2 in range(c):
            if L[r][j2] == v: return False
        for i2 in range(r):
            if L[i2][c] == v: return False
        return True
    def solve(pos):
        if pos == 16:
            results.append([row[:] for row in L]); return
        r, c = pos // 4, pos % 4
        for v in range(1, 5):
            if ok(r, c, v):
                L[r][c] = v; solve(pos + 1); L[r][c] = 0
    solve(0)
    return results

# ══════════════════════════════════════════════════════════════════════
random.seed(42)
print("=" * 65)
print("  JACOBSON-MATTHEWS   UNBIASED   CROSS-TABULATION")
print("=" * 65)
T0 = time.time()

# ── PART 0: Validation ───────────────────────────────────────────────
print("\n>>> PART 0: JM Validation (n=8, 20 samples, 5000 steps)")
t0 = time.time()
ok_ct = 0; distinct = set()
for _ in range(20):
    L = jm_sample(8, 5000)
    if is_latin_square(L): ok_ct += 1
    distinct.add(tuple(L.flatten()))
print(f"  Valid: {ok_ct}/20  Distinct: {len(distinct)}/20  Time: {time.time()-t0:.1f}s")
sys.stdout.flush()

# ── PART 1: n=4 EXHAUSTIVE ──────────────────────────────────────────
print("\n>>> PART 1: n=4 EXHAUSTIVE (all Latin squares)")
t0 = time.time()
all4 = all_ls_4()
print(f"  Enumerated: {len(all4)} LS")

rows_4 = []
for ls in all4:
    L = np.array(ls, dtype=int); A = get_A(L); d = bareiss(A)
    B = A % 2; B2 = B @ B
    rows_4.append({
        'd': int(d), 'v2': v2(d),
        'kd': kernel_dim_F2(A),
        'kd2': kernel_dim_F2(B2),
        'n2': (int(d) % 16 == 0) if d != 0 else None,
        'ok1': one_in_ker(A),
    })

# Full cross-tab
ct = Counter((r['kd'], str(r['n2']), r['d']==0) for r in rows_4)
print("  Cross-tab (kdim, 16|det, det==0):")
for k in sorted(ct): print(f"    kd={k[0]} 16|d={k[1]} z={k[2]} -> {ct[k]}")

# Excluding det=0
nz = [r for r in rows_4 if r['d'] != 0]
print(f"\n  Non-zero det: {len(nz)}/{len(rows_4)}")
ct2 = Counter((r['kd'], r['n2']) for r in nz)
print("  kdim_B vs 16|det (det!=0):")
for k in sorted(ct2): print(f"    kd={k[0]} 16|d={k[1]} -> {ct2[k]}")

ct3 = Counter((r['kd2'], r['n2']) for r in nz)
print("  kdim_B^2 vs 16|det (det!=0):")
for k in sorted(ct3): print(f"    kd2={k[0]} 16|d={k[1]} -> {ct3[k]}")

oc = sum(1 for r in rows_4 if r['ok1'])
print(f"\n  1 in ker(B): {oc}/{len(rows_4)} ({100*oc/len(rows_4):.1f}%)")
v2d = Counter(r['v2'] for r in nz)
print(f"  v2(det) det!=0: {dict(sorted(v2d.items()))}")
d4 = sum(1 for r in nz if int(r['d']) % 4 == 0)
print(f"  4|det: {d4}/{len(nz)} ({100*d4/len(nz):.1f}%)")
# Check suff/nec for criterion kd>=2 <=> 16|det
suff_f = sum(1 for r in nz if r['kd'] >= 2 and not r['n2'])
nec_f = sum(1 for r in nz if r['kd'] < 2 and r['n2'])
print(f"  Criterion kd>=2 <=> 16|det:  suff_fail={suff_f}  nec_fail={nec_f}")
print(f"  Time: {time.time()-t0:.1f}s")
sys.stdout.flush()

# ── PART 2: n=8 JM ──────────────────────────────────────────────────
N, S, NS = 8, 300, 3000
print(f"\n>>> PART 2: n={N} JM ({S} samples, {NS} steps)")
t0 = time.time()
rows_8 = []
for trial in range(S):
    L = jm_sample(N, NS); A = get_A(L); d = bareiss(A)
    rows_8.append({'d': int(d), 'v2': v2(d), 'kd': kernel_dim_F2(A),
                    'n2': (int(d) % 64 == 0) if d != 0 else None,
                    'ok1': one_in_ker(A)})
    if (trial+1) % 100 == 0:
        print(f"    {trial+1}/{S} ..."); sys.stdout.flush()

nz8 = [r for r in rows_8 if r['d'] != 0]; z8 = len(rows_8) - len(nz8)
ct8 = Counter((r['kd'], r['n2']) for r in nz8)
print(f"  det=0: {z8}/{len(rows_8)}")
print("  kdim_B vs 64|det (det!=0):")
for k in sorted(ct8): print(f"    kd={k[0]} 64|d={k[1]} -> {ct8[k]}")

suff_f8 = sum(1 for r in nz8 if r['kd'] >= 2 and not r['n2'])
nec_f8 = sum(1 for r in nz8 if r['kd'] < 2 and r['n2'])
print(f"  Criterion kd>=2 <=> 64|det:  suff_fail={suff_f8}  nec_fail={nec_f8}")
oc8 = sum(1 for r in rows_8 if r['ok1'])
print(f"  1 in ker(B): {oc8}/{len(rows_8)} ({100*oc8/len(rows_8):.1f}%)")
v2d8 = Counter(r['v2'] for r in nz8)
print(f"  v2(det) det!=0: {dict(sorted(v2d8.items()))}")
nd8 = sum(1 for r in nz8 if int(r['d']) % 8 == 0)
print(f"  8|det: {nd8}/{len(nz8)} ({100*nd8/len(nz8):.1f}%)")
print(f"  Time: {time.time()-t0:.1f}s")
sys.stdout.flush()

# ── PART 3: n=12 JM ─────────────────────────────────────────────────
N, S, NS = 12, 150, 5000
print(f"\n>>> PART 3: n={N} JM ({S} samples, {NS} steps)")
t0 = time.time()
rows_12 = []
for trial in range(S):
    L = jm_sample(N, NS); A = get_A(L); d = bareiss(A)
    rows_12.append({'d': int(d), 'v2': v2(d), 'kd': kernel_dim_F2(A),
                     'n2': (int(d) % 144 == 0) if d != 0 else None,
                     'ok1': one_in_ker(A)})
    if (trial+1) % 50 == 0:
        print(f"    {trial+1}/{S} ..."); sys.stdout.flush()

nz12 = [r for r in rows_12 if r['d'] != 0]; z12 = len(rows_12) - len(nz12)
ct12 = Counter((r['kd'], r['n2']) for r in nz12)
print(f"  det=0: {z12}/{len(rows_12)}")
print("  kdim_B vs 144|det (det!=0):")
for k in sorted(ct12): print(f"    kd={k[0]} 144|d={k[1]} -> {ct12[k]}")

suff_f12 = sum(1 for r in nz12 if r['kd'] >= 2 and not r['n2'])
nec_f12 = sum(1 for r in nz12 if r['kd'] < 2 and r['n2'])
print(f"  Criterion kd>=2 <=> 144|det:  suff_fail={suff_f12}  nec_fail={nec_f12}")
oc12 = sum(1 for r in rows_12 if r['ok1'])
print(f"  1 in ker(B): {oc12}/{len(rows_12)} ({100*oc12/len(rows_12):.1f}%)")
v2d12 = Counter(r['v2'] for r in nz12)
print(f"  v2(det) det!=0: {dict(sorted(v2d12.items()))}")
nd12 = sum(1 for r in nz12 if int(r['d']) % 12 == 0)
print(f"  12|det: {nd12}/{len(nz12)} ({100*nd12/len(nz12):.1f}%)")
print(f"  Time: {time.time()-t0:.1f}s")
sys.stdout.flush()

# ── PART 4: Theorem 4.3 (n ≡ 2 mod 4) ──────────────────────────────
for n_t, samp, stp in [(6, 300, 2000), (10, 100, 3000)]:
    n2 = n_t * n_t
    print(f"\n>>> PART 4: Th4.3 n={n_t} ({samp} JM, {stp} steps)")
    t0 = time.time(); fail = 0; v2_min = 999; vd = Counter()
    for trial in range(samp):
        L = jm_sample(n_t, stp); A = get_A(L); d = bareiss(A)
        vv = v2(d); vd[vv] += 1
        if vv < v2_min: v2_min = vv
        if d != 0 and int(d) % n2 != 0:
            fail += 1
            if fail <= 3: print(f"    COUNTEREX: det={d} v2={vv}")
    print(f"  Failures: {fail}/{samp}  min_v2={v2_min}")
    print(f"  v2 dist: {dict(sorted(vd.items()))}")
    print(f"  Time: {time.time()-t0:.1f}s")
    sys.stdout.flush()

# ── PART 5: Base column invariance (n=8) ─────────────────────────────
print(f"\n>>> PART 5: Base column invariance (n=8, 10 JM samples)")
t0 = time.time(); ok5 = 0
for trial in range(10):
    L = jm_sample(8, 3000); dets = []
    for base in range(8):
        cols = [j for j in range(8) if j != base]
        Ab = np.array([[L[i, c] - L[i, base] for c in cols]
                        for i in range(7)], dtype=int)
        dets.append(abs(int(bareiss(Ab))))
    if len(set(dets)) == 1:
        ok5 += 1
    else:
        print(f"    Trial {trial}: DIFFER {dets}")
print(f"  All |det| same: {ok5}/10  Time: {time.time()-t0:.1f}s")

print(f"\n{'='*65}")
print(f"  TOTAL: {time.time()-T0:.1f}s")
print(f"{'='*65}")
