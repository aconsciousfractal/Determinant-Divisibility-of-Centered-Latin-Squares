"""
jm_paper_verify.py — Check the *actual* paper claims with unbiased JM sampling.

Paper's Theorem (Unified divisibility):
  (n^2 / gcd(n,2)) | det(E_std)   where det(E_std) = n * det(A)
  => For even n: (n/2) | det(A)         [UNIVERSAL]
  => For odd n:  n | det(A)             [UNIVERSAL]

Theorem (F2-rank characterization, n ≡ 2 mod 4):
  n^2 | det(E_std) <=> rank_F2(B) < n-1 <=> dim ker(B) >= 1 <=> 2 | det(A)
  (Since odd part of n always divides det(A), only prime 2 can obstruct.)

Open question (n ≡ 0 mod 4):
  1 in ker(B) always, so 2 | det(A) always.
  But n | det(A) requires 2^a | det(A) where a = v2(n) >= 2.
  Paper says "requires finer 2-adic analysis".
"""

import numpy as np
from collections import Counter
import random, time, sys

# ── Helpers (same as jm_crosscheck.py) ───────────────────────────────

def kernel_dim_F2(M_int):
    B = np.array(M_int % 2, dtype=int)
    rows, cols = B.shape
    pr = 0
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
    return np.array([[L[i,j] - L[i,n-1] for j in range(n-1)]
                      for i in range(n-1)], dtype=int)

def jm_sample(n, steps):
    cube = [[[0]*n for _ in range(n)] for _ in range(n)]
    for i in range(n):
        for j in range(n):
            cube[i][j][(i+j)%n] = 1
    proper = True; ir = ic = ik = 0; count = 0
    while count < steps or not proper:
        if proper:
            i = random.randint(0, n-1); j = random.randint(0, n-1)
            row_ij = cube[i][j]; k_cur = 0
            while row_ij[k_cur] != 1: k_cur += 1
            k = random.randint(0, n-2)
            if k >= k_cur: k += 1
            count += 1
        else:
            i, j, k = ir, ic, ik
            cands = [x for x in range(n) if cube[i][j][x] == 1]
            k_cur = cands[random.randint(0, len(cands)-1)]
        jc = [jj for jj in range(n) if jj != j and cube[i][jj][k] == 1]
        j_p = jc[random.randint(0, len(jc)-1)]
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
                if cube[i][j][k] == 1: L[i,j] = k+1; break
    return L

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
random.seed(2026)
print("=" * 70)
print("  PAPER CLAIMS VERIFICATION (CORRECTED DIVISIBILITY)")
print("=" * 70)
T0 = time.time()

# ── PART 1: Universal theorem h|det(A), h=n/2 for even n ────────────
print("\n>>> PART 1: Paper's Theorem — h|det(A) universal (h=n/2, even n)")
print("    (Also: n|det(A) for odd n)")

# n=4 exhaustive
all4 = all_ls_4()
fail_h4 = 0; fail_n4 = 0
for ls in all4:
    L = np.array(ls, dtype=int); d = int(bareiss(get_A(L)))
    if d != 0:
        if d % 2 != 0: fail_h4 += 1   # h=2
        if d % 4 != 0: fail_n4 += 1   # n=4
print(f"  n=4 exhaustive (576 LS): h=2|det fails={fail_h4}, n=4|det fails={fail_n4}")

# JM samples for various n
for n, samp, stp in [(6, 300, 2000), (8, 300, 3000),
                      (10, 100, 3000), (12, 150, 5000)]:
    h = n // 2
    t0 = time.time()
    fail_h = 0; fail_n = 0; det0 = 0
    for _ in range(samp):
        L = jm_sample(n, stp); d = int(bareiss(get_A(L)))
        if d == 0: det0 += 1; continue
        if d % h != 0: fail_h += 1
        if d % n != 0: fail_n += 1
    print(f"  n={n:2d} JM ({samp:3d} samples): h={h}|det fails={fail_h}, "
          f"n={n}|det fails={fail_n}, det=0: {det0}  [{time.time()-t0:.1f}s]")
    sys.stdout.flush()

# ── PART 2: F2 characterization for n ≡ 2 mod 4 ────────────────────
print("\n>>> PART 2: F2-rank characterization (n ≡ 2 mod 4)")
print("    n|det(A) <=> dim ker(B) >= 1")

for n, samp, stp in [(6, 500, 2000), (10, 200, 3000)]:
    t0 = time.time()
    ct = Counter()
    for _ in range(samp):
        L = jm_sample(n, stp); A = get_A(L); d = int(bareiss(A))
        kd = kernel_dim_F2(A)
        n_div = (d % n == 0) if d != 0 else True  # 0 is divisible by anything
        ct[(kd, n_div)] += 1

    print(f"\n  n={n} ({samp} JM samples, {stp} steps)  [{time.time()-t0:.1f}s]")
    print(f"  Cross-tab (kdim, {n}|det):")
    suf = nec = 0
    for k in sorted(ct):
        print(f"    kdim={k[0]}, {n}|det={k[1]} -> {ct[k]}")
        if k[0] >= 1 and not k[1]: suf += 1  # kdim>=1 but n∤det
        if k[0] < 1 and k[1]: nec += 1       # kdim=0 but n|det
    # Paper says: n|det <=> kdim>=1 (for n≡2 mod 4)
    # Equivalently: 2|det <=> rank < n-1
    suf_ct = sum(ct[k] for k in ct if k[0] >= 1 and not k[1])
    nec_ct = sum(ct[k] for k in ct if k[0] < 1 and k[1])
    print(f"  Criterion kdim>=1 <=> {n}|det: suff_fail={suf_ct} nec_fail={nec_ct}")
    sys.stdout.flush()

# ── PART 3: n ≡ 0 mod 4 — detailed cross-tab ────────────────────────
print("\n>>> PART 3: n ≡ 0 mod 4 — kdim vs v2(det) cross-tabulation")

# n=4 exhaustive
print("\n  n=4 EXHAUSTIVE (all 576):")
ct4_kd_v2 = Counter()
ct4_kd_ndiv = Counter()
for ls in all4:
    L = np.array(ls, dtype=int); A = get_A(L); d = int(bareiss(A))
    kd = kernel_dim_F2(A)
    ct4_kd_v2[(kd, v2(d))] += 1
    ct4_kd_ndiv[(kd, d != 0 and d % 4 == 0)] += 1  # n=4, n|det
print("  (kdim, v2) distribution:")
for k in sorted(ct4_kd_v2): print(f"    kdim={k[0]}, v2={k[1]} -> {ct4_kd_v2[k]}")
print("  (kdim, 4|det):")
for k in sorted(ct4_kd_ndiv): print(f"    kdim={k[0]}, 4|det={k[1]} -> {ct4_kd_ndiv[k]}")

# n=8 JM — KEY TEST
N8 = 500; S8 = 5000
print(f"\n  n=8 JM ({N8} samples, {S8} steps):")
t0 = time.time()
ct8_kd_v2 = Counter()
ct8_kd_8div = Counter()
ct8_kd_4div = Counter()
for trial in range(N8):
    L = jm_sample(8, S8); A = get_A(L); d = int(bareiss(A))
    kd = kernel_dim_F2(A)
    ct8_kd_v2[(kd, v2(d))] += 1
    if d != 0:
        ct8_kd_8div[(kd, d % 8 == 0)] += 1   # n|det
        ct8_kd_4div[(kd, d % 4 == 0)] += 1   # h|det
    if (trial+1) % 100 == 0: print(f"    {trial+1}/{N8}..."); sys.stdout.flush()
print(f"  [{time.time()-t0:.1f}s]")

print("  (kdim, 4|det) — paper's guaranteed h=4:")
for k in sorted(ct8_kd_4div): print(f"    kdim={k[0]}, 4|det={k[1]} -> {ct8_kd_4div[k]}")
h4_fail = sum(ct8_kd_4div[k] for k in ct8_kd_4div if not k[1])
print(f"  4|det(A) failures: {h4_fail}")

print("  (kdim, 8|det) — n²|det(E_std) question:")
for k in sorted(ct8_kd_8div): print(f"    kdim={k[0]}, 8|det={k[1]} -> {ct8_kd_8div[k]}")
suf8 = sum(ct8_kd_8div[k] for k in ct8_kd_8div if k[0] >= 2 and not k[1])
nec8 = sum(ct8_kd_8div[k] for k in ct8_kd_8div if k[0] < 2 and k[1])
print(f"  Criterion kdim>=2 <=> 8|det: suff_fail={suf8} nec_fail={nec8}")

print("  (kdim, v2) detailed:")
for k in sorted(ct8_kd_v2):
    if ct8_kd_v2[k] >= 2:
        print(f"    kdim={k[0]}, v2={k[1]} -> {ct8_kd_v2[k]}")

# n=12 JM
N12 = 200; S12 = 5000
print(f"\n  n=12 JM ({N12} samples, {S12} steps):")
t0 = time.time()
ct12_kd_v2 = Counter()
ct12_kd_12div = Counter()
ct12_kd_6div = Counter()
for trial in range(N12):
    L = jm_sample(12, S12); A = get_A(L); d = int(bareiss(A))
    kd = kernel_dim_F2(A)
    ct12_kd_v2[(kd, v2(d))] += 1
    if d != 0:
        ct12_kd_12div[(kd, d % 12 == 0)] += 1
        ct12_kd_6div[(kd, d % 6 == 0)] += 1
    if (trial+1) % 50 == 0: print(f"    {trial+1}/{N12}..."); sys.stdout.flush()
print(f"  [{time.time()-t0:.1f}s]")

print("  (kdim, 6|det) — paper's guaranteed h=6:")
for k in sorted(ct12_kd_6div): print(f"    kdim={k[0]}, 6|det={k[1]} -> {ct12_kd_6div[k]}")
h6_fail = sum(ct12_kd_6div[k] for k in ct12_kd_6div if not k[1])
print(f"  6|det(A) failures: {h6_fail}")

print("  (kdim, 12|det) — n²|det(E_std) question:")
for k in sorted(ct12_kd_12div): print(f"    kdim={k[0]}, 12|det={k[1]} -> {ct12_kd_12div[k]}")

# ── PART 4: Paper's sharpness claim at n=8 ──────────────────────────
# Paper says ~9.5% have n²/2 | det(E_std) but n² ∤ det(E_std)
# i.e., 32 | det(E_std) but 64 ∤ det(E_std)
# i.e., 4|det(A) but 8∤det(A)  (since det(E_std) = 8*det(A))
# wait: det(E_std) = n*det(A) = 8*det(A)
# 32|det(E_std) = 32|8*det(A) = 4|det(A)
# 64|det(E_std) = 64|8*det(A) = 8|det(A)
print("\n>>> PART 4: Sharpness at n=8")
# Use the n=8 data already collected
if ct8_kd_8div:
    total_nz = sum(ct8_kd_8div.values())
    sharp = sum(ct8_kd_8div[k] for k in ct8_kd_8div if k[1] == False)  # 8∤det but 4|det (if 4|det is universal)
    h_ok = sum(ct8_kd_4div[k] for k in ct8_kd_4div if k[1] == True)
    sharp_pct = 100 * sharp / total_nz if total_nz else 0
    print(f"  Total det!=0: {total_nz}")
    print(f"  4|det(A): {h_ok}/{total_nz} ({100*h_ok/total_nz:.1f}%)")
    print(f"  4|det but 8∤det: {sharp}/{total_nz} ({sharp_pct:.1f}%)")
    print(f"  (Paper claims ~9.5% for this)")

# ── PART 5: Corollary correctness check ─────────────────────────────
# Paper's Corollary says "when n is even, n²|det(E_std) iff rank<n-1"
# For n≡0 mod 4, rank<n-1 is ALWAYS TRUE (1∈ker), but n²|det is NOT always.
# This means the Corollary is INCORRECT for n≡0 mod 4.
print("\n>>> PART 5: Corollary correctness check")
print("  The Corollary states: 'When n is even, n²|det(E_std) iff rank_F2(B)<n-1'")
print("  For n≡0 mod 4: rank_F2(B)<n-1 ALWAYS (since 1∈ker)")
if ct8_kd_8div:
    n2_ok = sum(ct8_kd_8div[k] for k in ct8_kd_8div if k[1])
    n2_fail = sum(ct8_kd_8div[k] for k in ct8_kd_8div if not k[1])
    print(f"  n=8: rank<7 always TRUE, but 8|det(A) (=n²|det(E_std)):")
    print(f"    TRUE: {n2_ok}  FALSE: {n2_fail}")
    if n2_fail > 0:
        print(f"  ==> Corollary is INCORRECT for n≡0 mod 4! "
              f"(rank<n-1 is necessary but not sufficient)")
    else:
        print(f"  ==> Corollary holds at n=8")

print(f"\n{'='*70}")
print(f"  TOTAL: {time.time()-T0:.1f}s")
print(f"{'='*70}")
