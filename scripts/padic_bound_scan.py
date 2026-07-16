"""
p-adic Sharp Bound Verification (the paper's "p-adic bound scan").

Research-phase script, added to the public repository 2026-07-09 as the
generator of `results/certified/padic_scan_output.txt` (the "p-adic
bound scan" item of the paper's Code-and-data-availability paragraph).
The campaign covers n in {6,8,9,10,12,15,16,18,20} (beyond the promised n <= 12),
with fixed per-order seeds (42..50); all primes p | n are scanned,
which includes p in {2, 3, 5} on the promised range.

Verifies the GENERALIZED sharp bound for ALL primes p | n
(paper `thm:sharp-v2` for p = 2 and `sec:p-adic` for odd p):

    v_p(det A)  >=  v_p(n/2) + max(0, k_p - 1)

where k_p = dim ker_{F_p}(A mod p).

Key insight: the proof of the sharp 2-adic bound uses ONLY:
  1.  A·1 = h·c  where h = n/2                   (row-sum identity)
  2.  d_{m} = h·(Q·c)_i  via adjugate identity    (=> h | d_m)
  3.  k_p invariant factors are p-divisible         (standard SNF)
  4.  Pigeonhole on the divisibility chain

NONE of these steps are 2-specific. The proof works for ALL primes.

Also checks the WEAKER bound v_p(n/p) + max(0, k_p - 1)
and shows the stronger bound v_p(n/2) is the correct one.
"""

import random
import sys
from fractions import Fraction
from collections import defaultdict
from math import gcd

# ── Jacobson-Matthews sampler ──────────────────────────────────────

def jm_sample(n, rng, steps=None):
    if steps is None:
        steps = 5 * n * n
    cube = [[[0]*n for _ in range(n)] for _ in range(n)]
    for i in range(n):
        for j in range(n):
            cube[i][j][(i+j) % n] = 1
    proper = True; ir = ic = ik = 0; count = 0
    while count < steps or not proper:
        if proper:
            i = rng.randint(0, n-1); j = rng.randint(0, n-1)
            k_cur = next(k for k in range(n) if cube[i][j][k] == 1)
            k = rng.randint(0, n-2)
            if k >= k_cur: k += 1
            count += 1
        else:
            i, j, k = ir, ic, ik
            cands = [x for x in range(n) if cube[i][j][x] == 1]
            k_cur = cands[rng.randint(0, len(cands)-1)]
        jc = [jj for jj in range(n) if jj != j and cube[i][jj][k] == 1]
        j_p = jc[rng.randint(0, len(jc)-1)]
        ic_ = [ii for ii in range(n) if ii != i and cube[ii][j][k] == 1]
        i_p = ic_[rng.randint(0, len(ic_)-1)]
        cube[i][j][k] += 1;     cube[i][j][k_cur] -= 1
        cube[i][j_p][k] -= 1;   cube[i][j_p][k_cur] += 1
        cube[i_p][j][k] -= 1;   cube[i_p][j][k_cur] += 1
        cube[i_p][j_p][k] += 1; cube[i_p][j_p][k_cur] -= 1
        if cube[i_p][j_p][k_cur] == -1:
            proper = False; ir, ic, ik = i_p, j_p, k_cur
        else:
            proper = True
    L = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(n):
            for k2 in range(n):
                if cube[i][j][k2] == 1:
                    L[i][j] = k2 + 1; break
    return L

# ── Exact integer determinant (Fraction-based Bareiss) ──────────────

def exact_det(M):
    n = len(M)
    A = [[Fraction(M[i][j]) for j in range(n)] for i in range(n)]
    sign = 1
    for col in range(n):
        pivot = None
        for row in range(col, n):
            if A[row][col] != 0:
                pivot = row; break
        if pivot is None:
            return 0
        if pivot != col:
            A[col], A[pivot] = A[pivot], A[col]
            sign *= -1
        for row in range(col + 1, n):
            if A[row][col] != 0:
                f = A[row][col] / A[col][col]
                for j in range(col, n):
                    A[row][j] -= f * A[col][j]
    result = Fraction(sign)
    for i in range(n):
        result *= A[i][i]
    return int(result)

# ── p-adic valuation ─────────────────────────────────────────────────

def vp(m, p):
    """p-adic valuation of integer m."""
    if m == 0:
        return float('inf')
    m = abs(m)
    v = 0
    while m % p == 0:
        m //= p
        v += 1
    return v

# ── Prime factorization (small n) ───────────────────────────────────

def prime_factors(n):
    """Return set of prime factors of n."""
    factors = set()
    d = 2
    while d * d <= n:
        while n % d == 0:
            factors.add(d)
            n //= d
        d += 1
    if n > 1:
        factors.add(n)
    return factors

# ── Kernel dimension over F_p ────────────────────────────────────────

def fp_kernel_dim(B, p):
    """Compute dim ker(B mod p) via Gaussian elimination over F_p."""
    m = len(B)
    n_cols = len(B[0])
    M = [[B[i][j] % p for j in range(n_cols)] for i in range(m)]
    pivot_count = 0
    row = 0
    for col in range(n_cols):
        found = None
        for r in range(row, m):
            if M[r][col] % p != 0:
                found = r
                break
        if found is None:
            continue
        M[row], M[found] = M[found], M[row]
        inv = pow(M[row][col], p - 2, p)  # Fermat inverse
        for j in range(n_cols):
            M[row][j] = (M[row][j] * inv) % p
        for r in range(m):
            if r != row and M[r][col] % p != 0:
                f = M[r][col]
                for j in range(n_cols):
                    M[r][j] = (M[r][j] - f * M[row][j]) % p
        pivot_count += 1
        row += 1
    rank = pivot_count
    return n_cols - rank

# ── Build A matrix ───────────────────────────────────────────────────

def build_A(L):
    n = len(L)
    m = n - 1
    A = [[L[i][j] - L[i][m] for j in range(m)] for i in range(m)]
    return A

# ── Main verification ────────────────────────────────────────────────

def verify_bound(n, num_samples, seed):
    """
    For order n, sample num_samples LS, check p-adic bound at all p | n.
    Returns dict of results.
    """
    rng = random.Random(seed)
    h = n // 2 if n % 2 == 0 else n  # h = n/2 for even, n for odd
    primes = sorted(prime_factors(n))
    
    print(f"\n{'='*70}")
    print(f"  n = {n},  primes = {primes},  h = n/2 = {h if n % 2 == 0 else 'N/A (odd)'}")
    print(f"  Samples = {num_samples}")
    print(f"{'='*70}")
    
    # Track results per prime
    results = {}
    for p in primes:
        results[p] = {
            'tests': 0,
            'violation_strong': 0,  # v_p(n/2) + max(0, k_p-1)
            'violation_weak': 0,    # v_p(n/p) + max(0, k_p-1) [paper conjecture]
            'tight_count': 0,       # cases where strong bound is exactly met
            'kp_dist': defaultdict(int),
            'vp_min_by_kp': {},     # min v_p observed per k_p
            'gap_sum': 0,           # sum of (actual - bound) for excess analysis
        }
    
    for trial in range(num_samples):
        L = jm_sample(n, rng)
        A = build_A(L)
        det_val = exact_det(A)
        
        for p in primes:
            kp = fp_kernel_dim(A, p)
            vp_det = vp(det_val, p)
            
            # Strong bound: v_p(n/2) + max(0, k_p - 1)
            if n % 2 == 0:
                vp_h = vp(n // 2, p)
            else:
                # For odd n, row sum = n * ((n+1-2L)/2), so n | row sum
                # The analogue of h is more subtle; use n for odd n
                vp_h = vp(n, p)
            
            strong_bound = vp_h + max(0, kp - 1)
            
            # Weak bound (paper conjecture): v_p(n/p) + max(0, k_p - 1)
            weak_bound = vp(n // p, p) + max(0, kp - 1)
            # Note: v_p(n/p) = v_p(n) - 1
            
            r = results[p]
            r['tests'] += 1
            r['kp_dist'][kp] += 1
            
            if vp_det < strong_bound:
                r['violation_strong'] += 1
                print(f"  *** STRONG VIOLATION at p={p}: n={n}, "
                      f"v_{p}(det)={vp_det} < {strong_bound} "
                      f"(v_{p}(h)={vp_h}, k_{p}={kp})")
            
            if vp_det < weak_bound:
                r['violation_weak'] += 1
            
            if vp_det == strong_bound:
                r['tight_count'] += 1
            
            gap = vp_det - strong_bound
            r['gap_sum'] += gap
            
            key = kp
            if key not in r['vp_min_by_kp'] or vp_det < r['vp_min_by_kp'][key]:
                r['vp_min_by_kp'][key] = vp_det
    
    # Print results
    for p in primes:
        r = results[p]
        vp_h_val = vp(n // 2, p) if n % 2 == 0 else vp(n, p)
        print(f"\n  Prime p = {p}:")
        print(f"    v_{p}(h) = v_{p}(n/2) = {vp_h_val}")
        print(f"    v_{p}(n/p) = {vp(n // p, p)}  [paper's weaker conjecture]")
        print(f"    Tests: {r['tests']}")
        print(f"    Strong bound violations (v_p(n/2)): {r['violation_strong']}")
        print(f"    Weak bound violations (v_p(n/p)):   {r['violation_weak']}")
        print(f"    Tight cases (strong bound exact):    {r['tight_count']}")
        print(f"    k_p distribution: {dict(sorted(r['kp_dist'].items()))}")
        print(f"    Min v_p by k_p:   {dict(sorted(r['vp_min_by_kp'].items()))}")
        
        # Show bound vs min observed
        print(f"    {'k_p':>4} | {'bound':>6} | {'min_vp':>6} | {'tight?':>6}")
        print(f"    {'----':>4}-+-{'------':>6}-+-{'------':>6}-+-{'------':>6}")
        for kp_val in sorted(r['vp_min_by_kp'].keys()):
            bound = vp_h_val + max(0, kp_val - 1)
            min_obs = r['vp_min_by_kp'][kp_val]
            tight = "YES" if min_obs == bound else f"+{min_obs - bound}"
            print(f"    {kp_val:>4} | {bound:>6} | {min_obs:>6} | {tight:>6}")
    
    return results

# ── Also verify the key lemma: 1 ∈ ker(A mod p) for odd p | n ────

def verify_kernel_membership(n, num_samples, seed):
    """Verify that 1 ∈ ker(A mod p) for all odd primes p | n."""
    rng = random.Random(seed)
    primes = [p for p in prime_factors(n) if p > 2]
    if not primes:
        return
    
    print(f"\n  Verifying 1 ∈ ker(A mod p) for odd primes {primes}:")
    for trial in range(num_samples):
        L = jm_sample(n, rng)
        A = build_A(L)
        m = len(A)
        for p in primes:
            # Check A·1 mod p = 0
            for i in range(m):
                row_sum = sum(A[i][j] for j in range(m)) % p
                if row_sum != 0:
                    print(f"  *** FAIL: A·1 mod {p} ≠ 0 at n={n}, row {i}")
                    return
    print(f"    All {num_samples} samples: 1 ∈ ker(A mod p) confirmed for all odd p. ✓")

# ── Also verify h | d_{m} directly via adj identity  ────────────────

def verify_h_divides_dm(n, num_samples, seed):
    """
    Verify h | d_{n-1} (largest invariant factor).
    We check: det(A) / gcd_of_all_(m-1)x(m-1)_minors is divisible by h.
    For efficiency, just check h | det(A) (weaker but necessary).
    """
    rng = random.Random(seed)
    h = n // 2
    violations = 0
    for trial in range(num_samples):
        L = jm_sample(n, rng)
        A = build_A(L)
        det_val = exact_det(A)
        if det_val != 0 and det_val % h != 0:
            violations += 1
            print(f"  *** h={h} does not divide det(A)={det_val}")
    print(f"\n  h = n/2 = {h} divides det(A): "
          f"{num_samples - violations}/{num_samples} passed "
          f"({'✓ ALL PASS' if violations == 0 else f'*** {violations} VIOLATIONS'})")


def run_campaign():
    """Full verification campaign across multiple orders."""
    
    configs = [
        # (n, samples, seed)
        (6,   300, 42),
        (8,   300, 43),
        (9,   200, 44),    # odd n, p=3
        (10,  200, 45),
        (12,  200, 46),
        (15,  150, 47),    # odd n, p=3,5
        (16,  100, 48),
        (18,  100, 49),
        (20,  100, 50),
    ]
    
    total_tests = 0
    total_violations = 0
    all_results = {}
    
    for n, samples, seed in configs:
        verify_kernel_membership(n, min(samples, 50), seed)
        if n % 2 == 0:
            verify_h_divides_dm(n, min(samples, 50), seed + 1000)
        results = verify_bound(n, samples, seed)
        all_results[n] = results
        for p, r in results.items():
            total_tests += r['tests']
            total_violations += r['violation_strong']
    
    # Grand summary
    print(f"\n{'='*70}")
    print(f"  GRAND SUMMARY")
    print(f"{'='*70}")
    print(f"  Total (n, p, sample) tests: {total_tests}")
    print(f"  Strong bound violations:     {total_violations}")
    print(f"  Verdict: {'✓ NO VIOLATIONS IN THE SAMPLED CAMPAIGN' if total_violations == 0 else '*** VIOLATIONS FOUND ***'}")
    
    # Compact table
    print(f"\n  {'n':>3} | {'p':>2} | {'v_p(n/2)':>8} | {'tests':>5} | {'viol':>4} | {'tight':>5} | k_p dist")
    print(f"  {'---':>3}-+-{'--':>2}-+-{'--------':>8}-+-{'-----':>5}-+-{'----':>4}-+-{'-----':>5}-+----------")
    for n, results in sorted(all_results.items()):
        for p, r in sorted(results.items()):
            vph = vp(n // 2, p) if n % 2 == 0 else vp(n, p)
            kd = dict(sorted(r['kp_dist'].items()))
            kd_str = ", ".join(f"{k}:{v}" for k, v in kd.items())
            print(f"  {n:>3} | {p:>2} | {vph:>8} | {r['tests']:>5} | "
                  f"{r['violation_strong']:>4} | {r['tight_count']:>5} | {kd_str}")
    
    return total_violations == 0


if __name__ == "__main__":
    success = run_campaign()
    sys.exit(0 if success else 1)
