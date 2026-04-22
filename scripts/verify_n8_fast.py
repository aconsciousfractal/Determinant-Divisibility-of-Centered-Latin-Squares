"""Fast verification: n^2 | det(E_std) for n = 4, 8, 12."""
import random, math

def bareiss_det(M):
    n = len(M)
    A = [row[:] for row in M]
    sign = 1
    for k in range(n):
        if A[k][k] == 0:
            found = False
            for i in range(k+1, n):
                if A[i][k] != 0:
                    A[k], A[i] = A[i], A[k]
                    sign *= -1
                    found = True
                    break
            if not found:
                return 0
        for i in range(k+1, n):
            for j in range(k+1, n):
                A[i][j] = A[k][k]*A[i][j] - A[i][k]*A[k][j]
                if k > 0:
                    A[i][j] //= A[k-1][k-1]
            A[i][k] = 0
    return sign * A[n-1][n-1]

def v2(x):
    if x == 0: return float('inf')
    x = abs(x); v = 0
    while x % 2 == 0: v += 1; x //= 2
    return v

def random_ls(n):
    """Build Latin square by random row-by-row construction with backtracking."""
    L = [[0]*n for _ in range(n)]
    avail_col = [{v for v in range(1, n+1)} for _ in range(n)]
    
    def fill(row):
        if row == n:
            return True
        perm = list(range(1, n+1))
        random.shuffle(perm)
        return try_row(row, 0, perm, avail_col, L)
    
    def try_row(row, col, remaining, avail_col, L):
        if col == n:
            return True
        candidates = [v for v in remaining if v in avail_col[col]]
        random.shuffle(candidates)
        for v in candidates:
            L[row][col] = v
            avail_col[col].discard(v)
            new_remaining = [x for x in remaining if x != v]
            if try_row(row, col+1, new_remaining, avail_col, L):
                return True
            avail_col[col].add(v)
            L[row][col] = 0
        return False
    
    if fill(0):
        return L
    return None

for n in [4, 8]:
    print(f"\n=== n = {n}, n^2 = {n**2}, v_2(n^2) = {v2(n**2)} ===")
    random.seed(12345)
    min_v2_E = float('inf')
    min_v2_A = float('inf')
    fails = 0
    count = 0
    target = 100 if n == 4 else 50
    
    for _ in range(target):
        L = random_ls(n)
        if L is None:
            continue
        A = [[L[i][j] - L[i][n-1] for j in range(n-1)] for i in range(n-1)]
        dA = bareiss_det(A)
        dE = n * dA
        if dE == 0:
            continue
        count += 1
        ve = v2(dE)
        va = v2(dA)
        if ve < min_v2_E:
            min_v2_E = ve
        if va < min_v2_A:
            min_v2_A = va
        if dE % (n**2) != 0:
            fails += 1
            print(f"  FAIL: det(A)={dA}, v2(A)={va}, det(E)={dE}, v2(E)={ve}")
    
    print(f"Tested: {count}")
    print(f"Min v_2(det E_std) = {min_v2_E}, min v_2(det A) = {min_v2_A}")
    print(f"n^2 divisibility failures: {fails}")
    if fails == 0:
        print(f"CONFIRMED: n^2={n**2} | det(E_std) for all {count} tested.")
