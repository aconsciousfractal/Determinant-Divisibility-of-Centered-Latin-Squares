"""Fast verification: n^2 | det(E_std) for n = 4, 8.
Uses isotopy operations on cyclic squares to generate diverse samples."""
import random

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

def make_cyclic(n):
    return [[(i+j) % n + 1 for j in range(n)] for i in range(n)]

def apply_isotopy(L, n, rng):
    """Apply random isotopy: row perm, col perm, symbol relabeling."""
    L2 = [row[:] for row in L]
    # Row permutation
    rows = list(range(n))
    rng.shuffle(rows)
    L2 = [L2[rows[i]] for i in range(n)]
    # Column permutation
    cols = list(range(n))
    rng.shuffle(cols)
    L2 = [[L2[i][cols[j]] for j in range(n)] for i in range(n)]
    # Symbol relabeling
    syms = list(range(1, n+1))
    rng.shuffle(syms)
    mapping = {old: syms[old-1] for old in range(1, n+1)}
    L2 = [[mapping[L2[i][j]] for j in range(n)] for i in range(n)]
    return L2

def is_valid(L, n):
    for i in range(n):
        if sorted(L[i]) != list(range(1, n+1)):
            return False
        if sorted(L[j][i] for j in range(n)) != list(range(1, n+1)):
            return False
    return True

# Also test non-cyclic starting squares
def make_back_circulant(n):
    """L[i][j] = (i*2 + j) mod n + 1 — different LS type."""
    return [[(i*2 + j) % n + 1 for j in range(n)] for i in range(n)]

# For n=4: enumerate ALL 576 reduced Latin squares
def enumerate_n4():
    """Generate all reduced LS of order 4 by brute force."""
    # First row fixed: [1,2,3,4]
    # First col fixed: [1,2,3,4]
    # Fill remaining 3x3 subgrid
    results = []
    for r1 in [[2,3,4],[2,4,3],[3,2,4],[3,4,2],[4,2,3],[4,3,2]]:
        row1 = [2] + list(r1)  # first col = 2
        # Wait, first col must be [1,2,3,4]. Row 0 = [1,2,3,4]. Row 1 starts with 2.
        # Actually for reduced: first row = first col = (1,2,...,n).
        # Row 0: [1,2,3,4], Row 1: [2,?,?,?], Row 2: [3,?,?,?], Row 3: [4,?,?,?]
        pass
    
    # Simpler: just enumerate
    from itertools import permutations
    n = 4
    first_row = list(range(1, n+1))
    squares = []
    # Row 1 starts with 2
    for r1 in permutations([1,3,4]):
        row1 = [2] + list(r1)
        # Check column constraints with row 0
        col_used = [{first_row[j], row1[j]} for j in range(n)]
        for r2 in permutations([1,2,4]):
            row2 = [3] + list(r2)
            ok = True
            for j in range(n):
                if row2[j] in col_used[j]:
                    ok = False
                    break
            if not ok:
                continue
            # Row 3 is determined
            row3 = [4] + [0]*3
            for j in range(1, n):
                remaining = {1,2,3,4} - {first_row[j], row1[j], row2[j]}
                if len(remaining) != 1:
                    ok = False
                    break
                row3[j] = remaining.pop()
            if not ok:
                continue
            if sorted(row3) == [1,2,3,4]:
                squares.append([first_row, row1, row2, row3])
    return squares

rng = random.Random(42)

# n=4 exhaustive
print("=== n = 4 (exhaustive reduced LS) ===")
n = 4
squares_4 = enumerate_n4()
print(f"Found {len(squares_4)} reduced Latin squares of order 4")
min_v2_E = float('inf')
min_v2_A = float('inf')
fails = 0
for L in squares_4:
    A = [[L[i][j] - L[i][n-1] for j in range(n-1)] for i in range(n-1)]
    dA = bareiss_det(A)
    dE = n * dA
    if dE == 0:
        print(f"  det(E_std) = 0!")
        continue
    ve = v2(dE)
    va = v2(dA)
    if ve < min_v2_E: min_v2_E = ve
    if va < min_v2_A: min_v2_A = va
    if dE % (n**2) != 0:
        fails += 1
        print(f"  FAIL: det(A)={dA}, v2(A)={va}, det(E)={dE}, v2(E)={ve}")

print(f"Min v_2(det E_std) = {min_v2_E}, v_2(n^2) = {v2(n**2)}")
print(f"Min v_2(det A) = {min_v2_A}, v_2(n) = {v2(n)}")
print(f"Failures: {fails}")

# n=8 via isotopy from cyclic
print(f"\n=== n = 8 (isotopy of cyclic, 200 samples) ===")
n = 8
base = make_cyclic(n)
min_v2_E = float('inf')
min_v2_A = float('inf')
fails = 0
count = 0
for trial in range(200):
    L = apply_isotopy(base, n, rng)
    assert is_valid(L, n), f"Invalid at trial {trial}"
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
        if trial < 10 or ve <= v2(n**2):
            print(f"  Trial {trial}: min v_2(E_std) = {ve}, det(A)={dA}, v_2(A)={va}")
    if va < min_v2_A: min_v2_A = va
    if dE % (n**2) != 0:
        fails += 1
        print(f"  FAIL trial {trial}: det(A)={dA}, v2(A)={va}, det(E)={dE}, v2(E)={ve}")

print(f"\nTested: {count}")
print(f"Min v_2(det E_std) = {min_v2_E}, v_2(n^2) = {v2(n**2)}")
print(f"Min v_2(det A) = {min_v2_A}, v_2(n) = {v2(n)}")
print(f"Failures: {fails}")

# Also test from back-circulant base (different isotopy class potentially)
print(f"\n=== n = 8 (isotopy of back-circulant, 200 samples) ===")
base2 = make_back_circulant(n)
if is_valid(base2, n):
    min_v2_E = float('inf')
    min_v2_A = float('inf')
    fails = 0
    count = 0
    for trial in range(200):
        L = apply_isotopy(base2, n, rng)
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
        if va < min_v2_A: min_v2_A = va
        if dE % (n**2) != 0:
            fails += 1
            print(f"  FAIL: det(A)={dA}, v2(A)={va}, det(E)={dE}, v2(E)={ve}")
    print(f"Tested: {count}")
    print(f"Min v_2(det E_std) = {min_v2_E}, v_2(n^2) = {v2(n**2)}")
    print(f"Failures: {fails}")
else:
    print("Back-circulant is not a valid LS for n=8")
