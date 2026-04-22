"""
Verify reviewer claim: det(E_std)=0 exists at n=4.
Check all 4 reduced Latin squares of order 4.
"""
from itertools import permutations

def is_latin_square(L, n):
    for i in range(n):
        if sorted(L[i]) != list(range(1, n+1)):
            return False
        col = [L[r][i] for r in range(n)]
        if sorted(col) != list(range(1, n+1)):
            return False
    return True

def generate_reduced_ls(n):
    """Generate all reduced LS of order n (first row=1..n, first col=1..n)."""
    results = []
    first_row = list(range(1, n+1))
    
    def fill(row, partial):
        if row == n:
            results.append([r[:] for r in partial])
            return
        # Row `row` must be a permutation of 1..n
        # with first_col constraint: partial[row][0] = row+1
        # and column constraints
        remaining = set(range(1, n+1))
        # Column 0 is fixed: partial[row][0] = row+1 already placed
        used_in_cols = [set() for _ in range(n)]
        for r in range(row):
            for c in range(n):
                used_in_cols[c].add(partial[r][c])
        
        # Generate valid rows
        def gen_row(col, current_row):
            if col == n:
                partial.append(current_row[:])
                fill(row + 1, partial)
                partial.pop()
                return
            for val in range(1, n+1):
                if val not in current_row and val not in used_in_cols[col]:
                    current_row.append(val)
                    used_in_cols[col].add(val)
                    gen_row(col + 1, current_row)
                    used_in_cols[col].discard(val)
                    current_row.pop()
        
        # First column is fixed
        start = [row + 1]
        used_in_cols[0].add(row + 1)
        gen_row(1, start)
        used_in_cols[0].discard(row + 1)
    
    fill(1, [first_row])
    return results

def bareiss(M):
    m = len(M); A = [row[:] for row in M]; sign = 1; prev = 1
    for i in range(m):
        if A[i][i] == 0:
            found = False
            for k in range(i+1, m):
                if A[k][i] != 0:
                    A[i], A[k] = A[k], A[i]; sign *= -1; found = True; break
            if not found:
                return 0
        for j in range(i+1, m):
            for k in range(i+1, m):
                A[j][k] = (A[j][k]*A[i][i] - A[j][i]*A[i][k]) // prev
            A[j][i] = 0
        prev = A[i][i]
    return sign * A[m-1][m-1]

print("All reduced Latin squares of order 4:")
print("=" * 60)
ls4 = generate_reduced_ls(4)
print(f"Count: {len(ls4)}")

for idx, L in enumerate(ls4):
    n = 4
    A = [[L[i][j] - L[i][n-1] for j in range(n-1)] for i in range(n-1)]
    det_A = bareiss(A)
    det_Estd = n * det_A
    print(f"\n  LS #{idx+1}: {L}")
    print(f"  A = {A}")
    print(f"  det(A) = {det_A}, det(E_std) = {det_Estd}")
    if det_A == 0:
        print(f"  *** SINGULAR! det=0 ***")

print("\n" + "=" * 60)
print("All reduced Latin squares of order 3:")
ls3 = generate_reduced_ls(3)
print(f"Count: {len(ls3)}")
for idx, L in enumerate(ls3):
    n = 3
    A = [[L[i][j] - L[i][n-1] for j in range(n-1)] for i in range(n-1)]
    det_A = bareiss(A)
    print(f"  LS #{idx+1}: det(A)={det_A}, det(E_std)={n*det_A}")

print("\n" + "=" * 60)
print("Check n=5: any singular LS?")
ls5 = generate_reduced_ls(5)
print(f"Count: {len(ls5)}")
singular_count = 0
for L in ls5:
    n = 5
    A = [[L[i][j] - L[i][n-1] for j in range(n-1)] for i in range(n-1)]
    det_A = bareiss(A)
    if det_A == 0:
        singular_count += 1
        print(f"  SINGULAR: {L}")
print(f"  Singular count at n=5: {singular_count}/{len(ls5)}")
