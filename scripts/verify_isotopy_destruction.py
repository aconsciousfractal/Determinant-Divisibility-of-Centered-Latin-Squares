"""verify_isotopy_destruction.py — the standalone script promised in the
paper's "Code and data availability" paragraph (Sec. 8).

Added to the public repository 2026-07-09 (the script was promised at
publication but missing until then — see CLAIM_LEDGER.md GAP-1).

Reproduces the isotopy non-invariance experiment of the Discussion
section, using the explicit n=10 counterexample of Example ex:n10:

  1. Verify the witness invariants stated in the paper's table:
     valid Latin square; doubly balanced parity (5 odd + 5 even per
     row and column); det(A) = 15 427 045 = 5 x 3 085 409 (odd);
     det(E_std) = 10 * det(A) = 154 270 450; 50 | det(E_std) but
     100 does NOT divide det(E_std); rank_F2(A mod 2) = 9 (full).
  2. Reduce L to standard form (first row = (1,...,10)) by the symbol
     relabeling gamma with gamma(L[0][j]) = j+1, and verify that the
     reduced square has EVEN det(A) — i.e. the counterexample property
     is destroyed by an isotopy, as claimed in the Discussion.
  3. Write the witness JSON certificate
     results/certified/ex_n10_witness.json (the "ex:n10 witness JSON
     certificate" of the availability paragraph).

All computations are exact (integer Bareiss determinant, F2 Gaussian
elimination). Deterministic: no randomness.
"""
import sys as _sys
try:
    _sys.stdout.reconfigure(encoding="utf-8")
    _sys.stderr.reconfigure(encoding="utf-8")
except Exception:
    pass

import json
import os

# Example ex:n10 of the paper (Sec. 6, "Explicit counterexample at n=10").
L10 = [
    [2, 5, 6, 4, 3, 8, 7, 9, 10, 1],
    [10, 6, 2, 7, 5, 1, 9, 8, 4, 3],
    [3, 4, 1, 10, 8, 5, 2, 7, 9, 6],
    [7, 1, 10, 8, 6, 9, 3, 4, 2, 5],
    [1, 7, 3, 6, 4, 10, 8, 2, 5, 9],
    [9, 2, 5, 3, 10, 6, 4, 1, 8, 7],
    [6, 9, 8, 1, 2, 7, 5, 10, 3, 4],
    [5, 8, 7, 2, 9, 4, 1, 3, 6, 10],
    [4, 3, 9, 5, 1, 2, 10, 6, 7, 8],
    [8, 10, 4, 9, 7, 3, 6, 5, 1, 2],
]


def is_latin_square(L):
    n = len(L)
    expected = set(range(1, n + 1))
    for i in range(n):
        if set(L[i]) != expected:
            return False
    for j in range(n):
        if set(L[i][j] for i in range(n)) != expected:
            return False
    return True


def compute_A(L):
    """The (n-1)x(n-1) difference matrix A (last column as reference)."""
    n = len(L)
    return [[L[i][j] - L[i][n - 1] for j in range(n - 1)]
            for i in range(n - 1)]


def bareiss_det(M):
    """Exact integer determinant via Bareiss elimination."""
    n = len(M)
    A = [row[:] for row in M]
    sign = 1
    prev = 1
    for i in range(n):
        if A[i][i] == 0:
            found = False
            for k in range(i + 1, n):
                if A[k][i] != 0:
                    A[i], A[k] = A[k], A[i]
                    sign *= -1
                    found = True
                    break
            if not found:
                return 0
        for j in range(i + 1, n):
            for k in range(i + 1, n):
                A[j][k] = (A[j][k] * A[i][i] - A[j][i] * A[i][k]) // prev
            A[j][i] = 0
        prev = A[i][i]
    return sign * A[n - 1][n - 1]


def rank_F2(M):
    """Rank of an integer matrix mod 2 via F2 Gaussian elimination."""
    B = [[x % 2 for x in row] for row in M]
    rows, cols = len(B), len(B[0])
    pr = 0
    for col in range(cols):
        pivot = None
        for row in range(pr, rows):
            if B[row][col] == 1:
                pivot = row
                break
        if pivot is None:
            continue
        B[pr], B[pivot] = B[pivot], B[pr]
        for row in range(rows):
            if row != pr and B[row][col] == 1:
                B[row] = [(a + b) % 2 for a, b in zip(B[row], B[pr])]
        pr += 1
    return pr


def standard_form(L):
    """Symbol relabeling gamma with gamma(L[0][j]) = j+1 (first row -> 1..n)."""
    n = len(L)
    gamma = {L[0][j]: j + 1 for j in range(n)}
    return [[gamma[x] for x in row] for row in L]


def main():
    n = 10
    checks = []

    def check(name, ok, detail=""):
        status = "PASS" if ok else "FAIL"
        print(f"  [{status}] {name}" + (f" — {detail}" if detail else ""))
        checks.append((name, bool(ok)))
        return ok

    print("=== ex:n10 witness invariants ===")
    check("valid Latin square", is_latin_square(L10))

    P = [[x % 2 for x in row] for row in L10]
    balanced = (all(sum(row) == 5 for row in P)
                and all(sum(P[i][j] for i in range(n)) == 5 for j in range(n)))
    check("doubly balanced parity (5 odd per row and column)", balanced)

    A = compute_A(L10)
    det_A = bareiss_det(A)
    check("det(A) = 15427045", det_A == 15427045, f"det(A) = {det_A}")
    check("det(A) odd (v_2 = 0)", det_A % 2 == 1)
    check("det(A) = 5 * 3085409", det_A == 5 * 3085409)

    det_E = n * det_A
    check("det(E_std) = 154270450", det_E == 154270450, f"det(E_std) = {det_E}")
    check("n^2/2 = 50 divides det(E_std)", det_E % 50 == 0)
    check("n^2 = 100 does NOT divide det(E_std)", det_E % 100 != 0)

    r = rank_F2(A)
    check("rank_F2(A mod 2) = 9 (full)", r == 9, f"rank = {r}")

    print("\n=== isotopy destruction (Discussion) ===")
    L_std = standard_form(L10)
    check("standard form is a Latin square", is_latin_square(L_std))
    check("standard form first row = (1,...,10)",
          L_std[0] == list(range(1, n + 1)))
    det_A_std = bareiss_det(compute_A(L_std))
    check("standard form has EVEN det(A) — counterexample destroyed",
          det_A_std % 2 == 0, f"det(A_std) = {det_A_std}")

    all_ok = all(ok for _, ok in checks)

    out_dir = os.path.join(os.path.dirname(__file__), "..", "results",
                           "certified")
    os.makedirs(out_dir, exist_ok=True)
    cert = {
        "certificate": "ex:n10 witness (paper Example ex:n10)",
        "generator": "scripts/verify_isotopy_destruction.py",
        "L": L10,
        "invariants": {
            "det_A": det_A,
            "det_E_std": det_E,
            "v2_det_A": 0,
            "rank_F2_A_mod_2": r,
            "doubly_balanced": balanced,
            "divisible_by_50": det_E % 50 == 0,
            "divisible_by_100": det_E % 100 == 0,
        },
        "isotopy_destruction": {
            "L_standard_form": L_std,
            "det_A_standard_form": det_A_std,
            "det_A_standard_form_even": det_A_std % 2 == 0,
        },
        "all_checks_pass": all_ok,
    }
    out_path = os.path.join(out_dir, "ex_n10_witness.json")
    with open(out_path, "w") as f:
        json.dump(cert, f, indent=2)
    print(f"\nWitness certificate written: {out_path}")

    print("ALL CHECKS PASS" if all_ok else "SOME CHECKS FAILED")
    return 0 if all_ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
