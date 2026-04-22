"""Core linear-algebra primitives used throughout the paper.

Conventions (matching the manuscript):

* ``L`` is an ``n x n`` Latin square with entries in ``{1, ..., n}``.
* ``E = L - (n+1)/2 * J`` is the centered matrix (entries are half-integers
  when ``n`` is even, integers when ``n`` is odd).
* ``P`` is the ``n x (n-1)`` basis matrix of the standard subspace
  ``V_std = ker(1^T)``, with columns ``b_i = e_i - e_{n-1}``.
* ``G = P^T P = I_{n-1} + J_{n-1}`` has ``det(G) = n``.
* ``A`` is the ``(n-1) x (n-1)`` integer matrix ``A_{ij} = L_{ij} - L_{i,n-1}``.
* ``E_std = P^T E P`` satisfies ``E_std = G * A`` (Lemma 2.4), hence
  ``det(E_std) = n * det(A)``.

All routines work over exact integer arithmetic using Python ``int``.
"""

from __future__ import annotations

from fractions import Fraction
from typing import Sequence

import numpy as np


Matrix = np.ndarray


def is_latin_square(L: Matrix) -> bool:
    """Return True iff ``L`` is a valid Latin square on ``{1, ..., n}``."""
    L = np.asarray(L, dtype=int)
    n = L.shape[0]
    if L.shape != (n, n):
        return False
    symbols = set(range(1, n + 1))
    for i in range(n):
        if set(L[i, :].tolist()) != symbols:
            return False
        if set(L[:, i].tolist()) != symbols:
            return False
    return True


def centered_matrix(L: Matrix) -> np.ndarray:
    """Return ``E = L - (n+1)/2 * J`` with exact ``Fraction`` entries."""
    L = np.asarray(L, dtype=int)
    n = L.shape[0]
    m = Fraction(n + 1, 2)
    E = np.empty_like(L, dtype=object)
    for i in range(n):
        for j in range(n):
            E[i, j] = Fraction(int(L[i, j])) - m
    return E


def difference_matrix(L: Matrix) -> np.ndarray:
    """Return the ``(n-1) x (n-1)`` integer matrix ``A_{ij} = L_{ij} - L_{i,n-1}``."""
    L = np.asarray(L, dtype=int)
    n = L.shape[0]
    last = L[:, -1]
    A = np.zeros((n - 1, n - 1), dtype=object)
    for i in range(n - 1):
        for j in range(n - 1):
            A[i, j] = int(L[i, j]) - int(last[i])
    return A


def gram_matrix(n: int) -> np.ndarray:
    """Return the Gram matrix ``G = I_{n-1} + J_{n-1}``, ``det(G) = n``."""
    G = np.eye(n - 1, dtype=object) + np.ones((n - 1, n - 1), dtype=object)
    return G


def gram_projected(L: Matrix) -> np.ndarray:
    """Return ``E_std = G * A`` (the Gram-projected matrix)."""
    A = difference_matrix(L)
    n = A.shape[0] + 1
    G = gram_matrix(n)
    return _imul(G, A)


def parity_pattern(L: Matrix) -> np.ndarray:
    """Return ``L mod 2`` as a 0/1 ``numpy`` array."""
    return (np.asarray(L, dtype=int) & 1).astype(np.int8)


# ---------------------------------------------------------------------------
# Bareiss algorithm over exact integers -- Sylvester's identity with
# integer-preserving Gaussian elimination (\cite{Bareiss68}).
# ---------------------------------------------------------------------------
def bareiss_det(M: np.ndarray) -> int:
    """Return ``det(M)`` exactly using the Bareiss algorithm.

    ``M`` must contain Python ``int`` / ``Fraction`` entries.  Runs in
    ``O(n^3)`` integer operations with no coefficient explosion.
    """
    M = np.asarray(M, dtype=object).copy()
    n, m = M.shape
    if n != m:
        raise ValueError("Bareiss determinant requires a square matrix")
    sign = 1
    prev = 1
    for k in range(n - 1):
        # Partial pivot if M[k,k] is zero.
        if M[k, k] == 0:
            for r in range(k + 1, n):
                if M[r, k] != 0:
                    M[[k, r]] = M[[r, k]]
                    sign = -sign
                    break
            else:
                return 0
        for i in range(k + 1, n):
            for j in range(k + 1, n):
                num = M[i, j] * M[k, k] - M[i, k] * M[k, j]
                # Exact division by the previous pivot (Bareiss guarantee).
                if isinstance(num, Fraction) or isinstance(prev, Fraction):
                    M[i, j] = num / prev
                else:
                    q, r = divmod(num, prev)
                    if r != 0:
                        M[i, j] = Fraction(num, prev)
                    else:
                        M[i, j] = q
            M[i, k] = 0
        prev = M[k, k]
    result = sign * M[n - 1, n - 1]
    if isinstance(result, Fraction):
        if result.denominator != 1:
            raise ArithmeticError(
                f"Non-integer determinant {result}; check matrix entries"
            )
        result = result.numerator
    return int(result)


def _imul(A: np.ndarray, B: np.ndarray) -> np.ndarray:
    """Exact matrix multiplication preserving the object dtype."""
    A = np.asarray(A, dtype=object)
    B = np.asarray(B, dtype=object)
    n, k = A.shape
    k2, m = B.shape
    assert k == k2
    out = np.zeros((n, m), dtype=object)
    for i in range(n):
        for j in range(m):
            s = 0
            for t in range(k):
                s = s + A[i, t] * B[t, j]
            out[i, j] = s
    return out


def det_A(L: Matrix) -> int:
    """Return ``det(A)`` where ``A`` is the difference matrix of ``L``."""
    return bareiss_det(difference_matrix(L))


def det_Estd(L: Matrix) -> int:
    """Return ``det(E_std) = n * det(A)``."""
    n = np.asarray(L).shape[0]
    return n * det_A(L)
