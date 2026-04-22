"""Linear algebra over ``F_2`` using ``numpy`` bitwise operations."""

from __future__ import annotations

from typing import List

import numpy as np


def _to_f2(M: np.ndarray) -> np.ndarray:
    return (np.asarray(M, dtype=np.int64) & 1).astype(np.int8)


def rank_f2(M: np.ndarray) -> int:
    """Return ``rank_{F_2}(M)`` via Gaussian elimination over ``F_2``."""
    A = _to_f2(M).copy()
    rows, cols = A.shape
    r = 0
    for c in range(cols):
        if r >= rows:
            break
        pivot = None
        for i in range(r, rows):
            if A[i, c] == 1:
                pivot = i
                break
        if pivot is None:
            continue
        if pivot != r:
            A[[r, pivot]] = A[[pivot, r]]
        for i in range(rows):
            if i != r and A[i, c] == 1:
                A[i] ^= A[r]
        r += 1
    return r


def ker_f2(M: np.ndarray) -> np.ndarray:
    """Return a basis of ``ker_{F_2}(M)`` as a ``(dim, cols)`` matrix."""
    A = _to_f2(M).copy()
    rows, cols = A.shape
    # Reduce to RREF.
    pivots: List[int] = []
    r = 0
    for c in range(cols):
        if r >= rows:
            break
        pivot = None
        for i in range(r, rows):
            if A[i, c] == 1:
                pivot = i
                break
        if pivot is None:
            continue
        if pivot != r:
            A[[r, pivot]] = A[[pivot, r]]
        for i in range(rows):
            if i != r and A[i, c] == 1:
                A[i] ^= A[r]
        pivots.append(c)
        r += 1
    free = [c for c in range(cols) if c not in pivots]
    basis = []
    for f in free:
        v = np.zeros(cols, dtype=np.int8)
        v[f] = 1
        for pi, pc in enumerate(pivots):
            v[pc] = A[pi, f]
        basis.append(v)
    if not basis:
        return np.zeros((0, cols), dtype=np.int8)
    return np.vstack(basis)


def adjugate_f2(B: np.ndarray) -> np.ndarray:
    """Return ``adj(B) mod 2`` as an ``n x n`` ``F_2`` matrix.

    Uses the minor definition: ``adj(B)_{ij} = (-1)^{i+j} * det(M_{ji})``
    where ``M_{ji}`` is the ``(n-1) x (n-1)`` submatrix obtained by deleting
    row ``j`` and column ``i``.  Signs are irrelevant over ``F_2``.
    """
    A = _to_f2(B).copy()
    n = A.shape[0]
    adj = np.zeros((n, n), dtype=np.int8)
    for i in range(n):
        for j in range(n):
            minor = np.delete(np.delete(A, j, axis=0), i, axis=1)
            adj[i, j] = _det_f2(minor)
    return adj


def _det_f2(M: np.ndarray) -> int:
    A = _to_f2(M).copy()
    n = A.shape[0]
    for c in range(n):
        pivot = None
        for i in range(c, n):
            if A[i, c] == 1:
                pivot = i
                break
        if pivot is None:
            return 0
        if pivot != c:
            A[[c, pivot]] = A[[pivot, c]]
        for i in range(c + 1, n):
            if A[i, c] == 1:
                A[i] ^= A[c]
    return 1
