"""Exhaustive enumeration of reduced Latin squares.

A *reduced* Latin square has first row = ``(1, 2, ..., n)`` and first
column = ``(1, 2, ..., n)``.  The enumeration uses a straightforward
backtracking with row/column/symbol pruning.
"""

from __future__ import annotations

from typing import Iterator, List

import numpy as np


def enumerate_reduced_latin(n: int) -> Iterator[np.ndarray]:
    """Yield every reduced Latin square of order ``n``."""
    L = np.zeros((n, n), dtype=int)
    for j in range(n):
        L[0, j] = j + 1
    for i in range(1, n):
        L[i, 0] = i + 1
    row_used: List[set] = [set(L[i].tolist()) - {0} for i in range(n)]
    col_used: List[set] = [set(L[:, j].tolist()) - {0} for j in range(n)]
    yield from _fill(L, 1, 1, row_used, col_used, n)


def _fill(L: np.ndarray, r: int, c: int, row_used, col_used, n: int) -> Iterator[np.ndarray]:
    if r >= n:
        yield L.copy()
        return
    nr, nc = (r, c + 1) if c + 1 < n else (r + 1, 1)
    for v in range(1, n + 1):
        if v in row_used[r] or v in col_used[c]:
            continue
        L[r, c] = v
        row_used[r].add(v)
        col_used[c].add(v)
        yield from _fill(L, nr, nc, row_used, col_used, n)
        row_used[r].discard(v)
        col_used[c].discard(v)
        L[r, c] = 0
