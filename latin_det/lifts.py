"""Lifts from a parity pattern to a Latin square.

Two constructions are provided:

* :func:`konig_lift` -- edge-coloring via K\"onig's theorem (constructive
  perfect-matching decomposition of the bipartite graph).  This is the
  lift used in Proposition 4.4 (Balanced realizability).
* :func:`mrv_backtrack` -- backtracking with minimum-remaining-values
  heuristic, useful when a specific labeled Latin square is needed.

Both return integer ``n x n`` matrices with entries in ``{1, ..., n}``.
"""

from __future__ import annotations

from typing import List, Optional

import numpy as np


def _perfect_matchings(P: np.ndarray) -> List[np.ndarray]:
    """Decompose the ``k``-regular bipartite graph with biadjacency ``P``
    into ``k`` perfect matchings.  Uses a naive augmenting-path routine;
    sufficient for ``n <= 50``.
    """
    P = P.astype(np.int8).copy()
    n = P.shape[0]
    k = int(P[0].sum())
    matchings: List[np.ndarray] = []
    for _ in range(k):
        match_cols = [-1] * n
        used_row = [False] * n
        used_col = [False] * n

        def try_match(r: int, visited: List[bool]) -> bool:
            for c in range(n):
                if P[r, c] == 1 and not visited[c]:
                    visited[c] = True
                    if match_cols[c] == -1:
                        match_cols[c] = r
                        return True
                    else:
                        # Try to reassign the previously-matched row.
                        if try_match(match_cols[c], visited):
                            match_cols[c] = r
                            return True
            return False

        # Maximum matching (Kuhn's algorithm).
        for r in range(n):
            visited = [False] * n
            try_match(r, visited)
        # Build the matching matrix and subtract it from P.
        M = np.zeros((n, n), dtype=np.int8)
        for c, r in enumerate(match_cols):
            if r == -1:
                raise RuntimeError(
                    "P is not k-regular or decomposition failed; "
                    "perfect matching not found"
                )
            M[r, c] = 1
            P[r, c] -= 1
        matchings.append(M)
    return matchings


def konig_lift(P: np.ndarray) -> np.ndarray:
    """Lift a doubly balanced ``n x n`` binary matrix ``P`` to a Latin
    square ``L`` with ``L mod 2 == P``.

    The construction: decompose ``P`` into ``n/2`` matchings (labels
    ``1, 3, 5, ...``) and ``J - P`` into ``n/2`` matchings (labels
    ``2, 4, 6, ...``).
    """
    P = np.asarray(P, dtype=np.int8)
    n = P.shape[0]
    if P.shape != (n, n):
        raise ValueError("P must be square")
    row_sum = P.sum(axis=1)
    col_sum = P.sum(axis=0)
    if not np.all(row_sum == n // 2) or not np.all(col_sum == n // 2):
        raise ValueError("P must be doubly balanced with row/col sum n/2")
    odd_matchings = _perfect_matchings(P)
    even_matchings = _perfect_matchings((1 - P).astype(np.int8))
    L = np.zeros((n, n), dtype=int)
    for r, M in enumerate(odd_matchings):
        L[M == 1] = 2 * r + 1
    for s, M in enumerate(even_matchings):
        L[M == 1] = 2 * s + 2
    return L


def mrv_backtrack(P: np.ndarray, rng: Optional[np.random.Generator] = None) -> Optional[np.ndarray]:
    """Complete a parity pattern ``P`` to a Latin square via MRV backtracking.

    Returns ``None`` if no completion exists (which cannot happen for
    doubly balanced ``P``: see Proposition 4.4).
    """
    P = np.asarray(P, dtype=np.int8)
    n = P.shape[0]
    rng = np.random.default_rng() if rng is None else rng

    L = np.zeros((n, n), dtype=int)
    row_used = [set() for _ in range(n)]
    col_used = [set() for _ in range(n)]

    def domain(r: int, c: int) -> List[int]:
        parity = int(P[r, c])
        all_vals = [v for v in range(1, n + 1) if v % 2 == parity]
        return [v for v in all_vals if v not in row_used[r] and v not in col_used[c]]

    def select_cell() -> Optional[tuple]:
        best = None
        best_d = None
        for r in range(n):
            for c in range(n):
                if L[r, c] == 0:
                    d = len(domain(r, c))
                    if d == 0:
                        return (r, c)
                    if best is None or d < best_d:
                        best, best_d = (r, c), d
        return best

    def backtrack() -> bool:
        cell = select_cell()
        if cell is None:
            return True
        r, c = cell
        dom = domain(r, c)
        rng.shuffle(dom)
        for v in dom:
            L[r, c] = v
            row_used[r].add(v)
            col_used[c].add(v)
            if backtrack():
                return True
            L[r, c] = 0
            row_used[r].discard(v)
            col_used[c].discard(v)
        return False

    if backtrack():
        return L
    return None
