"""Samplers for Latin squares and doubly balanced binary matrices.

This module provides two reproducible samplers used to generate the
empirical tables referenced in Sections 5-7:

* :func:`random_latin_square` -- row-wise random completion with
  backtracking.  Each call returns a *valid* Latin square; the
  distribution is **not** claimed to be uniform.  The paper's empirical
  remarks (Remark~5.X, Section~7) explicitly flag this.  The function
  seeds a dedicated :class:`numpy.random.Generator` so outputs are
  bit-reproducible.

* :func:`switch_chain_balanced` -- Markov chain on doubly balanced
  ``n x n`` binary matrices; moves are 2x2 checkerboard flips.  Samples
  preserve all row and column sums.

Back-compat shims (:func:`jacobson_matthews`,
:func:`jacobson_matthews_samples`) keep the earlier API working; they
now delegate to :func:`random_latin_square`.
"""

from __future__ import annotations

from typing import Iterator, List, Optional

import numpy as np


# ---------------------------------------------------------------------------
# Random Latin square via row-wise backtracking.
# ---------------------------------------------------------------------------
def random_latin_square(n: int, seed: int = 0) -> np.ndarray:
    """Return one random ``n x n`` Latin square with entries ``{1, ..., n}``.

    The algorithm fills the square row by row; within each row the
    columns are processed left to right, drawing the symbol uniformly
    from the set of symbols still legal at that cell.  If a row
    dead-ends, the entire square is restarted with a fresh attempt.
    """
    rng = np.random.default_rng(seed)
    return _random_latin_square_from_rng(n, rng)


def random_latin_squares(n: int, samples: int, seed: int = 0) -> Iterator[np.ndarray]:
    """Yield ``samples`` independent random Latin squares of order ``n``."""
    rng = np.random.default_rng(seed)
    for _ in range(samples):
        yield _random_latin_square_from_rng(n, rng)


def _random_latin_square_from_rng(n: int, rng: np.random.Generator) -> np.ndarray:
    for _attempt in range(10_000):
        L = np.zeros((n, n), dtype=int)
        col_used: List[set] = [set() for _ in range(n)]
        ok = True
        for r in range(n):
            row = _fill_row(n, col_used, rng)
            if row is None:
                ok = False
                break
            for c, v in enumerate(row):
                L[r, c] = v
                col_used[c].add(v)
        if ok:
            return L
    raise RuntimeError(
        f"Failed to construct a random Latin square of order {n} after 10000 tries"
    )


def _fill_row(n: int, col_used: List[set], rng: np.random.Generator) -> Optional[List[int]]:
    row: List[int] = [0] * n
    used_in_row: set = set()

    def recurse(c: int) -> bool:
        if c == n:
            return True
        options = [v for v in range(1, n + 1) if v not in col_used[c] and v not in used_in_row]
        if not options:
            return False
        rng.shuffle(options)
        for v in options:
            row[c] = v
            used_in_row.add(v)
            if recurse(c + 1):
                return True
            used_in_row.discard(v)
            row[c] = 0
        return False

    if recurse(0):
        return row
    return None


# ---------------------------------------------------------------------------
# Back-compat aliases matching the Section 7 naming in earlier drafts.
# ---------------------------------------------------------------------------
def jacobson_matthews(n: int, steps: int | None = None, seed: int = 0) -> np.ndarray:
    """Back-compat shim.  Returns a random Latin square; see
    :func:`random_latin_square`.  The ``steps`` argument is ignored.
    """
    del steps
    return random_latin_square(n, seed=seed)


def jacobson_matthews_samples(
    n: int,
    samples: int,
    steps_per_sample: int | None = None,
    seed: int = 0,
) -> Iterator[np.ndarray]:
    """Back-compat alias for :func:`random_latin_squares`."""
    del steps_per_sample
    return random_latin_squares(n, samples, seed=seed)


# ---------------------------------------------------------------------------
# Switch-chain MCMC on doubly balanced binary matrices R(n, n/2).
# ---------------------------------------------------------------------------
def switch_chain_balanced(
    n: int,
    samples: int,
    burnin: int = 50_000,
    thin: int = 200,
    seed: int = 0,
) -> Iterator[np.ndarray]:
    """Yield ``samples`` doubly balanced ``n x n`` binary matrices.

    Initial state: circulant indicator of ``{0, 1, ..., n/2 - 1}`` shifted
    cyclically (deterministic starting point for reproducibility).
    """
    if n % 2 != 0:
        raise ValueError("switch_chain_balanced requires n even")
    k = n // 2
    rng = np.random.default_rng(seed)
    P = np.zeros((n, n), dtype=np.int8)
    for i in range(n):
        for s in range(k):
            P[i, (i + s) % n] = 1
    for _ in range(burnin):
        _switch_step(P, rng)
    for _ in range(samples):
        for _ in range(thin):
            _switch_step(P, rng)
        yield P.copy()


def _switch_step(P: np.ndarray, rng: np.random.Generator) -> None:
    n = P.shape[0]
    for _ in range(10):
        i1, i2 = rng.choice(n, size=2, replace=False)
        j1, j2 = rng.choice(n, size=2, replace=False)
        a, b = P[i1, j1], P[i1, j2]
        c, d = P[i2, j1], P[i2, j2]
        if (a, b, c, d) == (1, 0, 0, 1):
            P[i1, j1] = 0
            P[i1, j2] = 1
            P[i2, j1] = 1
            P[i2, j2] = 0
            return
        if (a, b, c, d) == (0, 1, 1, 0):
            P[i1, j1] = 1
            P[i1, j2] = 0
            P[i2, j1] = 0
            P[i2, j2] = 1
            return
