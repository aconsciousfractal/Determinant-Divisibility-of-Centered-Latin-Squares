"""Smith Normal Form utilities (integer)."""

from __future__ import annotations

from typing import List, Tuple

import numpy as np

try:
    from sympy import Matrix as _SymMatrix
    from sympy.matrices.normalforms import smith_normal_form as _sym_snf
except ImportError as exc:  # pragma: no cover
    raise SystemExit(
        "latin_det.snf requires sympy>=1.12; install via `pip install sympy`"
    ) from exc


def smith_normal_form(M: np.ndarray) -> np.ndarray:
    """Return the Smith Normal Form of the integer matrix ``M``.

    Wraps :func:`sympy.matrices.normalforms.smith_normal_form`.
    """
    M = np.asarray(M, dtype=object)
    S = _sym_snf(_SymMatrix(M.tolist()))
    return np.array(S.tolist(), dtype=object)


def invariant_factors(M: np.ndarray) -> List[int]:
    """Return the diagonal invariant factors ``d_1 | d_2 | ... | d_r``."""
    S = smith_normal_form(M)
    r = min(S.shape)
    return [int(S[i, i]) for i in range(r)]


def v_p(x: int, p: int) -> float:
    """Return ``v_p(x)``; ``v_p(0)`` is ``+inf``."""
    if x == 0:
        return float("inf")
    v = 0
    x = abs(int(x))
    while x % p == 0:
        x //= p
        v += 1
    return v
