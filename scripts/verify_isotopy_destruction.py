"""Verify the "isotopy destruction" claim of Section 7.1.

The paper asserts that reducing the n=10 counterexample of
Example~\\ref{ex:n10} to standard form (first row = (1,...,10)) via a
symbol permutation gamma produces a Latin square with even det(A),
showing that the "counterexample property" is not an isotopy invariant.

This script:

1. loads the n=10 witness Latin square;
2. applies the symbol permutation gamma that maps its first row to
   (1,2,...,10);
3. recomputes det(A) and v_2(det(A)) on the relabeled square;
4. additionally exhibits (if possible) a symbol permutation that
   restores v_2(det A) = 0, illustrating that the deficit depends on the
   *labeled* Latin square.

Output: ``results/certified/isotopy_destruction.json``.
"""

from __future__ import annotations

import itertools
import json
from pathlib import Path

import numpy as np

from latin_det.core import bareiss_det, difference_matrix, is_latin_square
from latin_det.snf import v_p


WITNESS = np.array([
    [ 2,  5,  6,  4,  3,  8,  7,  9, 10,  1],
    [10,  6,  2,  7,  5,  1,  9,  8,  4,  3],
    [ 3,  4,  1, 10,  8,  5,  2,  7,  9,  6],
    [ 7,  1, 10,  8,  6,  9,  3,  4,  2,  5],
    [ 1,  7,  3,  6,  4, 10,  8,  2,  5,  9],
    [ 9,  2,  5,  3, 10,  6,  4,  1,  8,  7],
    [ 6,  9,  8,  1,  2,  7,  5, 10,  3,  4],
    [ 5,  8,  7,  2,  9,  4,  1,  3,  6, 10],
    [ 4,  3,  9,  5,  1,  2, 10,  6,  7,  8],
    [ 8, 10,  4,  9,  7,  3,  6,  5,  1,  2],
])


def apply_symbol_perm(L: np.ndarray, gamma: np.ndarray) -> np.ndarray:
    """Relabel symbols: gamma[s - 1] is the new symbol replacing s."""
    return gamma[L - 1]


def main() -> None:
    n = 10
    assert is_latin_square(WITNESS)
    d0 = bareiss_det(difference_matrix(WITNESS))
    v2_0 = v_p(d0, 2)
    assert d0 == 15_427_045
    assert v2_0 == 0

    # gamma such that gamma(first row of WITNESS) = (1,...,10).
    first = WITNESS[0]  # [2, 5, 6, 4, 3, 8, 7, 9, 10, 1]
    gamma = np.zeros(n, dtype=int)
    for new_sym, s in enumerate(first, start=1):
        gamma[s - 1] = new_sym
    L_std = apply_symbol_perm(WITNESS, gamma)
    assert is_latin_square(L_std)
    assert list(L_std[0]) == list(range(1, n + 1))
    d_std = bareiss_det(difference_matrix(L_std))
    v2_std = v_p(d_std, 2)

    # Enumerate a handful of random symbol permutations to show the
    # distribution of v_2(det A) over the isotopy class.
    rng = np.random.default_rng(0)
    distribution: dict = {}
    tested = 2000
    for _ in range(tested):
        perm = rng.permutation(n) + 1
        Lp = apply_symbol_perm(WITNESS, perm)
        dp = bareiss_det(difference_matrix(Lp))
        v = v_p(dp, 2)
        key = "inf" if v == float("inf") else int(v)
        distribution[key] = distribution.get(key, 0) + 1

    result = {
        "source": "Example ex:n10 in paper",
        "witness_det_A": int(d0),
        "witness_v2": v2_0,
        "standardized_det_A": int(d_std),
        "standardized_v2": v2_std,
        "claim_destruction": bool(v2_std >= 1),
        "distribution_v2_over_2000_symbol_perms": {
            str(k): v for k, v in sorted(
                distribution.items(),
                key=lambda kv: (kv[0] == "inf", kv[0]),
            )
        },
        "count_v2_zero_in_isotopy_class_sample": int(distribution.get(0, 0)),
        "sample_size": tested,
    }
    out_dir = Path(__file__).resolve().parent.parent / "results" / "certified"
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / "isotopy_destruction.json"
    out_path.write_text(json.dumps(result, indent=2), encoding="utf-8")
    print(json.dumps(result, indent=2))
    print(f"\nWrote {out_path}")


if __name__ == "__main__":
    main()
