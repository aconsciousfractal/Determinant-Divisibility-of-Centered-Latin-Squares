"""latin_det — reproducibility toolkit for *Determinant Divisibility of
Centered Latin Squares* (Babanskyy, 2026).

Each submodule implements one deductively-grounded routine used in the
paper.  Deterministic seeds are propagated through every sampler so that
the CSV/JSON datasets under ``results/certified/`` reproduce bit-for-bit.
"""

from .core import (
    centered_matrix,
    difference_matrix,
    gram_matrix,
    gram_projected,
    bareiss_det,
    is_latin_square,
    parity_pattern,
)
from .f2 import rank_f2, ker_f2, adjugate_f2
from .snf import smith_normal_form, invariant_factors
from .samplers import jacobson_matthews, switch_chain_balanced
from .lifts import konig_lift, mrv_backtrack
from .enumerate import enumerate_reduced_latin

__version__ = "1.0.0"

__all__ = [
    "centered_matrix",
    "difference_matrix",
    "gram_matrix",
    "gram_projected",
    "bareiss_det",
    "is_latin_square",
    "parity_pattern",
    "rank_f2",
    "ker_f2",
    "adjugate_f2",
    "smith_normal_form",
    "invariant_factors",
    "jacobson_matthews",
    "switch_chain_balanced",
    "konig_lift",
    "mrv_backtrack",
    "enumerate_reduced_latin",
]
