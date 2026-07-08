# Reproduce

Companion documentation added 2026-07-08, after publication (tag
`v1.0-submission`); it does not modify the paper.

## Environment

- Python 3.11+ with `numpy>=1.24`, `sympy>=1.12` (`pip install -r requirements.txt`).
- On Windows, set `PYTHONUTF8=1` (some scripts print Unicode symbols and
  crash on the default cp1252 console encoding otherwise).

## Smoke Replay (verified 2026-07-08)

```bash
PYTHONUTF8=1 python scripts/sandpile_n6_exhaustive.py
```

Exhaustively analyzes all 9408 reduced Latin squares of order 6.
Observed output (2026-07-08, ~1-2 min; excerpt, lines reordered, ASCII transliteration):

```text
Reduced LS with 36 | det(E_std): 8832
Reduced LS with 18 | det(E_std) but 36 not| det(E_std): 576
Distinct |det|: 197
Results saved to n6_results.txt
```

The 576 squares failing 36-divisibility at n=6 are the paper's
n ≡ 2 (mod 4) counterexample phenomenon at the smallest order
(paper Section 7.1, exhaustive classification at n=6).

## Long-Tier Checks (documented, not re-run in the 2026-07-08 retrofit)

These are Jacobson-Matthews / switch-chain sampling runs whose runtime
exceeds a bounded session (both exceeded a 2-4 minute cap when attempted on
2026-07-08 and were intentionally not run to completion):

```bash
PYTHONUTF8=1 python scripts/jm_paper_verify.py    # Sections 3-5 sampling cross-check
PYTHONUTF8=1 python scripts/verify_n8_jm.py       # n=8 witness: 32|det, 64 not|det
```

Sampling runs are consistency checks and witness generators only; they are
not proofs (see paper Section 8, Computational methods).

## Known Gap (recorded 2026-07-08)

The paper's "Code and data availability" paragraph (Section 8) describes a
curated `latin_det` package with certified datasets, seeds,
`scripts/verify_isotopy_destruction.py`, and a
`results/certified/SHA256SUMS` manifest. Those curated artifacts are not
present in this repository snapshot; the repository currently contains the
research-phase scripts and captured outputs. See
`CLAIM_LEDGER.md` (row C9/GAP-1) and `PUBLIC_CLAIM_BOUNDARY.md`. Closing
this gap (building the promised certified package) is an owner decision.
