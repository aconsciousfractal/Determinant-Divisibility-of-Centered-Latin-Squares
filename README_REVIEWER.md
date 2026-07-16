# Reviewer Guide

Companion documentation added 2026-07-08, after publication (tag
`v1.0-submission`); it does not modify the paper.

## Ten-Minute Path

1. Read the abstract in
   `paper/Determinant Divisibility of Centered Latin Squares.pdf`.
2. Read `CLAIM_LEDGER.md` — three theorem groups (C1-C6), one certified
   finite census (C8), sampling evidence (C7/C9), and the resolved
   package gap with its four documented deviations (GAP-1).
3. Run the smoke replay:
   `PYTHONUTF8=1 python scripts/sandpile_n6_exhaustive.py`
   (~1-2 min; expects 9408 = 8832 + 576 at order 6).

## Thirty-Minute Path

4. Read Sections 3-4 of the paper (unified divisibility; binary criterion)
   against `scripts/jm_paper_verify.py`'s docstring, which restates the
   claims being sampled.
5. Skim Section 7 (counterexamples) and compare with the 576 failing
   squares in `results/n6_results.txt`.
6. Read `PUBLIC_CLAIM_BOUNDARY.md` for what may and may not be quoted.

## Main Claims

- Unified divisibility theorem (all orders), F2-rank/adjugate complete
  binary criterion (even orders), universal counterexample family
  (n ≡ 2 mod 4, n ≥ 6) — proved in the paper.
- p-adic lower bounds proved; their tightness is sampled evidence only.
- n=6: exhaustive, replayable census (certified finite).

## Known Limits

- Sampling (JM, switch-chain) is not exact-uniform and proves nothing.
- Long sampling scripts are multi-minute-to-hour runs; the n=6 census
  was re-verified 2026-07-08; the certified-package generators were all
  run/replayed 2026-07-09 (`REPRODUCE.md`).
- GAP-1 (resolved 2026-07-09): the certified package now exists
  (`results/certified/`, SHA256SUMS) with four documented deviations
  from the paper's availability paragraph — the n=10 table's real seed
  is 20260412 (not the stated "seed 42"), the n=8 corpus is a
  regeneration to spec (original lost; 10.56% observed vs "~11%"),
  "latin_det" denotes this repo's layout, not an installable package,
  and the p-adic scan samples n ∈ {6,8,9,10,12,15,16,18,20} rather than
  exactly the promised n ≤ 12 (`CLAIM_LEDGER.md`).
