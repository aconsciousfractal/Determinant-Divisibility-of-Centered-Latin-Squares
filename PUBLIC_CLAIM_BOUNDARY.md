# Public Claim Boundary

Companion documentation added 2026-07-08, after publication (tag
`v1.0-submission`); it does not modify the paper.

## Can Say

- The paper proves the unified divisibility theorem, the complete
  F2-rank/adjugate binary criterion for even orders, and the existence of
  universal counterexamples for every n ≡ 2 (mod 4), n ≥ 6
  (ledger C1-C5).
- The p-adic lower bounds are proved (C6).
- The n=6 classification is an exhaustive exact computation over all 9408
  reduced Latin squares and replays from this repository (C8).

## Must Not Say

- That the p-adic bounds are proved TIGHT in general — tightness is
  sampling evidence up to n ≤ 20 only (C7).
- That Jacobson-Matthews or switch-chain samples are exactly uniform — the
  paper explicitly does not claim this (Sec. 8).
- That sampling tables prove any theorem — they are consistency checks and
  witness sources only.
- That this repository currently contains the curated certified package
  (`latin_det`, certified datasets, SHA256SUMS manifest, isotopy script)
  described in the paper's availability paragraph — it does not (GAP-1),
  until the owner closes that gap.

## Scope Notes

- All theorem claims are about the centered operator E = L − ((n+1)/2)J on
  the standard subspace, per the paper's conventions (Sec. 2).
- Nothing here extends the claims beyond the paper; this file only fences
  how the existing claims may be quoted.
