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
  sampling evidence on the nine sampled orders
  n ∈ {6,8,9,10,12,15,16,18,20} only (C7).
- That Jacobson-Matthews or switch-chain samples are exactly uniform — the
  paper explicitly does not claim this (Sec. 8).
- That sampling tables prove any theorem — they are consistency checks and
  witness sources only.
- That the certified package matches the paper's availability paragraph
  word-for-word. The package was BUILT 2026-07-09 (owner-directed) and
  the ledger's GAP-1 row records four documented deviations that may
  not be papered over when quoting it: (1) the n=10 switch-chain table
  was generated with seed 20260412, not the paper's "seed 42"; (2) the
  n=8 10,000-sample corpus is a REGENERATION to the paper's spec (the
  original did not survive), not the original run; (3) "the latin_det
  package" is this repository's scripts/ + results/ layout, not a
  separately installable Python package; (4) the p-adic scan samples
  n ∈ {6,8,9,10,12,15,16,18,20}, not exactly the promised n ≤ 12.

## Scope Notes

- All theorem claims are about the centered operator E = L − ((n+1)/2)J on
  the standard subspace, per the paper's conventions (Sec. 2).
- Nothing here extends the claims beyond the paper; this file only fences
  how the existing claims may be quoted.
