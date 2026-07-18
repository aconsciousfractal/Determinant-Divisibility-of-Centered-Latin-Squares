# Claim Ledger

Companion documentation added 2026-07-08, after publication (tag
`v1.0-submission`); it does not modify the paper. Source locators refer to
`paper/Determinant Divisibility of Centered Latin Squares.tex` (LaTeX labels
in parentheses). Levels: `theorem-in-paper` (proved in the paper),
`certified-finite` (exhaustive exact computation), `evidence` (sampling or
bounded computation, not proof), `blocked` (promised artifact missing),
`resolved` (a previously blocked gap, now closed — row kept for the audit
trail).

The added `CL` column uses the official PAPP taxonomy
(`Gate-Disciplined-Computational-Mathematics/docs/CLAIM_LEVELS.md`); added
2026-07-17 (ledger A-1). It maps the existing bespoke `Level` words without
changing any claim: `theorem-in-paper`→`CL5` (internal theorem),
`certified-finite`→`CL3` (certified finite result), `evidence`→`CL1`
(computational evidence), `blocked`→`CLB`, `resolved`→`CLM` (metadata /
audit-trail row). Existing `Level` words, IDs, statements, locators, and
scope caveats are unchanged.

| ID | Level | CL | Statement (scoped) | Source locator | Evidence in repo |
| --- | --- | --- | --- | --- | --- |
| C1 | theorem-in-paper | CL5 | Unified divisibility: n²/gcd(n,2) divides det(E_std) for every Latin square of order n ≥ 2; for odd n, n² divides det(E_std). | Sec. 3 (`thm:unified`) | Proof in paper; exhaustive n=4 (576 squares) and n=6 (9408) cross-checks; n=6 replayed 2026-07-08 (`REPRODUCE.md`). |
| C2 | theorem-in-paper | CL5 | F2-rank characterization, n ≡ 2 (mod 4): n² divides det(E_std) iff rank_F2(A mod 2) < n−1. | Sec. 4 (`thm:F2-characterization`) | Proof in paper; `verify_unified_criterion.py`, `verify_even_n.py` (research-phase checks). |
| C3 | theorem-in-paper | CL5 | Adjugate criterion, n ≡ 0 (mod 4): the improved divisibility holds iff adj(A mod 2)·1 = 0 over F2. | Sec. 4 (`thm:adjugate-criterion`) | Proof in paper. |
| C4 | theorem-in-paper | CL5 | The parity pattern L mod 2 is a complete invariant for the n²-divisibility obstruction (obstruction is entirely 2-adic). | Sec. 4 (`cor:parity-complete`, even n) and Sec. 5 (`cor:obstruction`) | Proof in paper. |
| C5 | theorem-in-paper | CL5 | Universal counterexamples: for every n ≡ 2 (mod 4), n ≥ 6, an explicit skip-one circulant parity family gives Latin squares with n² not dividing det(E_std). | Sec. 4.1 (`thm:universal-counterexamples`) | Proof in paper; n=6 instance visible in the exhaustive census (576 failing squares, replayed 2026-07-08). |
| C6 | theorem-in-paper | CL5 | p-adic lower bounds: v_p(det A) ≥ v_p(n)+max(0,k_p−1) for odd p dividing n; v_2(det A) ≥ v_2(n/2)+max(0,k_2−1) for even n. | Sec. 6 (`thm:sharp-v2`, `cor:sharp-vp`) | Proof in paper. |
| C7 | evidence | CL1 | TIGHTNESS of the C6 bounds across the sampled orders n ∈ {6,8,9,10,12,15,16,18,20}. The bounds are proved; their empirical tightness is sampling evidence only. | Abstract; Sec. 6 | `scripts/padic_bound_scan.py` re-run 2026-07-09: 2700 (n,p,sample) tests over the sampled orders n ∈ {6,8,9,10,12,15,16,18,20}, for the primes p | n (p ∈ {2,3,5}), 0 violations; equality with the bound is attained at every (n,p) pair, but NOT at every k_p (the artifact's per-k_p tables show non-tight tails, e.g. k_p = 4 rows; for n ≡ 2 mod 4, p = 2, k ≥ 1 tightness is impossible by the paper's own rank characterization). The paper's remark after `cor:sharp-vp` now states this at the (n,p) level (reworded 2026-07-15 from "equality was observed in every sampled case", which over-stated per-sample tightness); the scan's k_p-resolved tables are the finer record. Output pinned in `results/certified/padic_scan_output.txt`. |
| C8 | certified-finite | CL3 | n=6 exhaustive classification: all 9408 reduced Latin squares enumerated; 8832 satisfy 36-divisibility, 576 fail; orbit classification as stated. | Sec. 7.1 (`prop:n6`, `prop:n6-orbits`) | `scripts/sandpile_n6_exhaustive.py` → `results/n6_results.txt`; replayed green 2026-07-08. |
| C9 | evidence | CL1 | Cyclic-family formula and cokernel (sandpile) statements as stated in the paper; JM sampling consistency at n=8 (10,000 samples, seed 0) and switch-chain at n=10 (4,000 samples; paper says "seed 42", actual seed 20260412 — GAP-1/D1). | Sec. 3 (`prop:cyclic`), Sec. 3 (`cor:sandpile-cyclic`), Sec. 8 | Formula/cokernel: proof in paper (theorem-in-paper). Tables now certified: n=10 in `results/phase12/` (byte-identical replay 2026-07-09), n=8 regenerated in `results/certified/jm_n8_corpus_summary.json` (GAP-1/D2). |
| GAP-1 | resolved | CLM | The curated `latin_det` package promised in the paper's "Code and data availability" paragraph: certified datasets (n∈{5,6} exhaustive, n=8 seed-0 table, n=10 switch-chain table, cyclic SNF table n≤12, p-adic scan), `scripts/verify_isotopy_destruction.py`, and `results/certified/SHA256SUMS`. Was MISSING at the 2026-07-08 audit (RT-9); BUILT 2026-07-09, owner-directed — full item-by-item map in `REPRODUCE.md`. FOUR DOCUMENTED DEVIATIONS from the paragraph's wording: (D1) the n=10 table's actual seed is 20260412 (`scripts/phase12_corpus.py`, replay byte-identical to the tracked `results/phase12/`); the paper's "seed 42" is incorrect. (D2) the n=8 10,000-sample corpus is a REGENERATION to the paper's spec (10,000 JM samples, 320 steps, seed 0) — the original run did not survive; observed 10.56% sharp vs the paper's explicitly approximate "~11%"; additionally the paper-proof's witness selection rule ("first dim-ker-1 sample") is sampler-implementation-specific — the regeneration certifies a verified sharpness witness instead. (D3) "`latin_det` package" = this repository's `scripts/` + `results/` layout; no separately installable package exists. (D4) the delivered p-adic scan samples the nine orders n ∈ {6,8,9,10,12,15,16,18,20}, not exactly the promised n ≤ 12: orders 2..5 are absent from the sampled scan (n = 5 is covered exhaustively by `results/proof_n5_data.json`), and orders 15,16,18,20 are extra (beyond n ≤ 12); the paper's Sec. 6 remark now states this exact sampled set and the (n,p)-level tightness. | Sec. 8, final paragraph | Built and replay-verified 2026-07-09; manifest `results/certified/SHA256SUMS`. **Paper text reconciled 2026-07-15** (red-team audit): §8 "Code and data availability" now reads seed `20260412` (D1), `10.56%` sharp (D2), "public open-source repository (scripts/+results/ trees)" instead of a `latin_det` package (D3), and the p-adic bound scan stated over its exact sampled set `n∈{6,8,9,10,12,15,16,18,20}` (D4, corrected 2026-07-15 — the earlier "n≤20" over-declared the scan domain; the script tests only those 9 orders) — so the four deviations are now closed in the manuscript, not only documented here. |

Red-team notes: `scripts/red_team_step1.py` and `scripts/red_team_verify.py`
are the research-phase adversarial recompute scripts; they are part of the
evidence trail, not proofs.

## Errata

- **E1 (2026-07-18, paper corrected).** `rem:consecutive` (Sec. 4.1) stated
  that the consecutive-interval circulant satisfies
  `rank_F2(P_{T_n}) = n/2`, "the minimum possible for a doubly balanced
  circulant". Both clauses were wrong: the remark's own gcd computation
  (degree `k−1` drop from full rank `n`) gives `n − (k−1) = n/2 + 1`, and
  direct GF(2) elimination confirms ranks 4, 6, 8, 10, 12 at
  `n = 6, 10, 14, 18, 22`; nor is `T_n` rank-minimal among doubly balanced
  circulants (the period-two support `{0,2,…,n−2}` has rank 2). The remark
  now states `n/2 + 1` with the arithmetic explicit, and the minimality
  clause is removed. **No theorem is affected**: C5's mechanism uses the
  skip-one family `S_n` (full rank, proved separately); `T_n` remains
  degenerate (far below `n−1`) under either value.
- **E2 (2026-07-18, manifest repinned).** `results/certified/SHA256SUMS`
  (built 2026-07-09) pinned CRLF disk bytes from a Windows
  `core.autocrlf=true` checkout for 8 of its 10 entries, while git stores
  those text artifacts with LF endings — so `sha256sum -c` could never
  pass on a fresh clone. Verified file-by-file that every mismatching hash
  equals the CRLF-translation of the stored blob (**no artifact content
  changed**). The manifest is repinned to the canonical LF blob bytes and
  `scripts/make_sha256sums.py` now normalizes CRLF→LF before hashing so a
  regeneration is checkout-independent.
