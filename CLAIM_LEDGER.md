# Claim Ledger

Companion documentation added 2026-07-08, after publication (tag
`v1.0-submission`); it does not modify the paper. Source locators refer to
`paper/Determinant Divisibility of Centered Latin Squares.tex` (LaTeX labels
in parentheses). Levels: `theorem-in-paper` (proved in the paper),
`certified-finite` (exhaustive exact computation), `evidence` (sampling or
bounded computation, not proof), `blocked` (promised artifact missing).

| ID | Level | Statement (scoped) | Source locator | Evidence in repo |
| --- | --- | --- | --- | --- |
| C1 | theorem-in-paper | Unified divisibility: n²/gcd(n,2) divides det(E_std) for every Latin square of order n ≥ 2; for odd n, n² divides det(E_std). | Sec. 3 (`thm:unified`) | Proof in paper; exhaustive n=4 (576 squares) and n=6 (9408) cross-checks; n=6 replayed 2026-07-08 (`REPRODUCE.md`). |
| C2 | theorem-in-paper | F2-rank characterization, n ≡ 2 (mod 4): n² divides det(E_std) iff rank_F2(A mod 2) < n−1. | Sec. 4 (`thm:F2-characterization`) | Proof in paper; `verify_unified_criterion.py`, `verify_even_n.py` (research-phase checks). |
| C3 | theorem-in-paper | Adjugate criterion, n ≡ 0 (mod 4): the improved divisibility holds iff adj(A mod 2)·1 = 0 over F2. | Sec. 4 (`thm:adjugate-criterion`) | Proof in paper. |
| C4 | theorem-in-paper | The parity pattern L mod 2 is a complete invariant for the n²-divisibility obstruction (obstruction is entirely 2-adic). | Sec. 4 (`cor:parity-complete`, even n) and Sec. 5 (`cor:obstruction`) | Proof in paper. |
| C5 | theorem-in-paper | Universal counterexamples: for every n ≡ 2 (mod 4), n ≥ 6, an explicit skip-one circulant parity family gives Latin squares with n² not dividing det(E_std). | Sec. 4.1 (`thm:universal-counterexamples`) | Proof in paper; n=6 instance visible in the exhaustive census (576 failing squares, replayed 2026-07-08). |
| C6 | theorem-in-paper | p-adic lower bounds: v_p(det A) ≥ v_p(n)+max(0,k_p−1) for odd p dividing n; v_2(det A) ≥ v_2(n/2)+max(0,k_2−1) for even n. | Sec. 6 (`thm:sharp-v2`, `cor:sharp-vp`) | Proof in paper. |
| C7 | evidence | TIGHTNESS of the C6 bounds across sampled instances up to n ≤ 20. The bounds are proved; their empirical tightness is sampling evidence only. | Abstract; Sec. 6 | JM sampling tables (research-phase outputs under `results/`); not re-run 2026-07-08. |
| C8 | certified-finite | n=6 exhaustive classification: all 9408 reduced Latin squares enumerated; 8832 satisfy 36-divisibility, 576 fail; orbit classification as stated. | Sec. 7.1 (`prop:n6`, `prop:n6-orbits`) | `scripts/sandpile_n6_exhaustive.py` → `results/n6_results.txt`; replayed green 2026-07-08. |
| C9 | evidence | Cyclic-family formula and cokernel (sandpile) statements as stated in the paper; JM sampling consistency at n=8 (10,000 samples, seed 0) and switch-chain at n=10 (4,000 samples, seed 42). | Sec. 3 (`prop:cyclic`), Sec. 3 (`cor:sandpile-cyclic`), Sec. 8 | Formula/cokernel: proof in paper (theorem-in-paper). The named seeded sample TABLES are not present as certified datasets in this snapshot — see GAP-1. |
| GAP-1 | blocked | The curated `latin_det` package promised in the paper's "Code and data availability" paragraph: certified datasets (n∈{5,6} exhaustive, n=8 seed-0 table, n=10 seed-42 table, cyclic SNF table n≤12, p-adic scan), `scripts/verify_isotopy_destruction.py`, and `results/certified/SHA256SUMS`. | Sec. 8, final paragraph | MISSING from this repository snapshot (verified 2026-07-08). Building the promised package is an owner decision; until then the paper's availability paragraph overstates the repo contents (RT-9). |

Red-team notes: `scripts/red_team_step1.py` and `scripts/red_team_verify.py`
are the research-phase adversarial recompute scripts; they are part of the
evidence trail, not proofs.
