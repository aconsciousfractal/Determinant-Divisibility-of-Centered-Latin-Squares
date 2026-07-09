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

## Certified Package (built 2026-07-09, closing GAP-1)

The curated package promised by the paper's "Code and data availability"
paragraph (Section 8) was built on 2026-07-09 (owner-directed). It maps
to the paragraph as follows; the four documented deviations (D1-D4) are
in `CLAIM_LEDGER.md` GAP-1.

| Promised item | Artifact | Generator | Verified 2026-07-09 |
| --- | --- | --- | --- |
| n=5 exhaustive | `results/proof_n5_data.json` | research-phase (all 56 reduced squares) | pinned in SHA256SUMS |
| n=6 exhaustive | `results/n6_results.txt` | `scripts/sandpile_n6_exhaustive.py` | replayed green 2026-07-08 |
| n=10 switch-chain table (4,000 samples) | `results/phase12/` (4 files) | `scripts/phase12_corpus.py` | replay BYTE-IDENTICAL (35 s); real seed 20260412, paper's "seed 42" is incorrect |
| n=8 table (10,000 samples, seed 0) | `results/certified/jm_n8_corpus_summary.json` | `scripts/jm_n8_corpus.py` | REGENERATED to paper spec (original lost); 10.56% sharp vs paper's "approximately 11%"; 0 bound violations |
| cyclic SNF table n ≤ 12 | `results/certified/cyclic_snf_table.json` | `scripts/cyclic_snf_table.py` | 11/11 rows match diag(1,n,...,n) |
| p-adic bound scan | `results/certified/padic_scan_output.txt` | `scripts/padic_bound_scan.py` | 2700 (n,p,sample) tests, n = 6..20 (NOT exactly the promised n ≤ 12 — deviation D4 in `CLAIM_LEDGER.md`), 0 violations |
| ex:n10 witness JSON | `results/certified/ex_n10_witness.json` | `scripts/verify_isotopy_destruction.py` | 12/12 script checks PASS (the paper table's 8 invariants + the Discussion destruction) |
| isotopy script | `scripts/verify_isotopy_destruction.py` | (standalone, deterministic) | standard-form reduction gives EVEN det(A) = 3554320 — counterexample destroyed, as the Discussion claims |
| SHA-256 manifest | `results/certified/SHA256SUMS` | `scripts/make_sha256sums.py` | see below |

Replay commands (all deterministic; Windows: keep `PYTHONUTF8=1`):

```bash
PYTHONUTF8=1 python scripts/verify_isotopy_destruction.py   # ~1 s
PYTHONUTF8=1 python scripts/cyclic_snf_table.py             # ~1 s
PYTHONUTF8=1 python scripts/padic_bound_scan.py > results/certified/padic_scan_output.txt   # ~2 min; stdout IS the artifact — use bash, not PowerShell (PowerShell redirection writes UTF-16 and breaks the manifest)
PYTHONUTF8=1 python scripts/phase12_corpus.py               # ~35 s, rewrites results/phase12/ byte-identically
PYTHONUTF8=1 python scripts/jm_n8_corpus.py                 # ~3 min
```

Manifest check from the repo root:

```bash
sha256sum -c results/certified/SHA256SUMS
```
