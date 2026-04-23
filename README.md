# Determinant Divisibility of Centered Latin Squares

Companion repository for the paper

> **Determinant Divisibility of Centered Latin Squares: A Unified Theorem and the F₂-Rank Obstruction**
> Oleksiy Babanskyy, 2026.

Contains the manuscript source, the experimental scripts used to produce
and cross-check every empirical claim, and the datasets those scripts
generated during the research phase.

## Layout

```
.
├── paper/                     # LaTeX source (.tex), bibliography (.bbl), compiled PDF
├── scripts/                   # Self-contained Python scripts used during the research
│   ├── jm_paper_verify.py     #   core verification of Sections 3–5 claims (JM sampling)
│   ├── jm_crosscheck.py       #   unbiased cross-tabulation via Jacobson–Matthews
│   ├── jm_sharp_bound.py      #   sharpness of the n²-divisibility bound
│   ├── prove_sharp_bound.py   #   n=6 sharpness proof data
│   ├── sandpile_verify.py     #   Theorem 7.1 (sandpile bridge) verification
│   ├── sandpile_n6_exhaustive.py   #   n=6 exhaustive sandpile census (Prop 6.1)
│   ├── sandpile_cor37_proof.py     #   Corollary 3.7 structural proof
│   ├── sandpile_proof_final.py     #   sandpile closed-form proof
│   ├── sandpile_expand.py          #   sandpile expansion diagnostics
│   ├── verify_n8_jm.py        #   n=8: explicit 32|det, 64∤det witness
│   ├── verify_n16.py          #   n=16 verification
│   ├── verify_ergodicity.py   #   Jacobson–Matthews ergodicity check on n=5
│   ├── verify_hankel_proof.py #   Hankel-structure verification
│   ├── verify_sign_formula.py #   sign formula check (cyclic family)
│   ├── verify_unified_criterion.py #   unified divisibility criterion
│   ├── verify_definitive.py        #   definitive criterion check
│   ├── verify_final_clean.py       #   final cleanup run
│   ├── verify_even_n.py /_2.py     #   even-n criterion
│   ├── verify_base_column.py       #   base-column invariant check
│   ├── verify_n8*.py               #   n=8 sharpness family
│   ├── deep_investigation*.py      #   iterative exploration phase
│   ├── explore_*.py                #   kernel / mixing / structural exploration
│   ├── check_*.py                  #   Latin-square sanity checks
│   ├── red_team_step1.py /_verify.py   #   independent adversarial verification
│   ├── snf_explore.py / snf_v2_profile.py   #   Smith normal form profiling
│   ├── run_T3.py                   #   theorem T3 check runner
│   └── ...
├── results/                   # Outputs captured during research
│   ├── n6_results.txt                #   exhaustive n=6 SNF census (9 408 squares)
│   ├── phase12/, phase13/, phase13c/ #   auxiliary investigation datasets
│   ├── phase14/, phase14b/, phase14c/
│   ├── phase14*_output.txt           #   consolidated run logs
│   ├── prefold_svg/                  #   figures (SVG)
│   └── proof_n5_data.json            #   n=5 proof-assisting data
├── LICENSE                    # MIT
├── requirements.txt
└── .gitignore
```

## Dependencies

Python ≥ 3.10, `numpy ≥ 1.24`, `sympy ≥ 1.12`.

```bash
python -m venv .venv
.venv\Scripts\activate          # Windows
source .venv/bin/activate       # Linux / macOS
pip install -r requirements.txt
```

## Reproducing the paper

All scripts under `scripts/` are standalone (stdlib + numpy + sympy).
Each script contains its own inlined Jacobson–Matthews sampler or
exhaustive enumerator as appropriate.  Typical usage:

```bash
python scripts/jm_paper_verify.py           # ≈ 1 min, covers the main claims
python scripts/sandpile_n6_exhaustive.py    # ≈ 1–2 s, n=6 census
python scripts/verify_n8_jm.py              # ≈ 5–10 s, n=8 sharpness witness
python scripts/jm_sharp_bound.py            # ≈ 1–2 min, full sharpness sweep
python scripts/red_team_verify.py           # adversarial re-run
```

The empirical frequencies quoted in the paper
(e.g. “≈11% of Jacobson–Matthews samples at n=8 satisfy
32 ∣ det E_std but 64 ∤ det E_std”, “347 samples with dim ker = 1
across n=8, 12, 16”, “73% conditional frequency”, “2 700 samples
across n ≤ 20”) are direct outputs of these scripts; re-running
them under the seeds printed in the script headers reproduces the
figures up to statistical noise.

## Key paper results

- **Theorem (unified divisibility)**: `n² / gcd(n,2) | det(E_std)` for every Latin square, every n ≥ 2.
- **Theorem (F₂-rank criterion, n ≡ 2 mod 4)**: `n² | det(E_std) ⇔ rank_{F₂}(A mod 2) < n − 1`.
- **Theorem (adjugate criterion, n ≡ 0 mod 4)**: `n² | det(E_std) ⇔ adj(A mod 2) · 𝟏 = 0 over F₂`.
- **Proposition (n=6 census)**: of 9 408 reduced Latin squares, 8 832 satisfy `36 | det(E_std)` and 576 have the anomalous 18-divisible class.
- **Universal counterexamples**: an explicit infinite family of skip-one circulant parity patterns yields counterexamples to `n² | det(E_std)` for every n ≡ 2 mod 4, n ≥ 6.

See `paper/Determinant Divisibility of Centered Latin Squares.pdf` for the full statements and proofs.

## Citation

Babanskyy, O. (2026). *Determinant Divisibility of Centered Latin Squares: A Unified Theorem and the F₂-Rank Obstruction*.

## License

MIT — see [LICENSE](LICENSE).

## Contact

Issues and pull requests welcome.
