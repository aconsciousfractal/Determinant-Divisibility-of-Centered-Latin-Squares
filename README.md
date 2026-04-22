# Determinant Divisibility of Centered Latin Squares — Reproducibility Repository

Companion code, certified datasets, and verification scripts for the
paper

> **Determinant Divisibility of Centered Latin Squares**,
> Oleksiy Babanskyy, 2026.

The paper studies the integer determinant of the centered Latin-square
matrix $E = L - \tfrac{n+1}{2}J$ and the "difference" matrix
$A_{ij} = L_{ij} - L_{i,n-1}$ attached to an $n \times n$ Latin
square~$L$ with entries in $\{1,\dots,n\}$.  It proves that
$n \mid \det(E_{\mathrm{std}})$ for every~$L$, classifies the stronger
divisibility patterns $n^k$ in terms of the $\mathbb{F}_2$-rank of
$A \bmod 2$, and exhibits an explicit order-10 Latin square with
$\det(E_{\mathrm{std}}) / n^2 \notin \mathbb{Z}$.

This repository contains everything needed to re-derive the empirical
content of Sections 5–7 from scratch.

## Repository layout

```
.
├── latin_det/              # Python package (core algorithms)
│   ├── core.py             #   exact linear algebra (Bareiss, Gram, parity)
│   ├── f2.py               #   F_2 linear algebra (rank, kernel, adjugate)
│   ├── snf.py              #   Smith normal form (sympy wrapper)
│   ├── samplers.py         #   random Latin squares, switch chain
│   ├── lifts.py            #   Konig lift + MRV backtracking
│   ├── enumerate.py        #   reduced-LS enumeration for n <= 6
│   └── __main__.py         #   CLI entry point
├── scripts/
│   └── verify_isotopy_destruction.py
├── results/certified/      # Datasets referenced in the paper
│   ├── n5_det_enum.csv
│   ├── n6_det_enum.csv
│   ├── n8_jm_seed0.csv
│   ├── n10_mcmc_4000.csv
│   ├── cyclic_snf.csv
│   ├── padic_bounds_table.csv
│   ├── n10_witness.json
│   ├── isotopy_destruction.json
│   └── SHA256SUMS
├── requirements.txt
└── LICENSE                 # MIT
```

## Installation

Requires Python >= 3.10, `numpy >= 1.24`, `sympy >= 1.12`.

```bash
python -m venv .venv
.venv\Scripts\activate          # Windows
source .venv/bin/activate       # Linux / macOS
pip install -r requirements.txt
```

## Reproducing the paper's tables

From the repository root (with `PYTHONPATH=.` on Windows or `pip
install -e .` elsewhere):

| Paper reference | Command | Runtime |
|---|---|---|
| Appendix A — n=5 census | `python -m latin_det enum --n 5 --out results/certified/n5_det_enum.csv` | <1 s |
| Proposition 6.1 — n=6 census (9 408 squares) | `python -m latin_det enum --n 6 --out results/certified/n6_det_enum.csv` | ≈1.2 s |
| Corollary 3.5 — cyclic SNF, n <= 12 | `python -m latin_det snf --max-n 12 --out results/certified/cyclic_snf.csv` | <1 s |
| Example 6.4 — n=10 witness certificate | `python -m latin_det witness --out results/certified/n10_witness.json` | <1 s |
| Section 7 — random Latin squares n=8 (10 000 samples) | `python -m latin_det jm --n 8 --samples 10000 --seed 0 --out results/certified/n8_jm_seed0.csv` | ≈25 s |
| Section 7 — switch-chain R(10,5) (4 000 samples) | `python -m latin_det switch --n 10 --samples 4000 --burnin 50000 --thin 200 --seed 42 --out results/certified/n10_mcmc_4000.csv` | ≈8 s |
| Theorem 5.X — p-adic bound scan | `python -m latin_det padic --max-n 12 --primes 2,3,5 --samples-per-n 200 --seed 0 --out results/certified/padic_bounds_table.csv` | ≈3 s |
| Section 7.1 — isotopy destruction | `python scripts/verify_isotopy_destruction.py` | ≈1 s |

All deterministic outputs are bit-reproducible; their SHA-256 digests
are recorded in `results/certified/SHA256SUMS`.

## Verifying the integrity of the certified datasets

```bash
cd results/certified
sha256sum -c SHA256SUMS        # Linux / macOS
Get-Content SHA256SUMS | ForEach-Object { $h,$f = $_ -split '  '; if ((Get-FileHash -Algorithm SHA256 $f).Hash.ToLower() -ne $h) { Write-Error "MISMATCH: $f" } else { "OK: $f" } }   # Windows
```

## Key claims that can be independently verified

1. **Theorem 3.2** (trivial divisibility): every row of `n8_jm_seed0.csv`
   and `n10_mcmc_4000.csv` has `det_Estd = n · det_A` exactly.
2. **Proposition 3.8** (cyclic formula): `cyclic_snf.csv` reports
   `|det A| = n^{n-2}` and invariant factors `(1, n, n, …, n)`.
3. **Proposition 6.1** (n=6 census): `n6_det_enum.csv` has 9 408 rows;
   among them 8 832 satisfy `36 | det_Estd` and 576 are the "anomalous"
   18-divisible class.
4. **Example 6.4** (n=10 witness): `n10_witness.json` certifies
   `det_Estd = 154 270 450`, `v_2(det A) = 0`, `v_5(det A) = 1`,
   `rank_F_2(A) = 9`, and `100 ∤ det_Estd`.
5. **Section 7.1** (isotopy non-invariance): see
   `isotopy_destruction.json` — standardising the witness raises
   `v_2(det A)` from 0 to 4; only 337 / 2 000 random symbol
   relabellings preserve the odd-determinant property.

## Licence

MIT — see [LICENSE](LICENSE).  Please cite the paper when using this
code or data.

## Contact

Oleksiy Babanskyy — `aconsciousfractal@github`.  Issues and pull
requests welcome.
