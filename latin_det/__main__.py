"""CLI entry point for the ``latin_det`` package.

All subcommands accept ``--seed`` (for reproducible samplers) and
``--out`` (to write a CSV / JSON certificate in ``results/certified/``).
Run ``python -m latin_det --help`` for the full list.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import sys
from pathlib import Path

import numpy as np

from .core import (
    bareiss_det,
    det_A,
    det_Estd,
    difference_matrix,
    is_latin_square,
    parity_pattern,
)
from .enumerate import enumerate_reduced_latin
from .f2 import adjugate_f2, ker_f2, rank_f2
from .lifts import konig_lift
from .samplers import jacobson_matthews_samples, switch_chain_balanced
from .snf import invariant_factors, v_p


def _hash(values) -> str:
    h = hashlib.sha256()
    for v in values:
        h.update(str(v).encode("utf-8"))
        h.update(b"\n")
    return h.hexdigest()


# ---------------------------------------------------------------------------
def cmd_enum(args: argparse.Namespace) -> None:
    n = args.n
    out = Path(args.out) if args.out else None
    rows = []
    for idx, L in enumerate(enumerate_reduced_latin(n)):
        A = difference_matrix(L)
        d = bareiss_det(A)
        r = rank_f2(A)
        rows.append(
            {
                "id": idx,
                "det_A": int(d),
                "det_Estd": int(n * d),
                "rank_F2_A": int(r),
                "v2_detA": v_p(d, 2),
                "v2_detEstd": v_p(n * d, 2),
            }
        )
    if out:
        out.parent.mkdir(parents=True, exist_ok=True)
        with out.open("w", newline="", encoding="utf-8") as fh:
            writer = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
            writer.writeheader()
            writer.writerows(rows)
        print(f"Wrote {len(rows)} rows to {out}")
    else:
        counts = {}
        for row in rows:
            key = row["det_Estd"]
            counts[key] = counts.get(key, 0) + 1
        for k in sorted(counts):
            print(f"det_Estd = {k:>10}  count = {counts[k]}")
    # Consistency checks for paper claims.
    if n == 5:
        assert len(rows) == 56, f"Expected 56 reduced LS at n=5, got {len(rows)}"
    if n == 6:
        assert len(rows) == 9408, f"Expected 9408 reduced LS at n=6, got {len(rows)}"


# ---------------------------------------------------------------------------
def cmd_jm(args: argparse.Namespace) -> None:
    n = args.n
    samples = args.samples
    out = Path(args.out) if args.out else None
    generator = jacobson_matthews_samples(n, samples, seed=args.seed)
    rows = []
    for idx, L in enumerate(generator):
        assert is_latin_square(L), "JM produced non-Latin square"
        A = difference_matrix(L)
        d = bareiss_det(A)
        r = rank_f2(A)
        rows.append(
            {
                "id": idx,
                "det_A": int(d),
                "rank_F2_A": int(r),
                "v2_detA": v_p(d, 2),
            }
        )
    if out:
        out.parent.mkdir(parents=True, exist_ok=True)
        with out.open("w", newline="", encoding="utf-8") as fh:
            writer = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
            writer.writeheader()
            writer.writerows(rows)
        print(f"Wrote {len(rows)} rows to {out}")
    # Summary statistics.  Theorem 3.2 asserts n | det(E_std) = n * det(A),
    # which is trivial; the informative signal is the distribution of
    # v_2(det A) and v_p(det A) for odd primes p | n.
    v2_hist: dict = {}
    for row in rows:
        v2_hist[row["v2_detA"]] = v2_hist.get(row["v2_detA"], 0) + 1
    print(f"n={n} samples={samples} seed={args.seed}")
    print("  v_2(det A) histogram: " + ", ".join(
        f"{k}:{v}" for k, v in sorted(v2_hist.items())
    ))
    r_hist: dict = {}
    for row in rows:
        r_hist[row["rank_F2_A"]] = r_hist.get(row["rank_F2_A"], 0) + 1
    print("  rank_F2(A) histogram: " + ", ".join(
        f"{k}:{v}" for k, v in sorted(r_hist.items())
    ))


# ---------------------------------------------------------------------------
def cmd_switch(args: argparse.Namespace) -> None:
    n = args.n
    samples = args.samples
    burnin = args.burnin
    thin = args.thin
    out = Path(args.out) if args.out else None
    gen = switch_chain_balanced(n, samples, burnin=burnin, thin=thin, seed=args.seed)
    rows = []
    rank_dist = {}
    for idx, P in enumerate(gen):
        # Build B = A mod 2 from the pattern P via Lemma 4.2.
        B = np.bitwise_xor(P[:-1, :-1], P[:-1, -1:].astype(np.int8))
        r = rank_f2(B)
        rank_dist[r] = rank_dist.get(r, 0) + 1
        rows.append({"id": idx, "rank_F2_B": int(r)})
    if out:
        out.parent.mkdir(parents=True, exist_ok=True)
        with out.open("w", newline="", encoding="utf-8") as fh:
            writer = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
            writer.writeheader()
            writer.writerows(rows)
        print(f"Wrote {len(rows)} rows to {out}")
    for r in sorted(rank_dist):
        print(f"  rank_F2(B) = {r}: count = {rank_dist[r]}")


# ---------------------------------------------------------------------------
def cmd_snf(args: argparse.Namespace) -> None:
    """SNF of the cyclic Latin square for every n up to ``args.max_n``."""
    out = Path(args.out) if args.out else None
    rows = []
    for n in range(2, args.max_n + 1):
        L = np.array([[(i + j) % n + 1 for j in range(n)] for i in range(n)])
        A = difference_matrix(L)
        factors = invariant_factors(A)
        d = bareiss_det(A)
        rows.append(
            {
                "n": n,
                "det_A": int(d),
                "det_Estd": int(n * d),
                "invariant_factors": ";".join(str(x) for x in factors),
                "expected": "1, " + ", ".join([str(n)] * max(0, n - 2)),
            }
        )
    if out:
        out.parent.mkdir(parents=True, exist_ok=True)
        with out.open("w", newline="", encoding="utf-8") as fh:
            writer = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
            writer.writeheader()
            writer.writerows(rows)
        print(f"Wrote {len(rows)} rows to {out}")
    for r in rows:
        print(r)


# ---------------------------------------------------------------------------
def cmd_witness(args: argparse.Namespace) -> None:
    """Certify the explicit n=10 counterexample from Example 6.4."""
    L = np.array([
        [ 2,  5,  6,  4,  3,  8,  7,  9, 10,  1],
        [10,  6,  2,  7,  5,  1,  9,  8,  4,  3],
        [ 3,  4,  1, 10,  8,  5,  2,  7,  9,  6],
        [ 7,  1, 10,  8,  6,  9,  3,  4,  2,  5],
        [ 1,  7,  3,  6,  4, 10,  8,  2,  5,  9],
        [ 9,  2,  5,  3, 10,  6,  4,  1,  8,  7],
        [ 6,  9,  8,  1,  2,  7,  5, 10,  3,  4],
        [ 5,  8,  7,  2,  9,  4,  1,  3,  6, 10],
        [ 4,  3,  9,  5,  1,  2, 10,  6,  7,  8],
        [ 8, 10,  4,  9,  7,  3,  6,  5,  1,  2],
    ])
    assert is_latin_square(L), "Witness is not a Latin square"
    A = difference_matrix(L)
    d = bareiss_det(A)
    r2 = rank_f2(A)
    P = parity_pattern(L)
    data = {
        "source": "Example ex:n10 (Babanskyy 2026)",
        "n": 10,
        "det_A": int(d),
        "det_Estd": int(10 * d),
        "rank_F2_A": int(r2),
        "v2_det_A": v_p(d, 2),
        "v5_det_A": v_p(d, 5),
        "doubly_balanced_parity": bool(
            all(int(P[i].sum()) == 5 for i in range(10))
            and all(int(P[:, j].sum()) == 5 for j in range(10))
        ),
        "claim_50_divides_detEstd": int(10 * d) % 50 == 0,
        "claim_100_does_NOT_divide_detEstd": int(10 * d) % 100 != 0,
    }
    out = Path(args.out) if args.out else None
    if out:
        out.parent.mkdir(parents=True, exist_ok=True)
        out.write_text(json.dumps(data, indent=2), encoding="utf-8")
        print(f"Wrote certificate to {out}")
    print(json.dumps(data, indent=2))

    # Hard consistency assertions against the printed manuscript values.
    assert d == 15_427_045, f"det(A) mismatch: expected 15427045, got {d}"
    assert r2 == 9, f"rank_F2(A) mismatch: expected 9, got {r2}"
    assert int(10 * d) % 100 != 0, "Witness unexpectedly satisfies 100 | det(Estd)"


# ---------------------------------------------------------------------------
def cmd_padic(args: argparse.Namespace) -> None:
    """Scan (n, p, k_p, v_p(det A)) across sampled JM instances."""
    out = Path(args.out) if args.out else None
    rows = []
    primes = [int(p) for p in args.primes.split(",")]
    for n in range(2, args.max_n + 1):
        for L in jacobson_matthews_samples(n, args.samples_per_n, seed=args.seed):
            A = difference_matrix(L)
            d = bareiss_det(A)
            for p in primes:
                if n % p != 0:
                    continue
                k_p = _ker_dim_mod_p(A, p)
                row_bound = v_p(n, p) if n % 2 == 1 else v_p(n // 2, p) if p == 2 else v_p(n, p)
                bound = row_bound + max(0, k_p - 1)
                rows.append(
                    {
                        "n": n,
                        "p": p,
                        "k_p": int(k_p),
                        "v_p_detA": v_p(d, p),
                        "bound": bound,
                        "slack": v_p(d, p) - bound,
                    }
                )
    if out:
        out.parent.mkdir(parents=True, exist_ok=True)
        with out.open("w", newline="", encoding="utf-8") as fh:
            writer = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
            writer.writeheader()
            writer.writerows(rows)
        print(f"Wrote {len(rows)} rows to {out}")
    violations = [r for r in rows if r["slack"] < 0]
    print(f"Bound violations: {len(violations)}/{len(rows)}")


def _ker_dim_mod_p(A: np.ndarray, p: int) -> int:
    """Dimension of the kernel of ``A`` modulo the prime ``p``."""
    A = np.asarray(A, dtype=object)
    n = A.shape[0]
    M = [[int(A[i, j]) % p for j in range(n)] for i in range(n)]
    rank = 0
    col = 0
    row = 0
    while row < n and col < n:
        # Find pivot.
        pivot = None
        for i in range(row, n):
            if M[i][col] % p != 0:
                pivot = i
                break
        if pivot is None:
            col += 1
            continue
        M[row], M[pivot] = M[pivot], M[row]
        inv = pow(M[row][col], -1, p)
        for j in range(col, n):
            M[row][j] = (M[row][j] * inv) % p
        for i in range(n):
            if i != row and M[i][col] % p != 0:
                f = M[i][col]
                for j in range(col, n):
                    M[i][j] = (M[i][j] - f * M[row][j]) % p
        rank += 1
        row += 1
        col += 1
    return n - rank


# ---------------------------------------------------------------------------
def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog="latin_det")
    sub = parser.add_subparsers(dest="cmd", required=True)

    p_enum = sub.add_parser("enum", help="Exhaustively enumerate reduced Latin squares")
    p_enum.add_argument("--n", type=int, required=True)
    p_enum.add_argument("--out", type=str, default=None)
    p_enum.set_defaults(func=cmd_enum)

    p_jm = sub.add_parser("jm", help="Jacobson-Matthews samples")
    p_jm.add_argument("--n", type=int, required=True)
    p_jm.add_argument("--samples", type=int, default=1000)
    p_jm.add_argument("--seed", type=int, default=0)
    p_jm.add_argument("--out", type=str, default=None)
    p_jm.set_defaults(func=cmd_jm)

    p_sw = sub.add_parser("switch", help="Switch-chain samples on R(n, n/2)")
    p_sw.add_argument("--n", type=int, required=True)
    p_sw.add_argument("--samples", type=int, default=4000)
    p_sw.add_argument("--burnin", type=int, default=50_000)
    p_sw.add_argument("--thin", type=int, default=200)
    p_sw.add_argument("--seed", type=int, default=42)
    p_sw.add_argument("--out", type=str, default=None)
    p_sw.set_defaults(func=cmd_switch)

    p_snf = sub.add_parser("snf", help="Cyclic Latin-square SNF table")
    p_snf.add_argument("--max-n", type=int, default=12)
    p_snf.add_argument("--out", type=str, default=None)
    p_snf.set_defaults(func=cmd_snf)

    p_w = sub.add_parser("witness", help="Certify the n=10 counterexample (Ex. 6.4)")
    p_w.add_argument("--out", type=str, default=None)
    p_w.set_defaults(func=cmd_witness)

    p_p = sub.add_parser("padic", help="p-adic bound scan across JM samples")
    p_p.add_argument("--max-n", type=int, default=20)
    p_p.add_argument("--primes", type=str, default="2,3,5")
    p_p.add_argument("--samples-per-n", type=int, default=300)
    p_p.add_argument("--seed", type=int, default=0)
    p_p.add_argument("--out", type=str, default=None)
    p_p.set_defaults(func=cmd_padic)

    args = parser.parse_args(argv)
    args.func(args)
    return 0


if __name__ == "__main__":
    sys.exit(main())
