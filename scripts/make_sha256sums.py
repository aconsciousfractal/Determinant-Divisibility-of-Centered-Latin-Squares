"""make_sha256sums.py — regenerate results/certified/SHA256SUMS, the
SHA-256 manifest promised in the paper's Code-and-data-availability
paragraph (Sec. 8). Added 2026-07-09 (see CLAIM_LEDGER.md GAP-1).

Pins the canonical (git-stored, LF-normalized) bytes of every
deterministic artefact of the certified package. Standard `sha256sum`
format with repo-root-relative paths; verify with
`sha256sum -c results/certified/SHA256SUMS` from the repository root on
an LF checkout (or re-run this script and diff — it normalizes CRLF, so
it is checkout-independent).
"""
import hashlib
import os

REPO_ROOT = os.path.normpath(os.path.join(os.path.dirname(__file__), ".."))

ARTIFACTS = [
    "results/proof_n5_data.json",
    "results/n6_results.txt",
    "results/phase12/summary.json",
    "results/phase12/full_rank_corpus.json",
    "results/phase12/deficient_controls.json",
    "results/phase12/features.csv",
    "results/certified/cyclic_snf_table.json",
    "results/certified/ex_n10_witness.json",
    "results/certified/padic_scan_output.txt",
    "results/certified/jm_n8_corpus_summary.json",
]


def main():
    lines = []
    for rel in ARTIFACTS:
        path = os.path.join(REPO_ROOT, rel)
        with open(path, "rb") as f:
            data = f.read()
        # Normalize CRLF -> LF before hashing: git stores these text
        # artifacts with LF line endings, but a Windows checkout with
        # core.autocrlf=true materializes CRLF on disk. Hashing raw disk
        # bytes there pins hashes no fresh (LF) clone can reproduce —
        # which is exactly what happened to the 2026-07-09 manifest
        # (repinned 2026-07-18; see CLAIM_LEDGER.md Errata).
        data = data.replace(b"\r\n", b"\n")
        lines.append(f"{hashlib.sha256(data).hexdigest()}  {rel}")
        print(lines[-1])
    out = os.path.join(REPO_ROOT, "results", "certified", "SHA256SUMS")
    with open(out, "w", newline="\n") as f:
        f.write("\n".join(lines) + "\n")
    print(f"\nManifest written: {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
