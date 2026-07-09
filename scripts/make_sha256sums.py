"""make_sha256sums.py — regenerate results/certified/SHA256SUMS, the
SHA-256 manifest promised in the paper's Code-and-data-availability
paragraph (Sec. 8). Added 2026-07-09 (see CLAIM_LEDGER.md GAP-1).

Pins the byte-exact contents of every deterministic artefact of the
certified package. Standard `sha256sum` format with repo-root-relative
paths; verify with `sha256sum -c results/certified/SHA256SUMS` from the
repository root (or re-run this script and diff).
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
        h = hashlib.sha256()
        with open(path, "rb") as f:
            for chunk in iter(lambda: f.read(1 << 20), b""):
                h.update(chunk)
        lines.append(f"{h.hexdigest()}  {rel}")
        print(lines[-1])
    out = os.path.join(REPO_ROOT, "results", "certified", "SHA256SUMS")
    with open(out, "w", newline="\n") as f:
        f.write("\n".join(lines) + "\n")
    print(f"\nManifest written: {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
