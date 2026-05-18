#!/usr/bin/env python3
"""per-genome size and contig count from fasta files under a genus tree"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


def fasta_length(fasta_path: Path) -> tuple[int, int]:
    total_bp = 0
    seq_count = 0
    with fasta_path.open() as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                seq_count += 1
            else:
                total_bp += len(line)
    return total_bp, seq_count


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--base", required=True, type=Path)
    p.add_argument("--out",  required=True, type=Path)
    return p.parse_args()


def main() -> None:
    args = parse_args()
    args.out.parent.mkdir(parents=True, exist_ok=True)

    rows = []
    for genus_dir in sorted(args.base.iterdir()):
        if not genus_dir.is_dir():
            continue
        genus = genus_dir.name
        for fasta_path in sorted(genus_dir.glob("*.fasta")):
            total_bp, seq_count = fasta_length(fasta_path)
            rows.append({
                "fasta_name":     fasta_path.stem,
                "genus":          genus,
                "genome_size_bp": total_bp,
                "genome_size_mb": total_bp / 1_000_000,
                "n_contigs":      seq_count,
                "fasta_path":     str(fasta_path),
            })

    with args.out.open("w", newline="") as fh:
        writer = csv.DictWriter(
            fh,
            fieldnames=["fasta_name", "genus", "genome_size_bp",
                        "genome_size_mb", "n_contigs", "fasta_path"],
        )
        writer.writeheader()
        writer.writerows(rows)

    print(f"wrote {len(rows):,} rows to {args.out}")


if __name__ == "__main__":
    main()
