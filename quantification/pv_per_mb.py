#!/usr/bin/env python3
"""join pv counts with genome lengths and emit per-genome and per-genus pv/mb"""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--pv-counts",       required=True, type=Path)
    p.add_argument("--lengths-main",    required=True, type=Path)
    p.add_argument("--lengths-control", required=True, type=Path)
    p.add_argument("--out-genome",      required=True, type=Path)
    p.add_argument("--out-genus",       required=True, type=Path)
    return p.parse_args()


def load_lengths(path: Path, source: str) -> pd.DataFrame:
    df = pd.read_csv(path)[["fasta_name", "genus", "genome_size_bp", "genome_size_mb"]].copy()
    df["source"] = source
    return df


def main() -> None:
    args = parse_args()
    args.out_genome.parent.mkdir(parents=True, exist_ok=True)
    args.out_genus.parent.mkdir(parents=True, exist_ok=True)

    pv = pd.read_csv(args.pv_counts)
    lengths = pd.concat(
        [
            load_lengths(args.lengths_main, "main"),
            load_lengths(args.lengths_control, "control"),
        ],
        ignore_index=True,
    )

    merged = pv.merge(
        lengths[["source", "fasta_name", "genome_size_bp", "genome_size_mb"]],
        on=["source", "fasta_name"],
        how="left",
    )
    if "genus_x" in merged.columns:
        merged = merged.rename(columns={"genus_x": "genus"}).drop(columns=["genus_y"])

    merged["ppv_per_mb"] = merged["PV_count"] / merged["genome_size_mb"]
    merged = merged[[
        "source", "genus", "fasta_name",
        "PV_count", "genome_size_bp", "genome_size_mb", "ppv_per_mb",
    ]].sort_values(["source", "genus", "fasta_name"])

    merged.to_csv(args.out_genome, index=False)
    print(f"per-genome: {args.out_genome} ({len(merged):,} rows)")

    genus_summary = (
        merged
        .groupby(["source", "genus"])
        .agg(
            n_genomes=("fasta_name",        "count"),
            mean_PV=("PV_count",            "mean"),
            sd_PV=("PV_count",              "std"),
            mean_genome_size_mb=("genome_size_mb", "mean"),
            mean_ppv_per_mb=("ppv_per_mb",  "mean"),
            sd_ppv_per_mb=("ppv_per_mb",    "std"),
        )
        .round(6)
        .reset_index()
        .sort_values(["source", "mean_ppv_per_mb"], ascending=[True, False])
    )
    genus_summary.to_csv(args.out_genus, index=False)
    print(f"per-genus:  {args.out_genus} ({len(genus_summary):,} rows)")


if __name__ == "__main__":
    main()
