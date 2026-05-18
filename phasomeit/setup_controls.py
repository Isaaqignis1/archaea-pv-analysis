#!/usr/bin/env python3
"""stage per-genus phasomeit input folders for the bacterial controls"""

from __future__ import annotations

import argparse
import shutil
from pathlib import Path

import pandas as pd


# species -> genus folder
GENUS_MAP: dict[str, str] = {
    "Campylobacter jejuni": "Campylobacter",
    "Campylobacter coli":   "Campylobacter",
    "Escherichia coli":     "Escherichia",
    "Brucella abortus":     "Brucella",
    "Brucella melitensis":  "Brucella",
}


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--prokka-dir", required=True, type=Path,
                   help="root of prokka outputs, one subdir per sample")
    p.add_argument("--mapping",    required=True, type=Path,
                   help="tsv with sample_accession and scientific_name")
    p.add_argument("--out-root",   required=True, type=Path)
    return p.parse_args()


def main() -> None:
    args = parse_args()
    if not args.prokka_dir.exists():
        raise SystemExit(f"prokka folder not found: {args.prokka_dir}")
    if not args.mapping.exists():
        raise SystemExit(f"mapping file not found: {args.mapping}")

    args.out_root.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(args.mapping, sep="\t")
    if not {"sample_accession", "scientific_name"}.issubset(df.columns):
        raise SystemExit("mapping needs sample_accession and scientific_name columns")

    df["genus_group"] = df["scientific_name"].map(GENUS_MAP)
    if df["genus_group"].isna().any():
        raise SystemExit(f"unmapped species:\n{df[df['genus_group'].isna()]}")

    total_copied = total_missing = 0

    for genus, subdf in df.groupby("genus_group"):
        genus_dir = args.out_root / genus
        genus_dir.mkdir(parents=True, exist_ok=True)

        gbk_list_path = genus_dir / f"{genus}_gbks.txt"
        missing_path  = genus_dir / "missing_gbks.txt"

        gbk_paths: list[str] = []
        missing:   list[str] = []

        for _, row in subdf.iterrows():
            sid = row["sample_accession"]
            gbk = args.prokka_dir / sid / f"{sid}.gbk"
            gbk_paths.append(str(gbk))
            if not gbk.exists():
                missing.append(str(gbk))
                continue
            shutil.copy2(gbk, genus_dir / gbk.name)
            total_copied += 1

        gbk_list_path.write_text("\n".join(gbk_paths) + "\n", encoding="utf-8")

        if missing:
            missing_path.write_text("\n".join(missing) + "\n", encoding="utf-8")
            total_missing += len(missing)
        elif missing_path.exists():
            missing_path.unlink()

        print(f"{genus}: listed={len(subdf)} missing={len(missing)} "
              f"copied={len(subdf) - len(missing)}")

    print(f"\ntotal copied: {total_copied}, missing: {total_missing}")
    print(f"output root: {args.out_root}")


if __name__ == "__main__":
    main()
