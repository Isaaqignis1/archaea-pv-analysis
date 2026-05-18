#!/usr/bin/env python3
"""pv count per fasta across the archaeal main and bacterial control sets"""

from __future__ import annotations

import argparse
import re
from pathlib import Path

import pandas as pd


def load_mapping(genus: str, gbk_dir: Path) -> dict[str, str]:
    gbk_file = gbk_dir / genus / f"{genus}_gbks.txt"
    if not gbk_file.exists():
        raise FileNotFoundError(f"missing mapping file: {gbk_file}")

    mapping: dict[str, str] = {}
    with gbk_file.open() as f:
        for i, line in enumerate(f, start=1):
            path = line.strip()
            if not path:
                continue
            mapping[f"strain_{i}"] = Path(path).stem
    return mapping


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--main-base",    required=True, type=Path)
    p.add_argument("--main-csv",     required=True, type=Path)
    p.add_argument("--control-base", required=True, type=Path)
    p.add_argument("--control-csv",  required=True, type=Path)
    p.add_argument("--out",          required=True, type=Path)
    return p.parse_args()


def main() -> None:
    args = parse_args()
    args.out.parent.mkdir(parents=True, exist_ok=True)

    datasets = [
        {"source": "main",    "gbk_dir": args.main_base,    "csv": args.main_csv},
        {"source": "control", "gbk_dir": args.control_base, "csv": args.control_csv},
    ]

    all_results: list[dict] = []
    for ds in datasets:
        if not ds["csv"].exists():
            print(f"skipping {ds['source']}: missing {ds['csv']}")
            continue

        df = pd.read_csv(ds["csv"])
        pv_cols = [c for c in df.columns if c.endswith("_pv_locus")]

        for genus, subdf in df.groupby("genus"):
            mapping = load_mapping(genus, ds["gbk_dir"])

            for col in pv_cols:
                m = re.search(r"strain_(\d+)_pv_locus$", col)
                if not m:
                    continue
                strain_key = f"strain_{m.group(1)}"
                if strain_key not in mapping:
                    continue

                all_results.append({
                    "source":     ds["source"],
                    "genus":      genus,
                    "fasta_name": mapping[strain_key],
                    "PV_count":   int(subdf[col].notna().sum()),
                })

    out_df = pd.DataFrame(all_results)
    out_df.to_csv(args.out, index=False)
    print(f"wrote {len(out_df):,} rows to {args.out}")


if __name__ == "__main__":
    main()
