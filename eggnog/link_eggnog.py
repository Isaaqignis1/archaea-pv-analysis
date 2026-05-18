#!/usr/bin/env python3
"""join phasomeit pv members to eggnog annotations by locus tag"""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


ANNOT_COLS = [
    "COG_category", "COG_category_primary",
    "COG_function", "broad_role",
    "arCOG_id", "COG_id",
    "Description", "Preferred_name", "annotated",
]


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--members", required=True, type=Path)
    p.add_argument("--eggnog",  required=True, type=Path)
    p.add_argument("--out",     required=True, type=Path)
    return p.parse_args()


def main() -> None:
    args = parse_args()
    args.out.parent.mkdir(parents=True, exist_ok=True)

    members = pd.read_csv(args.members)
    eggnog  = pd.read_csv(args.eggnog, sep="\t")

    # one row per locus tag, eggnog can repeat
    eggnog_slim = (
        eggnog.dropna(subset=["locus_tag"])
              .drop_duplicates(subset=["locus_tag"])
              [["locus_tag", *ANNOT_COLS]]
    )

    joined = members.merge(eggnog_slim, on="locus_tag", how="left")
    joined.to_csv(args.out, index=False)

    matched   = joined["annotated"].fillna(False).astype(bool).sum()
    unmatched = len(joined) - matched
    print(f"wrote {len(joined):,} rows to {args.out}")
    print(f"  matched   {matched:,}")
    print(f"  unmatched {unmatched:,}")


if __name__ == "__main__":
    main()
