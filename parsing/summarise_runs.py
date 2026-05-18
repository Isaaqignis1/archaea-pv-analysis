#!/usr/bin/env python3
"""per-genus phasomeit run completeness summary"""

from __future__ import annotations

import argparse
import csv
import re
from pathlib import Path


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--root", required=True, type=Path)
    p.add_argument("--out",  required=True, type=Path)
    return p.parse_args()


def main() -> None:
    args = parse_args()
    rows = []
    for genus_dir in sorted(p for p in args.root.iterdir() if p.is_dir()):
        genus = genus_dir.name
        n_genomes = len(list(genus_dir.glob("*.gbk")))

        summary_dir = genus_dir / "summary_tracts"
        groups_dir  = summary_dir / "groups"
        strains_dir = summary_dir / "strains"

        run_success    = summary_dir.exists()
        n_strain_pages = len(list(strains_dir.glob("*.html"))) if strains_dir.exists() else 0

        n_group_html_files = 0
        n_unique_groups    = 0
        if groups_dir.exists():
            ignore = {"index.html", "tracts.html", "Not in group.html"}
            group_files = [f.name for f in groups_dir.glob("*.html") if f.name not in ignore]
            n_group_html_files = len(group_files)

            unique_groups: set[str] = set()
            for name in group_files:
                group_id = name[:-5].split("_")[0]
                if re.fullmatch(r"\d+", group_id):
                    unique_groups.add(group_id)
            n_unique_groups = len(unique_groups)

        rows.append({
            "genus":              genus,
            "run_success":        run_success,
            "n_genomes_gbk":      n_genomes,
            "n_strain_pages":     n_strain_pages,
            "n_group_html_files": n_group_html_files,
            "n_unique_groups":    n_unique_groups,
        })

    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w", newline="") as fh:
        writer = csv.DictWriter(
            fh,
            fieldnames=["genus", "run_success", "n_genomes_gbk",
                        "n_strain_pages", "n_group_html_files", "n_unique_groups"],
        )
        writer.writeheader()
        writer.writerows(rows)

    print(f"wrote {args.out}")


if __name__ == "__main__":
    main()
