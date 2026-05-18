#!/usr/bin/env python3
"""stage per-genus phasomeit input folders from <Genus>_gffs.txt lists"""

from __future__ import annotations

import argparse
import re
import shutil
from pathlib import Path


def genus_from_filename(p: Path) -> str:
    name = p.stem
    if name.endswith("_gffs"):
        return name[:-5]
    return re.sub(r"_gff(s)?$", "", name)


def gff_to_gbk(gff_path: str) -> str:
    return re.sub(r"\.gff$", ".gbk", gff_path)


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--in-dir",   required=True, type=Path,
                   help="directory with <Genus>_gffs.txt files")
    p.add_argument("--out-root", required=True, type=Path,
                   help="root for per-genus folders")
    return p.parse_args()


def main() -> None:
    args = parse_args()
    if not args.in_dir.exists():
        raise SystemExit(f"input folder not found: {args.in_dir}")
    args.out_root.mkdir(parents=True, exist_ok=True)

    txt_files = sorted(args.in_dir.glob("*_gffs.txt"))
    if not txt_files:
        raise SystemExit(f"no '*_gffs.txt' files in {args.in_dir}")

    total_copied = total_missing = 0

    for txt in txt_files:
        genus = genus_from_filename(txt)
        genus_dir = args.out_root / genus
        genus_dir.mkdir(parents=True, exist_ok=True)

        gbk_list_path = genus_dir / f"{genus}_gbks.txt"
        missing_path  = genus_dir / "missing_gbks.txt"

        gbk_paths: list[str] = []
        missing:   list[str] = []

        with txt.open() as f:
            for line in f:
                line = line.strip()
                if line:
                    gbk_paths.append(gff_to_gbk(line))

        gbk_list_path.write_text(
            "\n".join(gbk_paths) + ("\n" if gbk_paths else ""),
            encoding="utf-8",
        )

        for gbk in gbk_paths:
            src = Path(gbk)
            if not src.exists():
                missing.append(str(src))
                continue
            shutil.copy2(src, genus_dir / src.name)
            total_copied += 1

        if missing:
            missing_path.write_text("\n".join(missing) + "\n", encoding="utf-8")
            total_missing += len(missing)
        elif missing_path.exists():
            missing_path.unlink()

        print(f"{genus}: listed={len(gbk_paths)} missing={len(missing)} "
              f"copied={len(gbk_paths) - len(missing)}")

    print(f"\ntotal copied: {total_copied}, missing: {total_missing}")
    print(f"output root: {args.out_root}")


if __name__ == "__main__":
    main()
