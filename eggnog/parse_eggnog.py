#!/usr/bin/env python3
"""parse eggnog-mapper output into combined annotation and per-genus broad-role tables"""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


COG_DESCRIPTIONS: dict[str, str] = {
    "J": "Translation",
    "A": "RNA processing and modification",
    "K": "Transcription",
    "L": "Replication, recombination and repair",
    "B": "Chromatin structure and dynamics",
    "D": "Cell cycle control and division",
    "Y": "Nuclear structure",
    "V": "Defense mechanisms",
    "T": "Signal transduction",
    "M": "Cell wall/membrane/envelope biogenesis",
    "N": "Cell motility",
    "Z": "Cytoskeleton",
    "W": "Extracellular structures",
    "U": "Intracellular trafficking and secretion",
    "O": "Post-translational modification and chaperones",
    "X": "Mobilome: prophages and transposons",
    "C": "Energy production and conversion",
    "G": "Carbohydrate transport and metabolism",
    "E": "Amino acid transport and metabolism",
    "F": "Nucleotide transport and metabolism",
    "H": "Coenzyme transport and metabolism",
    "I": "Lipid transport and metabolism",
    "P": "Inorganic ion transport and metabolism",
    "Q": "Secondary metabolite biosynthesis and catabolism",
    "R": "General function prediction only",
    "S": "Function unknown",
}

BROAD_ROLE: dict[str, str] = {
    "C": "Core cellular Machinery",
    "H": "Core cellular Machinery",
    "I": "Core cellular Machinery",
    "F": "Core cellular Machinery",
    "J": "Core cellular Machinery",
    "L": "Core cellular Machinery",
    "O": "Core cellular Machinery",
    "D": "Core cellular Machinery",
    "U": "Core cellular Machinery",
    "P": "Adaptive Metabolism",
    "G": "Adaptive Metabolism",
    "Q": "Adaptive Metabolism",
    "E": "Adaptive Metabolism",
    "M": "Surface & Defence",
    "N": "Surface & Defence",
    "V": "Surface & Defence",
    "T": "Environmental Response",
    "K": "Environmental Response",
    "R": "Poorly Characterised",
    "S": "Poorly Characterised",
    "A": "Other",
    "B": "Other",
    "Y": "Other",
    "Z": "Other",
    "W": "Other",
    "X": "Other",
}


def parse_annotations(filepath: Path, domain: str) -> pd.DataFrame:
    print(f"parsing {domain}: {filepath}")

    rows: list[list[str]] = []
    header: list[str] | None = None

    with filepath.open() as f:
        for line in f:
            if line.startswith("##"):
                continue
            if line.startswith("#query") and header is None:
                header = line.lstrip("#").strip().split("\t")
            elif line.startswith("#"):
                continue
            else:
                rows.append(line.strip().split("\t"))

    df = pd.DataFrame(rows, columns=header)
    df.rename(columns={"#query": "query"}, inplace=True)

    df["sample_id"] = df["query"].str.extract(r"^([^_]+)")
    df["locus_tag"] = df["query"].str.extract(r"^[^_]+_(.+)$")

    if domain == "Archaea":
        df["arCOG_id"] = df["eggNOG_OGs"].str.extract(r"(arCOG\d+)@2157")
        df["COG_id"]   = df["eggNOG_OGs"].str.extract(r"(COG\d+)@1\|root")
    else:
        df["arCOG_id"] = None
        df["COG_id"]   = df["eggNOG_OGs"].str.extract(r"(COG\d+)@1\|root")

    keep = [
        "query", "sample_id", "locus_tag", "eggNOG_OGs",
        "arCOG_id", "COG_id",
        "COG_category", "Description", "Preferred_name",
        "GOs", "KEGG_ko", "PFAMs",
    ]
    df = df[[c for c in keep if c in df.columns]]
    df["domain"] = domain

    df["COG_category_primary"] = df["COG_category"].str.strip().str[0]
    df["COG_function"] = df["COG_category_primary"].map(COG_DESCRIPTIONS).fillna("Unknown")
    df["broad_role"]   = df["COG_category_primary"].map(BROAD_ROLE).fillna("Unknown")
    df["annotated"]    = df["COG_category"].notna() & (df["COG_category"].str.strip() != "-")

    return df


def load_sample_genus_map(path: Path) -> pd.DataFrame:
    sep = "\t" if path.suffix.lower() in {".tsv", ".tab"} else ","
    df = pd.read_csv(path, sep=sep)

    if "sample_id" in df.columns:
        pass
    elif "ID" in df.columns:
        df = df.rename(columns={"ID": "sample_id"})
    elif "sample_accession" in df.columns:
        df = df.rename(columns={"sample_accession": "sample_id"})
    else:
        raise SystemExit(f"{path}: need sample_id / ID / sample_accession")

    if "genus" not in df.columns:
        raise SystemExit(f"{path}: need a genus column")

    return df[["sample_id", "genus"]]


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--archaea",          required=True, type=Path)
    p.add_argument("--bacteria",         required=True, type=Path)
    p.add_argument("--outdir",           required=True, type=Path)
    p.add_argument("--sample-genus-map", type=Path, default=None,
                   help="optional sample-to-genus mapping for the per-genus view")
    return p.parse_args()


def main() -> None:
    args = parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    archaea_df  = parse_annotations(args.archaea,  "Archaea")
    bacteria_df = parse_annotations(args.bacteria, "Bacteria")
    combined    = pd.concat([archaea_df, bacteria_df], ignore_index=True)

    all_path = args.outdir / "all_annotations.tsv"
    combined.to_csv(all_path, sep="\t", index=False)
    print(f"saved {all_path}  ({len(combined):,} rows)")

    summary = (
        combined[combined["annotated"]]
        .groupby(["sample_id", "domain",
                  "COG_category_primary", "COG_function", "broad_role"])
        .size()
        .reset_index(name="gene_count")
    )
    summary_path = args.outdir / "per_sample_COG_summary.tsv"
    summary.to_csv(summary_path, sep="\t", index=False)
    print(f"saved {summary_path}  ({len(summary):,} rows)")

    archaea_out = archaea_df[[
        "query", "sample_id", "locus_tag", "arCOG_id", "COG_id",
        "COG_category", "COG_category_primary", "COG_function",
        "broad_role", "Description", "Preferred_name", "annotated",
    ]]
    archaea_out_path = args.outdir / "archaea_annotations_clean.tsv"
    archaea_out.to_csv(archaea_out_path, sep="\t", index=False)
    print(f"saved {archaea_out_path}  ({len(archaea_out):,} rows)")

    if args.sample_genus_map is not None and args.sample_genus_map.exists():
        mapping = load_sample_genus_map(args.sample_genus_map)
        joined  = combined.merge(mapping, on="sample_id", how="inner")
        joined  = joined[joined["annotated"]]

        per_genus = (
            joined
            .groupby(["domain", "genus", "broad_role"])
            .size()
            .reset_index(name="gene_count")
        )

        strain_counts = (
            mapping.groupby("genus")["sample_id"]
                   .nunique()
                   .reset_index(name="n_strains")
        )
        per_genus = per_genus.merge(strain_counts, on="genus", how="left")
        per_genus["per_strain"] = per_genus["gene_count"] / per_genus["n_strains"]
        per_genus = per_genus.sort_values(["domain", "genus", "broad_role"])

        out = args.outdir / "per_genus_broad_role_normalised.tsv"
        per_genus.to_csv(out, sep="\t", index=False)
        print(f"saved {out}  ({len(per_genus):,} rows)")
    else:
        print("skipped per-genus normalisation (pass --sample-genus-map to enable)")

    for domain, grp in combined.groupby("domain"):
        total = len(grp)
        annotated = grp["annotated"].sum()
        print(f"  {domain}: {total:,} proteins, {annotated:,} annotated "
              f"({100 * annotated / total:.1f}%)")


if __name__ == "__main__":
    main()
