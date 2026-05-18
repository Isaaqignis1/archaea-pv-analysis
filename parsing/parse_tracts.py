#!/usr/bin/env python3
"""extract ssr tract rows from phasomeit html, combining archaeal and bacterial sets"""

from __future__ import annotations

import argparse
import csv
import re
import statistics
from collections import Counter
from html import unescape as html_unescape
from pathlib import Path


def strip_tags(html: str) -> str:
    html = re.sub(r"(?is)<script.*?>.*?</script>", " ", html)
    html = re.sub(r"(?is)<style.*?>.*?</style>",   " ", html)
    html = re.sub(r"(?is)<[^>]+>", " ", html)
    return re.sub(r"\s+", " ", html_unescape(html)).strip()


def safe_int(val):
    try:
        return int(str(val).strip())
    except (ValueError, TypeError):
        return None


def norm(s: str) -> str:
    return re.sub(r"[^a-z0-9]", "", s.lower())


def extract_strain_accession(html: str) -> str:
    # filepath line is the most reliable source of the accession
    fp = re.search(r"Filepath:\s*([^\s<]+)", html)
    if fp:
        stem = Path(fp.group(1).strip()).stem
        if stem and stem.lower() != "strain":
            return stem
    sn = re.search(r"Strain name:\s*([^<\n\r]+)", html)
    if sn:
        name = sn.group(1).strip()
        if name and name.lower() != "strain":
            return name
    return ""


def parse_tract_type(raw: str):
    m = re.match(r"^([A-Za-z]+)(\d+)$", raw.strip())
    if m:
        return m.group(1).upper(), int(m.group(2))
    return raw.strip(), None


def pv_class_from_offset(offset):
    if offset is None: return "unknown"
    if offset == 0:    return "intragenic"
    if offset < 0:     return "upstream"
    return "downstream"


def group_num_from_cell(cell_html: str):
    m = re.search(r"groups/(\d+)\.html", cell_html)
    if m:
        return int(m.group(1))
    return safe_int(strip_tags(cell_html))


MEMBER_BLOCK_RE = re.compile(
    r'<p\s+id="([^"]+)">'
    r'.*?'
    r'white-space\s*:\s*pre[^>]*>'
    r'(.*?)'
    r'</div>',
    re.DOTALL | re.IGNORECASE,
)
LINE_OFFSET_RE = re.compile(r"^\+(\d+)\s+")


def count_residues_in_seq_block(seq_block_html: str) -> int:
    total = 0
    for line in seq_block_html.split("<br />"):
        plain = re.sub(r"<[^>]+>", "", line)
        plain = html_unescape(plain)
        plain = LINE_OFFSET_RE.sub("", plain)
        total += len(plain.replace(" ", "").replace("\n", "").replace("\r", ""))
    return total


def extract_tract_aa_pos(seq_block_html: str):
    gene_length = count_residues_in_seq_block(seq_block_html)
    gene_length = gene_length if gene_length > 0 else None

    tract_aa_pos = None
    for line in seq_block_html.split("<br />"):
        if "#0033dd" not in line.lower():
            continue
        plain_line = html_unescape(re.sub(r"<[^>]+>", "", line))
        offset_m = LINE_OFFSET_RE.match(plain_line)
        if not offset_m:
            offset_m = LINE_OFFSET_RE.match(plain_line.lstrip())
        line_offset = int(offset_m.group(1)) if offset_m else 0

        before_span  = re.split(r"<span[^>]*#0033dd", line, maxsplit=1)[0]
        before_plain = html_unescape(re.sub(r"<[^>]+>", "", before_span))
        before_plain = LINE_OFFSET_RE.sub("", before_plain)
        residues_before = len(before_plain.replace(" ", ""))

        tract_aa_pos = line_offset + residues_before + 1
        break

    return tract_aa_pos, gene_length


def parse_group_html(html: str) -> dict:
    result = {}
    for m in MEMBER_BLOCK_RE.finditer(html):
        locus_tag      = m.group(1)
        seq_block_html = m.group(2)
        pos, length    = extract_tract_aa_pos(seq_block_html)
        result[locus_tag] = {
            "tract_aa_pos":       pos,
            "gene_length_aa_seq": length,
        }
    return result


def build_group_lookup(genus_dir: Path) -> dict:
    groups_dir = genus_dir / "summary_tracts" / "groups"
    if not groups_dir.is_dir():
        return {}

    lookup: dict[int, dict] = {}
    group_re = re.compile(r"^(\d+)\.html$")

    for html_file in sorted(groups_dir.glob("*.html")):
        m = group_re.match(html_file.name)
        if not m:
            continue
        group_num = int(m.group(1))

        try:
            html = html_file.read_text(errors="ignore")
        except Exception:
            continue

        members = parse_group_html(html)
        if not members:
            continue

        positions = [v["tract_aa_pos"]       for v in members.values() if v["tract_aa_pos"]       is not None]
        lengths   = [v["gene_length_aa_seq"] for v in members.values() if v["gene_length_aa_seq"] is not None]

        lookup[group_num] = {
            "median_tract_aa_pos":       statistics.median(positions) if positions else None,
            "median_gene_length_aa_seq": statistics.median(lengths)   if lengths   else None,
            "members": members,
        }

    return lookup


STRAIN_COL_NAMES = {
    "tractno":        "tract_no",
    "contig":         "contig",
    "location":       "location_bp",
    "tract":          "tract_type_raw",
    "onlength":       "on_length_raw",
    "genegroup":      "group_num_raw",
    "offsetfromgene": "offset_raw",
    "gene":           "gene",
    "length":         "gene_length_aa",
    "function":       "function",
}
STRAIN_REQUIRED = {"tractno", "tract", "location", "genegroup"}


def parse_strain_html(html: str, domain: str, genus: str,
                      fname: str, group_lookup: dict) -> list[dict]:
    rows_out: list[dict] = []
    accession = extract_strain_accession(html)
    tables = re.findall(r"(?is)<table\b.*?>.*?</table>", html)

    for tbl in tables:
        raw_rows = re.findall(r"(?is)<tr\b.*?>.*?</tr>", tbl)
        if len(raw_rows) < 2:
            continue

        header_cells = re.findall(r"(?is)<t[dh]\b.*?>(.*?)</t[dh]>", raw_rows[0])
        header_norm  = [norm(strip_tags(c)) for c in header_cells]

        if not STRAIN_REQUIRED.issubset(set(header_norm)):
            continue

        col_map = {
            STRAIN_COL_NAMES[nh]: i
            for i, nh in enumerate(header_norm) if nh in STRAIN_COL_NAMES
        }

        for row_html in raw_rows[1:]:
            cells_raw  = re.findall(r"(?is)<t[dh]\b.*?>(.*?)</t[dh]>", row_html)
            cells_text = [strip_tags(c) for c in cells_raw]

            def get_text(field: str) -> str:
                idx = col_map.get(field)
                if idx is None or idx >= len(cells_text):
                    return ""
                return cells_text[idx].strip()

            def get_raw(field: str) -> str:
                idx = col_map.get(field)
                if idx is None or idx >= len(cells_raw):
                    return ""
                return cells_raw[idx]

            tract_raw        = get_text("tract_type_raw")
            tract_unit, rlen = parse_tract_type(tract_raw)
            on_str           = get_text("on_length_raw")
            on_length        = safe_int(on_str) if on_str.lower() != "none" else None
            group_num        = group_num_from_cell(get_raw("group_num_raw"))
            offset           = safe_int(get_text("offset_raw"))
            pv_cls           = pv_class_from_offset(offset)
            gene_len_aa      = safe_int(get_text("gene_length_aa"))
            locus_tag        = get_text("gene")

            tract_aa_pos = gene_length_aa_seq = None
            if group_num is not None and group_num in group_lookup:
                grp = group_lookup[group_num]
                if locus_tag in grp["members"]:
                    tract_aa_pos       = grp["members"][locus_tag]["tract_aa_pos"]
                    gene_length_aa_seq = grp["members"][locus_tag]["gene_length_aa_seq"]
                else:
                    tract_aa_pos       = grp["median_tract_aa_pos"]
                    gene_length_aa_seq = grp["median_gene_length_aa_seq"]

            rows_out.append({
                "domain":             domain,
                "genus":              genus,
                "strain_file":        fname,
                "strain_accession":   accession,
                "tract_no":           safe_int(get_text("tract_no")),
                "contig":             get_text("contig"),
                "location_bp":        safe_int(get_text("location_bp")),
                "tract_type":         tract_raw,
                "tract_unit":         tract_unit,
                "tract_repeat_len":   rlen,
                "on_length":          on_length,
                "group_num":          group_num,
                "offset_from_gene":   offset,
                "pv_class":           pv_cls,
                "gene":               locus_tag,
                "gene_length_aa":     gene_len_aa,
                "function":           get_text("function"),
                "tract_aa_pos":       tract_aa_pos,
                "gene_length_aa_seq": gene_length_aa_seq,
            })

        if rows_out:
            break  # one valid table per strain file

    return rows_out


def crawl_base(base: Path, domain: str, all_rows: list) -> tuple[int, int]:
    if not base.exists():
        print(f"  base not found: {base}")
        return 0, 0

    n_genera = n_files = 0
    for genus_dir in sorted(p for p in base.iterdir() if p.is_dir()):
        strains_dir = genus_dir / "summary_tracts" / "strains"
        if not strains_dir.is_dir():
            continue

        genus = genus_dir.name
        print(f"  {genus} ...", end=" ", flush=True)
        group_lookup = build_group_lookup(genus_dir)

        strain_files = sorted(strains_dir.glob("*.html"))
        genus_rows = 0
        for html_file in strain_files:
            try:
                html = html_file.read_text(errors="ignore")
            except Exception as exc:
                print(f"\n    could not read {html_file}: {exc}")
                continue
            rows = parse_strain_html(html, domain, genus, html_file.name, group_lookup)
            all_rows.extend(rows)
            genus_rows += len(rows)
            n_files += 1

        n_genera += 1
        print(f"{len(strain_files)} strain files, {genus_rows:,} tracts, "
              f"{len(group_lookup)} groups indexed")

    return n_genera, n_files


FIELDNAMES = [
    "domain", "genus", "strain_file", "strain_accession",
    "tract_no", "contig", "location_bp",
    "tract_type", "tract_unit", "tract_repeat_len", "on_length",
    "group_num", "offset_from_gene", "pv_class",
    "gene", "gene_length_aa", "function",
    "tract_aa_pos", "gene_length_aa_seq",
]


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--archaea-base",  required=True, type=Path)
    p.add_argument("--controls-base", required=True, type=Path)
    p.add_argument("--outdir",        required=True, type=Path)
    return p.parse_args()


def main() -> None:
    args = parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    all_rows: list[dict] = []

    print("archaea")
    arch_genera, arch_files = crawl_base(args.archaea_base, "Archaea", all_rows)

    print("\nbacteria")
    ctrl_genera, ctrl_files = crawl_base(args.controls_base, "Bacteria", all_rows)

    out_path = args.outdir / "phasomeit_tract_data.csv"
    with out_path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=FIELDNAMES, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(all_rows)

    print(f"\noutput: {out_path}")
    print(f"total rows: {len(all_rows):,}")
    print(f"  archaea  {sum(1 for r in all_rows if r['domain']=='Archaea'):,}  "
          f"({arch_genera} genera, {arch_files} strain files)")
    print(f"  bacteria {sum(1 for r in all_rows if r['domain']=='Bacteria'):,}  "
          f"({ctrl_genera} genera, {ctrl_files} strain files)")

    for k, v in sorted(Counter(r["pv_class"] for r in all_rows).items()):
        with_pos = sum(
            1 for r in all_rows
            if r["pv_class"] == k and r["tract_aa_pos"] is not None
        )
        print(f"  {k:15s} {v:6,}  (tract_aa_pos filled: {with_pos:,})")


if __name__ == "__main__":
    main()
