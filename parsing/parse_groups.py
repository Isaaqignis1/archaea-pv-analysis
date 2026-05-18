#!/usr/bin/env python3
"""parse phasomeit group html pages into group_summary, annotated_functions, members and pairwise csvs"""

from __future__ import annotations

import argparse
import csv
import re
from html import unescape
from pathlib import Path


SKIP_FILES = {"index.html", "tracts.html", "Not in group.html"}

# cell background colour -> pv status
COLOUR_STATUS = {
    "#00b000": "intragenic_pv",
    "#ff9900": "regulatory_pv",
    "#cccccc": "non_pv_homologue",
}

SUMMARY_HEADERS = {
    "group", "name", "likelyfunction",
    "pvingene", "totalpv", "totalgenes",
}

COLOUR_RE     = re.compile(r"background-color\s*:\s*(#[0-9a-fA-F]{6})", re.I)
HREF_RE       = re.compile(r"href\s*=\s*[\"']#([^\"']+)[\"']")
PLAIN_LOCI_RE = re.compile(r"\b([A-Z]{2,}_\d{5})\b")
SUBFILE_RE    = re.compile(r"^\d+_\d+\.html$")
GROUP_RE      = re.compile(r"^(\d+)\.html$")


def unescape_text(html: str) -> str:
    html = re.sub(r"(?is)<script.*?>.*?</script>", " ", html)
    html = re.sub(r"(?is)<style.*?>.*?</style>",   " ", html)
    html = re.sub(r"(?is)<.*?>", " ", html)
    return re.sub(r"\s+", " ", unescape(html)).strip()


def norm(text: str) -> str:
    return re.sub(r"[^a-z0-9]+", "", text.lower())


def find_tables(html: str) -> list[str]:
    return re.findall(r"(?is)<table\b.*?>.*?</table>", html)


def parse_table(table_html: str) -> list[list[str]]:
    rows = []
    for row_html in re.findall(r"(?is)<tr\b.*?>.*?</tr>", table_html):
        cells = re.findall(r"(?is)<t[dh]\b.*?>(.*?)</t[dh]>", row_html)
        rows.append([unescape_text(c) for c in cells])
    return rows


def parse_table_raw(table_html: str) -> list[list[str]]:
    rows = []
    for row_html in re.findall(r"(?is)<tr\b.*?>.*?</tr>", table_html):
        cells = re.findall(r"(?is)<t[dh]\b.*?>(.*?)</t[dh]>", row_html)
        rows.append(cells)
    return rows


def safe_int(x) -> int | str:
    try:
        return int(str(x).strip())
    except Exception:
        return ""


def decode_strain_cell(raw_html: str) -> dict:
    result = {"pv_status": "absent", "pv_locus": "", "all_loci": ""}

    colour_m = COLOUR_RE.search(raw_html)
    if colour_m:
        colour = colour_m.group(1).lower()
        result["pv_status"] = COLOUR_STATUS.get(colour, f"unknown:{colour}")
        href_m = HREF_RE.search(raw_html)
        if href_m:
            result["pv_locus"] = href_m.group(1)

    all_loci = re.findall(r"href=[\"']#([A-Z]+_\d+)[\"']", raw_html)
    plain    = PLAIN_LOCI_RE.findall(unescape_text(raw_html))
    seen, combined = set(), []
    for locus in all_loci + plain:
        if locus not in seen:
            seen.add(locus)
            combined.append(locus)
    result["all_loci"] = ",".join(combined)
    return result


def parse_summary_table(table_html: str, genus: str, file: str,
                        page_title: str) -> list[dict]:
    raw_rows  = parse_table_raw(table_html)
    text_rows = parse_table(table_html)
    if len(text_rows) < 2:
        return []

    scalar_map: dict[str, int] = {}
    strain_cols: list[tuple[int, str]] = []
    for i, h in enumerate(text_rows[0]):
        nh = norm(h)
        if   nh == "group":          scalar_map["group_num"]       = i
        elif nh == "name":           scalar_map["group_name"]      = i
        elif nh == "likelyfunction": scalar_map["likely_function"] = i
        elif nh == "pvingene":       scalar_map["pv_in_gene"]      = i
        elif nh == "totalpv":        scalar_map["total_pv"]        = i
        elif nh == "totalgenes":     scalar_map["total_genes"]     = i
        elif i >= 6:
            strain_cols.append((i, h.strip()))

    if len(scalar_map) < 5:
        return []

    out_rows: list[dict] = []
    for raw_row, text_row in zip(raw_rows[1:], text_rows[1:]):
        def get(field: str) -> str:
            idx = scalar_map.get(field)
            if idx is None or idx >= len(text_row):
                return ""
            return text_row[idx]

        pv_in_gene  = safe_int(get("pv_in_gene"))
        total_pv    = safe_int(get("total_pv"))
        total_genes = safe_int(get("total_genes"))

        regulatory_pv: int | str = ""
        if isinstance(pv_in_gene, int) and isinstance(total_pv, int):
            regulatory_pv = total_pv - pv_in_gene

        row = {
            "genus":             genus,
            "file":              file,
            "page_title":        page_title,
            "group_num":         get("group_num"),
            "group_name":        get("group_name"),
            "likely_function":   get("likely_function"),
            "pv_in_gene":        pv_in_gene,
            "total_pv":          total_pv,
            "total_genes":       total_genes,
            "regulatory_pv":     regulatory_pv,
            "has_intragenic_pv": 1 if isinstance(pv_in_gene, int) and pv_in_gene > 0 else 0,
            "has_any_pv":        1 if isinstance(total_pv, int) and total_pv > 0 else 0,
        }

        for col_idx, strain_label in strain_cols:
            raw_cell = raw_row[col_idx] if col_idx < len(raw_row) else ""
            decoded  = decode_strain_cell(raw_cell)
            safe_label = re.sub(r"[^a-z0-9]+", "_", strain_label.lower()).strip("_")
            row[f"strain_{safe_label}_status"]   = decoded["pv_status"]
            row[f"strain_{safe_label}_pv_locus"] = decoded["pv_locus"]
            row[f"strain_{safe_label}_all_loci"] = decoded["all_loci"]

        out_rows.append(row)

    return out_rows


def parse_annotated_functions(table_html: str, genus: str, file: str,
                              page_title: str) -> list[dict]:
    # phasomeit spells occuring with one r
    text_rows = parse_table(table_html)
    if len(text_rows) < 2:
        return []

    headers = [norm(h) for h in text_rows[0]]

    def col(names: list[str]) -> int | None:
        for n in names:
            if n in headers:
                return headers.index(n)
        return None

    func_col  = col(["function", "likelyfunction", "annotatedfunctions"])
    count_col = col(["occuring", "occurring", "count", "n"])
    score_col = col(["score"])

    if func_col is None:
        return []

    out = []
    for row in text_rows[1:]:
        if not any(row):
            continue
        out.append({
            "genus":              genus,
            "file":               file,
            "page_title":         page_title,
            "annotated_function": row[func_col]  if func_col  is not None and func_col  < len(row) else "",
            "occuring":           row[count_col] if count_col is not None and count_col < len(row) else "",
            "score":              row[score_col] if score_col is not None and score_col < len(row) else "",
        })
    return out


def is_annotated_functions_table(text_rows: list[list[str]]) -> bool:
    if not text_rows:
        return False
    headers = {norm(h) for h in text_rows[0]}
    return (
        ("function" in headers or "annotatedfunctions" in headers)
        and ("occuring" in headers or "occurring" in headers or "count" in headers)
    )


MEMBER_BLOCK_RE = re.compile(
    r'<p\s+id="([^"]+)">\s*<b>[^:]+:\s*[^<]+</b>\s*'
    r'<a\s+href="([^"]+)">\[tract entry\]</a>.*?'
    r'Function:\s*(.*?)</div>'
    r'.*?white-space\s*:\s*pre[^>]*>(.*?)</div>',
    re.DOTALL,
)


def extract_tract_info(seq_block_html: str) -> dict:
    tract_spans = re.findall(
        r"<span[^>]*background-color\s*:\s*#0033dd[^>]*>(.*?)</span>",
        seq_block_html, re.I,
    )
    tract_residues = "".join(tract_spans)

    tract_aa_pos: int | str = ""
    for line in seq_block_html.split("<br />"):
        if "#0033dd" not in line.lower():
            continue
        offset_m = re.search(r"\+(\d+)\s+", line)
        if offset_m:
            line_offset = int(offset_m.group(1))
            seq_start   = re.search(r"\+\d+\s+", line)
            seq_part    = line[seq_start.end():] if seq_start else line
            before      = re.split(r"<span", seq_part)[0]
            before_plain = re.sub(r"<[^>]+>", "", before)
            tract_aa_pos = line_offset + len(before_plain) + 1
        break

    return {
        "tract_residues":  tract_residues,
        "tract_length_aa": len(tract_residues),
        "tract_aa_pos":    tract_aa_pos,
    }


def parse_members(html: str, genus: str, file: str,
                  page_title: str, group_num: str) -> list[dict]:
    out = []
    for m in MEMBER_BLOCK_RE.finditer(html):
        locus_id        = m.group(1)
        tract_entry_url = m.group(2)
        function        = re.sub(r"\s+", " ", m.group(3)).strip()
        seq_block_html  = m.group(4)

        tract_info      = extract_tract_info(seq_block_html)
        strain_file_m   = re.search(r"/strains/(\d+)\.html", tract_entry_url)
        strain_file_idx = strain_file_m.group(1) if strain_file_m else ""

        out.append({
            "genus":           genus,
            "file":            file,
            "page_title":      page_title,
            "group_num":       group_num,
            "locus_tag":       locus_id,
            "strain_file_idx": strain_file_idx,
            "function":        function,
            "tract_residues":  tract_info["tract_residues"],
            "tract_length_aa": tract_info["tract_length_aa"],
            "tract_aa_pos":    tract_info["tract_aa_pos"],
        })
    return out


PAIRWISE_RE = re.compile(
    r'(\w+)\s+vs:\s+(\w+)\s+in\s+<a\s+href="([^"]+)">(.*?)</a>[^<]*<br/>'
    r'\s*<span[^>]*>Gene length:\s*(\d+)bp\s*/\s*(\d+)aa\s+PV:\s*(Yes|No)</span>'
    r'\s*<br\s*/>\s*<span[^>]*>Function:\s*(.*?)</span>'
    r'.*?Score:\s*([\d.]+)\s+bits:\s*([\d.]+)\s+e-value:\s*(\S+)<br\s*/>'
    r'\s*length:\s*(\d+)\s+gaps:\s*(\d+)\s+id:\s*(\d+)\s+'
    r'positives:\s*(\d+)\s+coverage:\s*([\d.]+)\s+query coverage\s+([\d.]+)',
    re.DOTALL,
)


def parse_pairwise_file(html: str, genus: str, file: str,
                        page_title: str, group_num: str) -> list[dict]:
    out = []
    for m in PAIRWISE_RE.finditer(html):
        (query, target, _strain_href, strain_name,
         bp, aa, pv_status, function,
         score, bits, evalue,
         length, gaps, identity, positives, coverage, query_cov) = m.groups()

        out.append({
            "genus":           genus,
            "file":            file,
            "page_title":      page_title,
            "group_num":       group_num,
            "query_locus":     query,
            "target_locus":    target,
            "strain_name":     strain_name.strip(),
            "gene_length_bp":  safe_int(bp),
            "gene_length_aa":  safe_int(aa),
            "pv_status":       pv_status,
            "function":        function.strip(),
            "blast_score":     score,
            "blast_bits":      bits,
            "blast_evalue":    evalue,
            "align_length":    safe_int(length),
            "align_gaps":      safe_int(gaps),
            "align_identity":  safe_int(identity),
            "align_positives": safe_int(positives),
            "align_coverage":  coverage,
            "align_query_cov": query_cov,
        })
    return out


def group_num_from_filename(name: str) -> str:
    m = re.match(r"^(\d+)", name)
    return m.group(1) if m else ""


def write_csv(path: Path, fieldnames: list[str], rows: list[dict]) -> None:
    if not rows:
        with path.open("w", newline="", encoding="utf-8") as fh:
            csv.DictWriter(fh, fieldnames=fieldnames,
                           extrasaction="ignore").writeheader()
        return

    # let dynamic strain_* columns through
    all_keys = list(dict.fromkeys(
        list(fieldnames) +
        [k for row in rows for k in row if k not in fieldnames]
    ))
    with path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=all_keys,
                                extrasaction="ignore", restval="")
        writer.writeheader()
        writer.writerows(rows)


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--base",   required=True, type=Path,
                   help="phasomeit genus root, one subdir per genus")
    p.add_argument("--outdir", required=True, type=Path)
    return p.parse_args()


def main() -> None:
    args = parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    group_summary_rows: list[dict]      = []
    annotated_function_rows: list[dict] = []
    member_rows: list[dict]             = []
    pairwise_rows: list[dict]           = []

    for genus_dir in sorted(p for p in args.base.iterdir() if p.is_dir()):
        groups_dir = genus_dir / "summary_tracts" / "groups"
        if not groups_dir.is_dir():
            continue

        for html_file in sorted(groups_dir.glob("*.html")):
            fname = html_file.name
            if fname in SKIP_FILES:
                continue
            try:
                html = html_file.read_text(errors="ignore")
            except Exception:
                continue

            genus      = genus_dir.name
            group_num  = group_num_from_filename(fname)
            title_m    = re.search(r"(?is)<title>(.*?)</title>", html)
            page_title = unescape_text(title_m.group(1)) if title_m else ""

            if SUBFILE_RE.match(fname):
                pairwise_rows.extend(
                    parse_pairwise_file(html, genus, fname, page_title, group_num)
                )
                continue

            if not GROUP_RE.match(fname):
                continue

            for table_html in find_tables(html):
                text_rows = parse_table(table_html)
                if not text_rows:
                    continue
                headers = {norm(h) for h in text_rows[0]}
                if len(headers & SUMMARY_HEADERS) >= 5:
                    group_summary_rows.extend(
                        parse_summary_table(table_html, genus, fname, page_title)
                    )
                elif is_annotated_functions_table(text_rows):
                    annotated_function_rows.extend(
                        parse_annotated_functions(table_html, genus, fname, page_title)
                    )

            member_rows.extend(
                parse_members(html, genus, fname, page_title, group_num)
            )

    write_csv(
        args.outdir / "phasomeit_group_summary.csv",
        ["genus", "file", "page_title",
         "group_num", "group_name", "likely_function",
         "pv_in_gene", "total_pv", "total_genes",
         "regulatory_pv", "has_intragenic_pv", "has_any_pv"],
        group_summary_rows,
    )
    write_csv(
        args.outdir / "phasomeit_annotated_functions.csv",
        ["genus", "file", "page_title",
         "annotated_function", "occuring", "score"],
        annotated_function_rows,
    )
    write_csv(
        args.outdir / "phasomeit_members.csv",
        ["genus", "file", "page_title", "group_num",
         "locus_tag", "strain_file_idx", "function",
         "tract_residues", "tract_length_aa", "tract_aa_pos"],
        member_rows,
    )
    write_csv(
        args.outdir / "phasomeit_pairwise.csv",
        ["genus", "file", "page_title", "group_num",
         "query_locus", "target_locus", "strain_name",
         "gene_length_bp", "gene_length_aa", "pv_status", "function",
         "blast_score", "blast_bits", "blast_evalue",
         "align_length", "align_gaps", "align_identity",
         "align_positives", "align_coverage", "align_query_cov"],
        pairwise_rows,
    )

    print(f"output: {args.outdir}")
    print(f"  group_summary       {len(group_summary_rows):,}")
    print(f"  annotated_functions {len(annotated_function_rows):,}")
    print(f"  member_rows         {len(member_rows):,}")
    print(f"  pairwise_rows       {len(pairwise_rows):,}")


if __name__ == "__main__":
    main()
