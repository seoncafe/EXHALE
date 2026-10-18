#!/usr/bin/env python3
"""Check every Fortran environment read against the survey's source excerpts.

The report is JSON so a review can inspect missing comments and mismatched
excerpts without depending on a temporary script or a previous console log.
"""

from __future__ import annotations

import argparse
import json
import re
from pathlib import Path


READ = re.compile(
    r"\bget_environment_variable\s*\(\s*&?\s*['\"]([^'\"]+)['\"]",
    re.IGNORECASE | re.DOTALL,
)
ROW = re.compile(r"^\| `([^`]+)` \|")
SOURCE_SUFFIXES = {".f90", ".F90", ".inc"}


def comment_text(line: str) -> str | None:
    stripped = line.lstrip()
    if stripped.startswith("!") and not stripped.startswith("!$omp"):
        return stripped[1:].strip()
    return None


def normalized(value: str) -> str:
    return " ".join(value.split())


def source_reads(tree: Path) -> dict[str, list[dict[str, object]]]:
    found: dict[str, list[dict[str, object]]] = {}
    for path in sorted((tree / "src").rglob("*")):
        if path.suffix not in SOURCE_SUFFIXES or not path.is_file():
            continue
        source = path.read_text(errors="replace")
        lines = source.splitlines()
        for match in READ.finditer(source):
            line_no = source.count("\n", 0, match.start()) + 1
            if "!" in lines[line_no - 1].split("get_environment_variable", 1)[0]:
                continue
            name = match.group(1)
            start = max(0, line_no - 81)
            named = [i for i in range(start, line_no - 1)
                     if (comment_text(lines[i]) or "").find(name) >= 0]
            excerpt = ""
            if named:
                i = named[-1]
                first = i
                while first > start and comment_text(lines[first - 1]) is not None:
                    first -= 1
                block = [comment_text(lines[k]) or "" for k in range(first, i + 1)]
                j = i + 1
                while j < line_no - 1 and comment_text(lines[j]) is not None:
                    block.append(comment_text(lines[j]) or "")
                    j += 1
                excerpt = normalized(" ".join(block))
            found.setdefault(name, []).append({
                "file": str(path.relative_to(tree)),
                "line": line_no,
                "named_comment": excerpt,
            })
    return found


def survey_rows(path: Path) -> dict[str, str]:
    lines = path.read_text().splitlines()
    start = next(i for i, line in enumerate(lines) if line.startswith("## 5."))
    end = next((i for i in range(start + 1, len(lines))
                if lines[i].startswith("## ") or
                lines[i].startswith("### item-named")), len(lines))
    rows: dict[str, str] = {}
    for line_no, line in enumerate(lines[start:end], start + 1):
        match = ROW.match(line)
        if not match:
            continue
        name = match.group(1)
        if name in rows:
            raise ValueError(f"duplicate survey row: {path}: {name}")
        cells = [cell.strip() for cell in line.strip().strip("|").split("|")]
        if len(cells) != 4:
            raise ValueError(f"expected four survey columns: {path}:{line_no}: {name}")
        if cells[2] not in {"yes", "no"}:
            raise ValueError(f"invalid documented field: {path}:{line_no}: {name}")
        rows[name] = cells[3]
    return rows


def audit(tree: Path, survey: Path) -> dict[str, object]:
    reads = source_reads(tree)
    rows = survey_rows(survey)
    absent_rows = sorted(set(reads) - set(rows))
    absent_reads = sorted(set(rows) - set(reads))
    missing_comments = []
    quote_mismatches = []
    for name in sorted(set(rows) & set(reads)):
        sites = reads[name]
        for site in sites:
            if not site["named_comment"]:
                missing_comments.append({"name": name, **site})
        quote = normalized(rows[name])
        if quote and not any(quote in str(site["named_comment"])
                             for site in sites):
            quote_mismatches.append({"name": name, "quote": rows[name],
                                     "sites": sites})
    return {
        "tree": str(tree), "survey": str(survey),
        "survey_rows": len(rows), "source_names": len(reads),
        "read_sites": sum(map(len, reads.values())),
        "absent_survey_rows": absent_rows,
        "absent_source_reads": absent_reads,
        "sites_without_named_comment": missing_comments,
        "quotes_not_in_named_comment": quote_mismatches,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("tree", type=Path)
    parser.add_argument("survey", type=Path)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    report = audit(args.tree.resolve(), args.survey.resolve())
    data = json.dumps(report, indent=2) + "\n"
    if args.output:
        args.output.write_text(data)
    else:
        print(data, end="")
    if any(report[key] for key in (
        "absent_survey_rows", "absent_source_reads",
        "sites_without_named_comment", "quotes_not_in_named_comment",
    )):
        raise SystemExit(1)


if __name__ == "__main__":
    main()
