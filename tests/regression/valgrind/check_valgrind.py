#!/usr/bin/env python3
"""Validate valgrind leak summary against regression thresholds."""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

SUMMARY_PATTERNS = {
    "definitely":
    re.compile(
        r"definitely lost:\s*([0-9,]+)\s+bytes\s+in\s+([0-9,]+)\s+blocks"),
    "indirectly":
    re.compile(
        r"indirectly lost:\s*([0-9,]+)\s+bytes\s+in\s+([0-9,]+)\s+blocks"),
    "possibly":
    re.compile(
        r"possibly lost:\s*([0-9,]+)\s+bytes\s+in\s+([0-9,]+)\s+blocks"),
    "reachable":
    re.compile(
        r"still reachable:\s*([0-9,]+)\s+bytes\s+in\s+([0-9,]+)\s+blocks"),
    "suppressed":
    re.compile(r"suppressed:\s*([0-9,]+)\s+bytes\s+in\s+([0-9,]+)\s+blocks"),
}


def _to_int(value: str) -> int:
    return int(value.replace(",", ""))


def parse_summary(text: str) -> dict[str, tuple[int, int]]:
    if "All heap blocks were freed -- no leaks are possible" in text:
        return {
            "definitely": (0, 0),
            "indirectly": (0, 0),
            "possibly": (0, 0),
            "reachable": (0, 0),
            "suppressed": (0, 0),
        }

    parsed: dict[str, tuple[int, int]] = {}
    for key, pattern in SUMMARY_PATTERNS.items():
        match = pattern.search(text)
        if match is None:
            raise ValueError(
                f"Could not find valgrind summary line for '{key}'.")
        parsed[key] = (_to_int(match.group(1)), _to_int(match.group(2)))

    return parsed


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--log",
                        required=True,
                        help="Path to valgrind log file")
    parser.add_argument(
        "--max-definite-bytes",
        type=int,
        default=0,
        help="Maximum allowed definitely-lost bytes",
    )
    parser.add_argument(
        "--max-indirect-bytes",
        type=int,
        default=0,
        help="Maximum allowed indirectly-lost bytes",
    )
    parser.add_argument(
        "--max-possible-bytes",
        type=int,
        default=0,
        help="Maximum allowed possibly-lost bytes",
    )
    parser.add_argument(
        "--max-reachable-bytes",
        type=int,
        default=1300000,
        help="Maximum allowed still-reachable bytes",
    )
    args = parser.parse_args()

    log_path = Path(args.log)
    if not log_path.is_file():
        print(f"Valgrind log not found: {log_path}", file=sys.stderr)
        return 2

    content = log_path.read_text(encoding="utf-8", errors="replace")

    try:
        summary = parse_summary(content)
    except ValueError as err:
        print(str(err), file=sys.stderr)
        return 2

    definitely_bytes, definitely_blocks = summary["definitely"]
    indirectly_bytes, indirectly_blocks = summary["indirectly"]
    possibly_bytes, possibly_blocks = summary["possibly"]
    reachable_bytes, reachable_blocks = summary["reachable"]
    suppressed_bytes, suppressed_blocks = summary["suppressed"]

    print("Valgrind summary:")
    print(
        f"  definitely lost: {definitely_bytes} bytes in {definitely_blocks} blocks"
    )
    print(
        f"  indirectly lost: {indirectly_bytes} bytes in {indirectly_blocks} blocks"
    )
    print(
        f"  possibly lost:   {possibly_bytes} bytes in {possibly_blocks} blocks"
    )
    print(
        f"  still reachable: {reachable_bytes} bytes in {reachable_blocks} blocks"
    )
    print(
        f"  suppressed:      {suppressed_bytes} bytes in {suppressed_blocks} blocks"
    )

    errors: list[str] = []
    if definitely_bytes > args.max_definite_bytes:
        errors.append("definitely-lost bytes exceeded threshold: "
                      f"{definitely_bytes} > {args.max_definite_bytes}")

    if indirectly_bytes > args.max_indirect_bytes:
        errors.append("indirectly-lost bytes exceeded threshold: "
                      f"{indirectly_bytes} > {args.max_indirect_bytes}")

    if possibly_bytes > args.max_possible_bytes:
        errors.append("possibly-lost bytes exceeded threshold: "
                      f"{possibly_bytes} > {args.max_possible_bytes}")

    if reachable_bytes > args.max_reachable_bytes:
        errors.append("still-reachable bytes exceeded threshold: "
                      f"{reachable_bytes} > {args.max_reachable_bytes}")

    if errors:
        print("Valgrind regression check failed:", file=sys.stderr)
        for msg in errors:
            print(f"  - {msg}", file=sys.stderr)
        return 1

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
