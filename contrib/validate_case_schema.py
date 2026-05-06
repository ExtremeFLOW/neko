#!/usr/bin/env python3
"""Validate Neko case files against the JSON Schema."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

try:
    import json5
except ImportError as exc:
    raise SystemExit(
        "Missing dependency: json5."
    ) from exc

try:
    import jsonschema
except ImportError as exc:
    raise SystemExit(
        "Missing dependency: jsonschema."
    ) from exc

try:
    from referencing import Registry, Resource
except ImportError as exc:
    raise SystemExit(
        "Missing dependency: referencing"
    ) from exc


DEFAULT_SCHEMA = (
    Path(__file__).resolve().parents[1] / "doc" / "schemas" / "case-file.schema.json"
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Validate Neko case files using a JSON5-compatible frontend and the "
            "repository schema."
        )
    )
    parser.add_argument(
        "paths",
        nargs="+",
        help="Case files or directories containing .case/.json files.",
    )
    parser.add_argument(
        "--schema",
        default=str(DEFAULT_SCHEMA),
        help="Path to the schema file.",
    )
    return parser.parse_args()


def iter_case_files(paths: list[str]) -> list[Path]:
    files: set[Path] = set()
    for raw_path in paths:
        path = Path(raw_path)
        if path.is_dir():
            files.update(path.rglob("*.case"))
            files.update(path.rglob("*.json"))
        else:
            files.add(path)
    return sorted(files)


def load_schema(path: Path) -> jsonschema.protocols.Validator:
    """Load the root schema and preload sibling schemas for local $ref resolution."""
    with path.open(encoding="utf-8") as handle:
        schema = json.load(handle)

    # Build a local schema registry from every JSON schema file under the schema
    # root directory, including nested folders such as `simcomps/`. Each file is
    # registered under its declared `$id`, so cross-file references resolve
    # locally without fetching anything externally.
    registry = Registry()
    for schema_file in path.parent.rglob("*.json"):
        with schema_file.open(encoding="utf-8") as handle:
            candidate = json.load(handle)
        schema_id = candidate.get("$id")
        if schema_id:
            registry = registry.with_resource(
                schema_id, Resource.from_contents(candidate)
            )

    # Pick the correct validator class for the root schema's declared draft,
    # verify that the schema itself is well-formed, then create a validator
    # instance that resolves `$ref`s through the preloaded registry above.
    validator_cls = jsonschema.validators.validator_for(schema)
    validator_cls.check_schema(schema)
    return validator_cls(schema, registry=registry)


def main() -> int:
    args = parse_args()
    schema_path = Path(args.schema).resolve()
    validator = load_schema(schema_path)
    case_files = iter_case_files(args.paths)

    if not case_files:
        print("No case files found.", file=sys.stderr)
        return 1

    failures = 0
    for case_file in case_files:
        try:
            with case_file.open(encoding="utf-8") as handle:
                data = json5.load(handle)
            validator.validate(data)
            print(f"OK   {case_file}")
        except jsonschema.ValidationError as exc:
            failures += 1
            print(f"FAIL {case_file}")
            print(f"     {exc.message}")
            if exc.json_path:
                print(f"     at {exc.json_path}")
        except Exception as exc:  # noqa: BLE001
            failures += 1
            print(f"FAIL {case_file}")
            print(f"     {exc}")

    if failures:
        print(f"{failures} file(s) failed validation.", file=sys.stderr)
        return 1

    print(f"Validated {len(case_files)} file(s).")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
