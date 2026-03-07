#!/usr/bin/env bash

set -euo pipefail

usage() {
  cat <<'EOF'
Usage: contrib/add_unit_test/add_unit_test.sh <test_name>

Creates a new serial pFUnit test from tests/unit/templates/serial and wires it
into tests/unit/Makefile.am, configure.ac, and tests/unit/.gitignore.

The test name must match: ^[a-z][a-z0-9_]*$
EOF
}

die() {
  printf 'error: %s\n' "$1" >&2
  exit 1
}

insert_before_pattern() {
  local file="$1"
  local pattern="$2"
  local new_line="$3"
  local tmp

  if grep -Fqx "$new_line" "$file"; then
    return 0
  fi

  tmp="$(mktemp "${TMPDIR:-/tmp}/add-unit-test.XXXXXX")"
  if ! awk -v pattern="$pattern" -v new_line="$new_line" '
    index($0, pattern) && !inserted {
      print new_line
      inserted = 1
    }
    { print }
    END {
      if (!inserted) {
        exit 2
      }
    }
  ' "$file" >"$tmp"; then
    rm -f "$tmp"
    die "could not update $file"
  fi

  mv "$tmp" "$file"
}

rename_template_files() {
  local test_dir="$1"
  local test_name="$2"
  local pf_file="test_${test_name}.pf"
  local binary_name="${test_name}_test"

  mv "$test_dir/test_serial.pf" "$test_dir/$pf_file"

  perl -0pi -e \
    "s/serial_test/${binary_name}/g; s/test_serial\\.pf/${pf_file}/g" \
    "$test_dir/Makefile.in"

  perl -0pi -e \
    "s/module test_serial/module test_${test_name}/g; \
     s/end module test_serial/end module test_${test_name}/g; \
     s/test_serial_template_passes/test_template_passes/g" \
    "$test_dir/$pf_file"
}

if [[ $# -ne 1 ]]; then
  usage
  exit 1
fi

test_name="$1"
if ! [[ "$test_name" =~ ^[a-z][a-z0-9_]*$ ]]; then
  usage
  die "test name must match ^[a-z][a-z0-9_]*$"
fi

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repo_root="$(cd "$script_dir/../.." && pwd)"
unit_dir="$repo_root/tests/unit"
template_dir="$unit_dir/templates/serial"
test_dir="$unit_dir/$test_name"

[[ -d "$template_dir" ]] || die "missing template directory: $template_dir"
[[ -e "$test_dir" ]] && die "target already exists: $test_dir"

cp -R "$template_dir" "$test_dir"

find "$test_dir" -maxdepth 1 -type f \
  \( -name 'Makefile' -o -name '*.o' -o -name '*.mod' -o -name '*.a' \
     -o -name '*.inc' -o -name '*.F90' -o -name '*.log' -o -name '*.trs' \
     -o -name 'serial_test' \) -delete

rename_template_files "$test_dir" "$test_name"

insert_before_pattern "$unit_dir/Makefile.am" \
  "templates/serial" \
  "	  ${test_name}\\"
insert_before_pattern "$unit_dir/Makefile.am" \
  "templates/serial/serial_test" \
  "	${test_name}/${test_name}_test\\"
insert_before_pattern "$unit_dir/Makefile.am" \
  "templates/serial/test_serial.pf" \
  "	${test_name}/test_${test_name}.pf\\"

insert_before_pattern "$repo_root/configure.ac" \
  "tests/unit/templates/serial/Makefile" \
  "        tests/unit/${test_name}/Makefile\\"

if ! grep -Fqx "${test_name}/${test_name}_test" "$unit_dir/.gitignore"; then
  printf '%s\n' "${test_name}/${test_name}_test" >>"$unit_dir/.gitignore"
fi

cat <<EOF
Created tests/unit/${test_name}
EOF
