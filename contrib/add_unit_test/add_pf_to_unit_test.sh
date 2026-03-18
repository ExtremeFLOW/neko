#!/bin/sh

set -eu

# Wrapper for adding a new .pf file to an existing unit-test suite.
# <suite_name> and <pf_name> must match ^[a-z][a-z0-9_]*$.
# A leading test_ and trailing .pf in <pf_name> are accepted and stripped.
# The normalized <pf_name> must be the bare stem.

usage() {
    echo "Usage: contrib/add_unit_test/add_file_to_unit_test.sh <suite_name> <pf_name>" >&2
}

if [ "$#" -ne 2 ]; then
    usage
    exit 1
fi

script_dir=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
suite_name=$1
pf_name=$2

case "$pf_name" in
    *.pf) pf_name=${pf_name%.pf} ;;
esac

case "$pf_name" in
    test_*) pf_name=${pf_name#test_} ;;
esac

exec python3 "${script_dir}/add_unit_test.py" add-file "$suite_name" "$pf_name"
