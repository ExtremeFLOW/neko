#!/bin/sh

set -eu

# Wrapper for creating a new unit-test suite.
# <suite_name> must match ^[a-z][a-z0-9_]*$.
# <is_parallel> must be one of: true, false, yes, no, 1, 0.

usage() {
    echo "Usage: contrib/add_unit_test/add_unit_test.sh <suite_name> <is_parallel>" >&2
}

if [ "$#" -ne 2 ]; then
    usage
    exit 1
fi

script_dir=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)

# add_unit_test.py needs Python 3.7 or newer.
if ! python3 -c 'import sys; sys.exit(0 if sys.version_info >= (3, 7) else 1)' 2>/dev/null; then
    echo "error: python3 must be version 3.7 or newer" >&2
    exit 1
fi

exec python3 "${script_dir}/add_unit_test.py" create-suite "$@"
