#!/usr/bin/env bash
set -euo pipefail

script_dir="$(cd "$(dirname "$0")" && pwd)"
case_dir="$(cd "$script_dir/.." && pwd)"
repo_dir="$(cd "$script_dir/../../.." && pwd)"

cd "$script_dir"
gmsh mach3_cylinder.geo -2
printf '2\nmach3_cylinder\n0\n' | "$repo_dir/contrib/gmsh2nek/gmsh2nek"
rea2nbin mach3_cylinder.re2 "$case_dir/mach3_cylinder.nmsh"
