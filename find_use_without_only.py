import os
import re
import argparse

def find_use_without_only(directory, extensions=('.f90', '.F90', '.f95', '.f')):
    use_line = re.compile(r'^\s*use\s+\w+', re.IGNORECASE)
    has_only = re.compile(r',\s*only\s*:', re.IGNORECASE)

    results = []

    for root, _, files in os.walk(directory):
        for fname in files:
            if fname.endswith(extensions):
                path = os.path.join(root, fname)
                with open(path, 'r', encoding='utf-8', errors='ignore') as f:
                    for lineno, line in enumerate(f, 1):
                        if use_line.search(line) and not has_only.search(line):
                            results.append((path, lineno, line.strip()))
    return results

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Detect `use` statements in Fortran without `only`.")
    parser.add_argument("directory", help="Directory to scan.")
    args = parser.parse_args()

    for path, lineno, line in find_use_without_only(args.directory):
        print(f"{path}:{lineno}: {line}")

