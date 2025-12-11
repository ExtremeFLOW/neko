#!/usr/bin/env python3
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

FORTRAN_EXTS = {".f90", ".F90", ".F", ".f"}

# Regex patterns
TYPE_START_RE = re.compile(r"^\s*type(?!\s*\()\b[^:]*::\s*([A-Za-z_][A-Za-z0-9_]*)", re.IGNORECASE)
TYPE_END_RE = re.compile(r"^\s*end\s*type\b", re.IGNORECASE)
COMPONENT_LINE_RE = re.compile(r"::(.*)$")
HAS_POINTER_ATTR_RE = re.compile(r"\bpointer\b", re.IGNORECASE)
HAS_ALLOCATABLE_ATTR_RE = re.compile(r"\ballocatable\b", re.IGNORECASE)
IS_PROC_BINDING_RE = re.compile(r"^\s*(procedure|generic|final)\b", re.IGNORECASE)
FREE_BINDING_LINE_RE = re.compile(r"^\s*procedure\b[^:]*::(?P<rhs>[^!]+)", re.IGNORECASE)
FREE_ANY_BINDING_RE = re.compile(r"\b(?P<name>free|free_base)\b\s*(?:=>\s*(?P<impl>[A-Za-z_][A-Za-z0-9_]*))?", re.IGNORECASE)
SUBROUTINE_START_RE = re.compile(r"^\s*subroutine\s+([A-Za-z_][A-Za-z0-9_]*)\s*\(", re.IGNORECASE)
SUBROUTINE_END_RE = re.compile(r"^\s*end\s*subroutine\b", re.IGNORECASE)
DEALLOC_RE_TMPL = r"\bdeallocate\s*\([^)]*\b{comp}\b"
NULLIFY_RE_TMPL = r"\bnullify\s*\([^)]*\b{comp}\b"
PTR_TO_NULL_RE_TMPL = r"\b{comp}\b\s*=>\s*null\s*\("

class TypeInfo:
    def __init__(self, name, start_line):
        self.name = name
        self.start_line = start_line
        self.end_line = None
        self.allocatables = []  # list of component names
        self.pointers = []      # list of component names
        self.free_map = {}  # mapping of {'free': impl_name, 'free_base': impl_name}
        self.free_declared = False

    def needs_free(self):
        return bool(self.allocatables or self.pointers)


def extract_component_names(decl_tail: str):
    """Extract component names from the tail of a declaration after '::'.
    Handles comma-separated lists with array declarators like 'a(:,:)', 'b(:)'.
    Splits only on commas that are not inside parentheses and returns the
    leading identifier of each entry.
    """
    # strip trailing comments
    decl_tail = decl_tail.split('!', 1)[0]
    items = []
    buf = []
    depth = 0
    for ch in decl_tail:
        if ch == '(':
            depth += 1
            buf.append(ch)
        elif ch == ')':
            depth = max(0, depth - 1)
            buf.append(ch)
        elif ch == ',' and depth == 0:
            item = ''.join(buf).strip()
            if item:
                items.append(item)
            buf = []
        else:
            buf.append(ch)
    tail_last = ''.join(buf).strip()
    if tail_last:
        items.append(tail_last)

    names = []
    for it in items:
        # remove any initialization part at top level
        eq_pos = it.find('=')
        if eq_pos != -1:
            it = it[:eq_pos]
        it = it.strip()
        # leading identifier
        m = re.match(r"\s*([A-Za-z_][A-Za-z0-9_]*)", it)
        if m:
            names.append(m.group(1))
    return names


def parse_types(file_text: str):
    lines = file_text.splitlines()
    types = []
    in_type = False
    current = None
    for i, line in enumerate(lines, start=1):
        if not in_type:
            m = TYPE_START_RE.match(line)
            if m:
                in_type = True
                current = TypeInfo(m.group(1), i)
                continue
        else:
            # inside type definition
            if TYPE_END_RE.match(line):
                current.end_line = i
                types.append(current)
                in_type = False
                current = None
                continue
            # procedure binding lines inside type
            if IS_PROC_BINDING_RE.match(line):
                m = FREE_BINDING_LINE_RE.match(line)
                if m:
                    rhs = m.group('rhs')
                    is_deferred = bool(re.search(r"\bdeferred\b", line, flags=re.IGNORECASE))
                    for mm in FREE_ANY_BINDING_RE.finditer(rhs):
                        nm = mm.group('name').lower()
                        impl = mm.group('impl')
                        # If this is a DEFERRED binding and no explicit impl is provided, skip recording
                        if is_deferred and not impl:
                            continue
                        current.free_map[nm] = (impl or nm)
                    if current.free_map:
                        current.free_declared = True
                continue
            # component declaration lines (with ::)
            if '::' in line and (HAS_POINTER_ATTR_RE.search(line) or HAS_ALLOCATABLE_ATTR_RE.search(line)):
                if IS_PROC_BINDING_RE.match(line):
                    continue
                m = COMPONENT_LINE_RE.search(line)
                if m:
                    names = extract_component_names(m.group(1))
                    if HAS_ALLOCATABLE_ATTR_RE.search(line):
                        current.allocatables.extend(names)
                    if HAS_POINTER_ATTR_RE.search(line):
                        current.pointers.extend(names)
    return types


def find_subroutine_body(lines, proc_name):
    start = None
    for i, line in enumerate(lines):
        m = SUBROUTINE_START_RE.match(line)
        if m and m.group(1).lower() == proc_name.lower():
            start = i
            break
    if start is None:
        return None, None
    # find end
    for j in range(start+1, len(lines)):
        if SUBROUTINE_END_RE.match(lines[j]):
            return start, j
    return start, len(lines)-1


def check_free_cleans(lines, start, end, allocs, ptrs):
    text = "\n".join(lines[start:end+1])
    missing = { 'deallocate': [], 'nullify': [] }
    for a in allocs:
        if not re.search(DEALLOC_RE_TMPL.format(comp=re.escape(a)), text, flags=re.IGNORECASE):
            missing['deallocate'].append(a)
    for p in ptrs:
        if not (re.search(NULLIFY_RE_TMPL.format(comp=re.escape(p)), text, flags=re.IGNORECASE) or
                re.search(PTR_TO_NULL_RE_TMPL.format(comp=re.escape(p)), text, flags=re.IGNORECASE)):
            missing['nullify'].append(p)
    return missing


def main():
    # Collect files
    files = []
    for p in ROOT.rglob('*'):
        if p.is_file() and p.suffix in FORTRAN_EXTS:
            files.append(p)
    files = sorted(files)

    flagged = []
    report_lines = ["# Free-check report", "", f"Scanned files: {len(files)}", ""]

    for f in files:
        try:
            txt = f.read_text(encoding='utf-8', errors='ignore')
        except Exception:
            continue
        types = parse_types(txt)
        if not types:
            continue
        lines = txt.splitlines()
        file_flagged = False
        details_added = False
        for tinfo in types:
            if not tinfo.needs_free():
                continue
            if not tinfo.free_declared:
                if not details_added:
                    report_lines.append(f"## {f.relative_to(ROOT)}")
                    details_added = True
                report_lines.append(f"- Type `{tinfo.name}`: MISSING type-bound `free` (has alloc/pointer components)")
                file_flagged = True
                continue
            # choose implementation name: prefer 'free' if present, else 'free_base'
            impl = None
            chosen = None
            # Prefer a binding that has an implementation resolved (i.e., not deferred-only)
            if 'free' in tinfo.free_map and tinfo.free_map['free']:
                chosen = 'free'
                impl = tinfo.free_map['free']
            elif 'free_base' in tinfo.free_map and tinfo.free_map['free_base']:
                chosen = 'free_base'
                impl = tinfo.free_map['free_base']
            else:
                # shouldn't happen because free_declared True implies free_map non-empty
                if not details_added:
                    report_lines.append(f"## {f.relative_to(ROOT)}")
                    details_added = True
                report_lines.append(f"- Type `{tinfo.name}`: has a procedure binding but neither 'free' nor 'free_base' was captured")
                file_flagged = True
                continue
            s, e = find_subroutine_body(lines, impl)
            if s is None:
                if not details_added:
                    report_lines.append(f"## {f.relative_to(ROOT)}")
                    details_added = True
                report_lines.append(f"- Type `{tinfo.name}`: {chosen} bound to `{impl}`, but implementation not found in file")
                file_flagged = True
                continue
            missing = check_free_cleans(lines, s, e, tinfo.allocatables, tinfo.pointers)
            if missing['deallocate'] or missing['nullify']:
                if not details_added:
                    report_lines.append(f"## {f.relative_to(ROOT)}")
                    details_added = True
                if missing['deallocate']:
                    report_lines.append(f"- Type `{tinfo.name}`: free missing deallocate for: {', '.join(missing['deallocate'])}")
                if missing['nullify']:
                    report_lines.append(f"- Type `{tinfo.name}`: free missing nullify for: {', '.join(missing['nullify'])}")
                file_flagged = True
        if file_flagged:
            flagged.append(f)

    # Write outputs
    tmp = ROOT / 'tmp_free_check'
    tmp.mkdir(parents=True, exist_ok=True)
    (tmp / 'flagged_free_types.txt').write_text("\n".join(str(p.relative_to(ROOT)) for p in flagged) + ("\n" if flagged else ""))
    (tmp / 'free_check_report.md').write_text("\n".join(report_lines) + "\n")

    print(f"Flagged files: {len(flagged)}")
    for p in flagged:
        print(f" - {p.relative_to(ROOT)}")
    print(f"Report: tmp_free_check/free_check_report.md")

if __name__ == '__main__':
    sys.exit(main())
