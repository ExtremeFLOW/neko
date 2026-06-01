# icem2re2

One-shot converter from an ICEM / ANSYS Fluent `.msh` mesh to a Nek5000 `.re2`
binary mesh.

Pipeline: `.msh` → `.rea` (via `mshconvert.py`) → `.re2` (via `pymech`).

## Attribution & license

- `mshconvert.py` is **derived from** [`mikaem/tools`](https://github.com/mikaem/tools)
  by Mikael Mortensen. The upstream repository carries **no explicit license
  notice**; it is included here on the assumption that its use in Neko is
  permitted by the original author. Local changes: a Python-3 compatibility
  port, `eval()` → `float()` in node-coordinate parsing (security fix), and BC
  plumbing through `convert()` / `scan_fluent_mesh()`. If Neko's contribution
  policy requires a formal license, please confirm with the upstream author and
  add an appropriate SPDX identifier.
- `rea2re2.py` and `icem2re2.py` are original to this contribution.
- `pymech` is a third-party dependency (BSD-3-Clause), installed via `pip` — it
  is not vendored here.

## Requirements

Python 3 with:

- `numpy`
- `scipy`
- `pymech`

Install with `pip`:

```
pip install numpy scipy pymech
```

## Usage

```
./icem2re2.py input.msh output.re2 --bcs bcs.json \
    [--periodic periodic.json] [--keep-rea]
```

- `input.msh`  — ANSYS Fluent `.msh` exported from ICEM CFD.
- `output.re2` — destination Neko `.re2` mesh.
- `--bcs`      — JSON file mapping Fluent zone ids (integers) to Nek boundary
                 condition letters.
- `--periodic` — JSON file declaring translational-periodic zone pairs.
                 See "Periodic boundaries" below.
- `--keep-rea` — keep the intermediate `.rea` next to `output.re2`.

### `bcs.json` format

Keys are Fluent zone ids; values are single-character Nek BC codes
(`v`, `W`, `o`, `P`, …).

Zone ids in Fluent `.msh` files are written in **hexadecimal**. You may give
the JSON key in either form — both are accepted:

- Hexadecimal, exactly as it appears in the `.msh` (`"0xd"`, `"0x10"`), or
- Decimal, after converting from hex (`"13"`, `"16"`).

For example, a zone written as `(13` (hex) in the `.msh` file is zone id
`19` in decimal — write either `"0x13"` or `"19"` in the JSON.

```json
{
    "0xd":  "v",
    "0xe":  "W",
    "0xf":  "W",
    "0x10": "o",
    "0x11": "W",
    "0x12": "W"
}
```

### Periodic boundaries

Neko encodes periodicity in the mesh itself (the case file says nothing about
it), so the `.re2` must record, for every periodic face, the partner element
and face it pairs with. `mshconvert` handles this natively when you pass it the
zone pairs and their translation vector.

Workflow:

1. In `bcs.json`, **omit** the periodic zones. They are not regular BCs;
   mshconvert auto-assigns them BC letter `'P'` with the right partner info.
2. In `periodic.json`, declare the pairs by **zone id** and the translation
   vector from the first zone to the second. mshconvert then walks node-by-node
   on each side under that translation to find partners.

Example (zone 13 ↔ zone 14 periodic in *z* with spacing 0.05, plus walls):

```json
// bcs.json — periodic zones 13 and 14 are NOT listed
{
    "12": "v",
    "15": "o",
    "16": "W",
    "17": "W"
}
```

```json
// periodic.json — translation goes from zone 13 to zone 14
[
    {"zones": [13, 14], "displacement": [0.0, 0.0, 0.05]}
]
```

The two zones must have the same number of nodes; if mshconvert can't find a
partner for every node under the given translation it'll report it.

Only translational periodicity is supported — rotational is not implemented.

### Example

```
./icem2re2.py fluent.msh fluent.re2 --bcs bcs.json
./icem2re2.py fluent.msh fluent.re2 --bcs bcs.json --periodic periodic.json
```

## Installing on your PATH

Symlink the wrapper under the short name

```
ln -s $(pwd)/icem2re2.py ~/.local/bin/icem2re2
```

after which plain `icem2re2 …` works anywhere.

## Files

- `icem2re2.py`     — CLI wrapper, single entry point.
- `mshconvert.py`   — Fluent `.msh` → Nek `.rea` library.
- `rea2re2.py`      — Nek `.rea` → `.re2` converter (also usable standalone).
