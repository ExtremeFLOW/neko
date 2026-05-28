# icem2re2

One-shot converter from an ICEM / ANSYS Fluent `.msh` mesh to a Nek5000 `.re2`
binary mesh.

Pipeline: `.msh` → `.rea` (via `mshconvert.py`) → `.re2` (via `pymech`).

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
./icem2re2.py input.msh output.re2 --bcs bcs.json [--keep-rea]
```

- `input.msh`  — ANSYS Fluent `.msh` exported from ICEM CFD.
- `output.re2` — destination Neko `.re2` mesh.
- `--bcs`      — JSON file mapping Fluent zone ids (integers) to Nek boundary
                 condition letters.
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

### Example

```
./icem2re2.py fluent.msh fluent.re2 --bcs bcs.json
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
