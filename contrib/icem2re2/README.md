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

Keys are Fluent zone ids (as they appear in the `.msh` file); values are the
single-letter Nek BC codes (`v`, `W`, `o`, `P`, …).

```json
{
    "13": "v",
    "14": "W",
    "15": "W",
    "16": "o",
    "17": "W",
    "18": "W"
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
