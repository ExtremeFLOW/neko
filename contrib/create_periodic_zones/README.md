## `create_periodic_zones`

`create_periodic_zones` reads a Neko `.nmsh` file, converts selected pairs of
labeled boundary zones into periodic zones, and writes a new `.nmsh`.

The utility assumes that each requested pair is related by a pure translation
and that the facets match one-to-one.

## Usage

```text
create_periodic_zones input.nmsh output.nmsh "(1,2),(3,4)" [--tol=value]
```

The zone-pair specification is flexible. The following examples are equivalent:

```text
"(1,2),(3,4)"
"1:2 3:4"
"1 2 3 4"
```

## Behaviour

- Existing periodic facets in the input mesh are preserved.
- Requested labeled zone pairs are removed from the labeled-zone table and
  replaced by periodic facets.
- Labeled zones not mentioned in the pair list are copied unchanged.
- Curved-element information is preserved.

## Recommended workflow

1. Prepare or convert the input mesh so that the candidate periodic boundaries
   are currently represented as labeled zones.
2. Run `create_periodic_zones`.
3. Run `mesh_checker` on the output mesh to verify the resulting zone data.

For a simple box mesh, `genmeshbox` may be the easier option when the desired
periodicity is known at mesh-generation time. `create_periodic_zones` is mainly
useful when periodicity needs to be added after the `.nmsh` already exists.
