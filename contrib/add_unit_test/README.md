# add_unit_test

Helpers for creating and extending pFUnit suites under `tests/unit`.

## Create a New Suite

Serial:

```bash
contrib/add_unit_test/add_unit_test.sh my_test false
```

MPI:

```bash
contrib/add_unit_test/add_unit_test.sh my_parallel_test true
```

This creates a new suite from the generic templates and wires it into:
- `tests/unit/Makefile.am`
- `configure.ac`
- `tests/unit/.gitignore`

## Add a `.pf` File to an Existing Suite

```bash
contrib/add_unit_test/add_file_to_unit_test.sh htable htable_cptr
```

This creates:

```text
tests/unit/htable/test_htable_cptr.pf
```

and adds it to:
- `tests/unit/htable/Makefile.in`
- `tests/unit/Makefile.am`