# Linting and Formatting tools

This repository contains a set of tools to help you lint and format your code.

## Linter

The [linter](lint.sh) script will help you lint your code. It uses the `flint`
tool to check your code based on the recommended style guide.

To use the linter, run the following command:

```sh
./lint.sh <path>
```

Where `<path>` is the path to the file or directory you want to lint. If you
don't provide a path, the script will lint all modified fortran files in the
repository.

The linter will provide a list of suggested improvements to your code. You can
also see the list of improvements in the report (`flinter-report.txt`) generated
by the linter.

## Formatter

The [formatter](format.sh) script will help you format your code. It uses the
`findent` tool to format your code based on the recommended style guide.

To use the formatter, run the following command:

```sh
./format.sh <path>
```

Where `<path>` is the path to the file or directory you want to format. If you
don't provide a path, the script will format all modified fortran files in the
repository.

Additional settings can be found by running `./format.sh --help`.
