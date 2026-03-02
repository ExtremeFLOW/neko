# Neko documentation guide for AI agents

## Overview

Neko's documentation consists of the following:

- A repo-level README.md file.
- README.md files inside individual folders in `examples` that document the
  example provided.
- AGENTS.md files, written for LLM-based agents, but potentially useful for a
  human reader.
- A CHANGELOG.md, which accumulates changes introduced in PRs.
- Docstrings for modules, types and procedures inside the source code, written
  in Doxygen style.
- Markdown pages located under `doc/pages`, comprising the user and developer
  guides, appendices, etc.

Doxygen is used to produce HTML pages combining the automatically parsed
docstrings and all the manually written pages, see `doc/Doxyfile.in`.

## Rules

- All documentation is written in British English.
- In source code, docstrings should respect the max line length of 80 columns.
  In markdown, this should be respected whenever possible, but exceptions are
  allowed, for example for large tables.
- In markdown files, section headers should contain a tag so that they can be
  referenced, for example, `{#user-file}`. AGENTS.md files are exempt from this
  rule.
- Whenever applicable, one should use a relevant link when referencing a section
  or type. For example, "see the [examples](@ref programming-examples) section".
- The same applies for types, for example `[space_t](#space::space_t)`.
- Shell commands, paths, and similar things should be written in `monospace`.

## Reviewing documentation

If you are asked to review documentation in a PR or for new code, you should
do the following.

- Make sure the documentation follows the rules above.
- Check for typos.
- Check for grammatical or syntactical mistakes.
- Check that each new procedure, type, and module has a docstring.
  - For procedures, each dummy argument should be documented with a @param
    annotation.
  - For types, the bare minimum is a docstring above the type declaration.
    Ideally, each procedure declaration inside the type should have a docstring
    matching the first sentence of the procedure's docstring at its
    implementation. Ideally, each component of the type should be documented.
    Sometimes a group of components can be documented by one docstring when it
    makes logical sense.
  - If you feel very confident, you can try to provide a suggestion for the
    docstrings. If you are unsure, just point out that they are missing.
- Check whether the new source code introduces new parameters for the Neko case
  file, or changes the behaviour or presence of existing ones. This is typically
  done with subroutines from the `json_utils` module or directly with the
  capabilities of `json_fortran`. If so, make sure that new parameters or the
  identified changes are documented. Make sure that the documentation follows
  what is actually implemented in the source.
- Make sure that the CHANGELOG.md is relevantly updated.
- Check that the copyright statement in the header of the source files is
  correctly updated to include the year of the reviewed contribution.
