# vsearch shell completion

Completion scripts for bash, zsh and fish, all generated from one YAML
spec so that adding or changing an option is a one-place change.

Installation instructions for users are in the main
[`README.md`](../README.md); this file is about maintaining the scripts.

## Layout

```
spec/vsearch.yaml               the source of truth: every option, once
scripts/spec.py                 typed loader (Spec, Option, Arg)
scripts/gen_bash.py             -> completions/vsearch
scripts/gen_zsh.py              -> completions/_vsearch
scripts/gen_fish.py             -> completions/vsearch.fish
scripts/extract_cli_options.py  reads the option tables out of src/cli.cc
scripts/check_option_drift.py   diffs those tables against the spec
scripts/check_completions.sh    drives every check ('make check')
completions/                    generated, and tracked in git
```

Each generated file already carries the name it is installed under, so
`make install` needs no rename: bash looks up a file named after the
command, zsh an `_`-prefixed autoloaded function, fish a
`<command>.fish`.

## The spec

- `enums` — reusable value sets (`mask_method`, `strand`).
- `file_types` — extension filters for input filenames, including the
  compressed forms wherever the underlying format supports them.
  `.faa` is deliberately absent: vsearch is nucleotide-only.
- `commands` — every option vsearch treats as an action, `--help` and
  `--version` included, since the C++ parser gives those command rows
  of their own.
- `options` — every option declared once, with its argument type, and
  the `command:` and `hidden:` flags.
- `command_options` — the exhaustive command-to-options mapping.

## Changing an option

1. **Edit `spec/vsearch.yaml`. Never edit `completions/`** — those are
   build products and any hand edit is lost on the next regeneration
   (and reported by `make check` before that).
2. A new action flag goes in `commands:` *as well as* `options:`, with
   `command: true`.
3. Add the option to every relevant entry in `command_options:`. The
   mapping is **exhaustive**: nothing is implied. If `--threads` is
   valid for a command, it must be listed there.
4. `make -C completion regenerate-completions`
5. `make -C completion check`
6. Verify by hand in at least bash and fish (recipes below).

## Checks

```sh
make -C completion check          # or plain 'make check' from the top
```

Three things, in order:

1. **Drift** — `spec/vsearch.yaml` against the `option_specs` and
   `valid_options` tables in `src/cli.cc`, which are the authority.
   This is the check that matters: a spec maintained by hand goes stale,
   and a stale spec does not merely under-suggest, it offers options the
   binary rejects.
2. **Staleness** — the tracked scripts against what the generators
   currently emit.
3. **Syntax** — each script against the shell it targets.

A missing dependency (Python, PyYAML, zsh, fish) is skipped with a
notice, so a contributor without fish installed is not blocked. `make
check STRICT=1` turns every skip into a failure; CI sets it, because a
skip there would quietly turn the whole check into a pass.

What the drift check cannot see, because `cli.cc` does not carry it:
descriptions, argument types beyond value-or-flag, the enum value lists,
and the `hidden:` set. Those stay human-maintained, and the manual pages
under [`../man`](../man) are their reference. **Do not use `vsearch
--help` as the reference**: it is a summary and it omits options
(the spec was originally bootstrapped from it and had to be corrected
against the C++ tables afterwards).

## `hidden: true`

vsearch accepts a few options that completion should not offer:
`--fulldp`, `--cons_truncate`, `--slots`, `--pattern`, `--band`,
`--hspw`, `--minhsp` and `--xdrop_nw`. They still appear in
`command_options` wherever the parser accepts them — the generators skip
them when emitting suggestions but still recognise them as valid. This
is an editorial list; the drift check only verifies that a hidden option
still exists in `cli.cc`.

## File-type filters

An input `FILENAME` argument sets `arg.ftype` to a key in `file_types:`
to filter by extension. Output filenames use `type: outfile` and never
filter — the user is naming a new path. Inputs with no useful extension
convention (`--labels`, `--label_words`) use `type: file` with no
`ftype`, which accepts any file.

## Verifying by hand

```sh
# bash: source the script and drive _vsearch through COMP_WORDS
bash -c 'source completions/vsearch
  COMP_WORDS=(vsearch --cluster_fast q.fa --i); COMP_CWORD=3
  _vsearch; printf "%s\n" "${COMPREPLY[@]}"'

# fish: use the built-in driver
fish -c 'source completions/vsearch.fish
  complete -C "vsearch --cluster_fast q.fa --i"'
```

zsh cannot be driven this way: `_arguments` only runs inside the real
completion system, so `zsh -n` (part of `make check`) plus an
interactive session is as far as it goes.

## Generator notes

### bash

One `_vsearch` function: scan `${COMP_WORDS[@]}` for a command flag,
then `case "$prev"` to complete a value (enums via `compgen -W`, files
via `_vsearch_filedir`, nothing for numeric and string options), and
otherwise `case "$cmd"` to offer that command's option list.
`_vsearch_filedir` uses bash-completion's `_filedir` when it is
available and falls back to `compgen` when it is not.

A bare `<Tab>` suggests options even without a leading dash. That is a
deliberate divergence from the older bash convention, because zsh and
fish both behave that way.

The script stays compatible with bash 3.2, which is what macOS ships:
no `declare -A`, `mapfile`, `readarray`, `coproc` or `${var,,}`, and
`compopt` is guarded.

### zsh

An `_arguments` dispatcher, one branch per command, with native value
descriptors: enums become `:label:(a b c)`, input files
`:file:_files -g 'glob'`, output files `:file:_files`, everything else a
bare `:label:`.

Escaping happens in two passes and they must not be combined:
`arg_desc()` escapes the characters that have meaning inside the
`[...]` description bracket (`[`, `]`, `:`, `\`), and `sh_squote()`
then wraps the finished spec in single quotes. The description and the
action are concatenated *first*, then quoted.

### fish

Declarative: one `complete -c vsearch` directive per (condition,
option) pair, with three helpers to find which command is on the line —
`__vsearch_current_command`, `__vsearch_no_command` and
`__vsearch_command_is`. Per option: `-f` for a flag, `-x -a '...'` for
an enum, `-rfa '(__fish_complete_suffix ...)'` for a filtered input
file, `-r` for any input file, `-rF` for an output file, `-x` for a
value with no suggestions.

## History

The scripts were developed in `torognes/vsearch_completion` and moved
here so the spec can track `src/cli.cc` commit by commit; that
repository is archived. Requested in
[issue #417](https://github.com/torognes/vsearch/issues/417).
