"""Generate the bash completion script for vsearch.

Writes to stdout. Invoked from completion/Makefile.am:
    cd scripts && python3 gen_bash.py > ../completions/vsearch
"""

from __future__ import annotations

import sys
from collections import defaultdict

import spec as spec_mod


HEADER = """\
# vsearch completion for bash
# Auto-generated from spec/vsearch.yaml — do not edit by hand.
#
# Install: place this file (named "vsearch") in a bash-completion lookup
# directory, e.g. /usr/share/bash-completion/completions/vsearch, or source
# it from ~/.bashrc. 'make install' does this for you.

shopt -s extglob 2>/dev/null

# Complete a filename argument, optionally restricted to an extension
# pattern (e.g. "fasta|fa|fna").  Uses _filedir from bash-completion when
# available; falls back to a compgen-based implementation otherwise.
_vsearch_filedir() {
    local exts="$1"
    if declare -F _filedir >/dev/null 2>&1; then
        if [[ -n "$exts" ]]; then
            _filedir "@($exts)"
        else
            _filedir
        fi
        return
    fi
    COMPREPLY=()
    if [[ -n "$exts" ]]; then
        COMPREPLY=( $(compgen -f -X "!*.@($exts)" -- "$cur") )
        COMPREPLY+=( $(compgen -d -S / -- "$cur") )
    else
        COMPREPLY=( $(compgen -f -- "$cur") )
    fi
    compopt -o filenames 2>/dev/null
}
"""

FOOTER = """\

complete -F _vsearch vsearch
"""


def shell_quote(s: str) -> str:
    """Quote a string for use inside bash double quotes."""
    return s.replace("\\", "\\\\").replace('"', '\\"').replace("`", "\\`").replace("$", "\\$")


def fmt_pipe(words: list[str]) -> str:
    return "|".join(words)


def main() -> None:
    s = spec_mod.load()

    # Group value-taking options by what their argument completion looks like.
    enum_groups: dict[str, list[str]] = defaultdict(list)   # enum_key -> [opt_name]
    file_groups: dict[str | None, list[str]] = defaultdict(list)  # ftype (or None) -> [opt_name]
    outfile_opts: list[str] = []
    novalue_opts: list[str] = []   # int/real/string/int_list/real_list

    for o in s.options:
        a = o.arg
        if a.type == "flag":
            continue
        if a.type == "enum":
            enum_groups[a.values].append(o.name)
        elif a.type == "file":
            file_groups[a.ftype].append(o.name)
        elif a.type == "outfile":
            outfile_opts.append(o.name)
        else:  # int, real, string, int_list, real_list
            novalue_opts.append(o.name)

    out = [HEADER, "_vsearch() {", "    local cur prev cmd opts i",
           "    COMPREPLY=()",
           '    cur="${COMP_WORDS[COMP_CWORD]}"',
           '    prev="${COMP_WORDS[COMP_CWORD-1]}"',
           "",
           "    # ---- Detect which command (if any) is on the line ----",
           '    cmd=""',
           "    for ((i=1; i<COMP_CWORD; i++)); do",
           '        case "${COMP_WORDS[i]}" in']

    out.append(f"            {fmt_pipe(s.commands)})")
    out.append('                cmd="${COMP_WORDS[i]}"')
    out.append("                break ;;")
    out.append("        esac")
    out.append("    done")
    out.append("")

    # ---- Argument-value completion based on $prev ----
    out.append("    # ---- Complete the value for a value-taking option ----")
    out.append('    case "$prev" in')

    # Enums
    for enum_key, names in sorted(enum_groups.items()):
        values = " ".join(s.enum_values(enum_key))
        out.append(f"        {fmt_pipe(sorted(names))})")
        out.append(f'            COMPREPLY=( $(compgen -W "{values}" -- "$cur") )')
        out.append("            return ;;")

    # Files grouped by ftype
    for ftype, names in sorted(file_groups.items(), key=lambda kv: (kv[0] or "~")):
        if ftype is None:
            ext_pat = ""
            comment = "any file"
        else:
            ext_pat = "|".join(s.extensions(ftype))
            comment = f"{ftype} input"
        out.append(f"        {fmt_pipe(sorted(names))})  # {comment}")
        out.append(f'            _vsearch_filedir "{ext_pat}"')
        out.append("            return ;;")

    # Outfiles
    if outfile_opts:
        out.append(f"        {fmt_pipe(sorted(outfile_opts))})  # output file")
        out.append('            _vsearch_filedir ""')
        out.append("            return ;;")

    # Numeric/string/list — no useful completion, but consume the next slot
    if novalue_opts:
        out.append(f"        {fmt_pipe(sorted(novalue_opts))})")
        out.append("            return ;;")

    out.append("    esac")
    out.append("")

    # ---- Option name suggestions ----
    # Reached when $prev didn't trigger value completion. We suggest options
    # regardless of whether $cur is empty or already starts with "-", so that
    # bare <Tab> shows the available choices (matching zsh / fish behaviour).
    out.append("    # ---- Suggest option names ----")
    out.append('    case "$cmd" in')

    for command in s.commands:
        visible = [o.name for o in s.options_for_command(command)]
        # Include the command itself in the list (so re-typing the command flag still completes).
        visible_set = sorted(set(visible))
        out.append(f"        {command})")
        out.append(f'            opts="{" ".join(visible_set)}";;')

    # No command yet: show all action commands (including --help / --version,
    # which are commands too). Short forms (-h, -v) are pulled in from any
    # command option's `short:` attribute.
    pre_opts: list[str] = list(s.commands)
    for cmd_name in s.commands:
        opt = next(o for o in s.options if o.name == cmd_name)
        if opt.short:
            pre_opts.append(opt.short)
    pre_opts = sorted(set(pre_opts))
    out.append("        *)")
    out.append(f'            opts="{" ".join(pre_opts)}";;')
    out.append("    esac")
    out.append('    COMPREPLY=( $(compgen -W "$opts" -- "$cur") )')
    out.append("}")
    out.append(FOOTER)

    sys.stdout.write("\n".join(out))


if __name__ == "__main__":
    main()
