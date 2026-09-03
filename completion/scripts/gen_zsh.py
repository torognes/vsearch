"""Generate the zsh completion script for vsearch.

Writes to stdout. Invoked from completion/Makefile.am:
    cd scripts && python3 gen_zsh.py > ../completions/_vsearch
"""

from __future__ import annotations

import sys

import spec as spec_mod


HEADER = """\
#compdef vsearch
# vsearch completion for zsh
# Auto-generated from spec/vsearch.yaml — do not edit by hand.
#
# Install: place this file (named "_vsearch") somewhere in your $fpath,
# e.g. /usr/local/share/zsh/site-functions/_vsearch or ~/.zsh/completions/_vsearch
# (and add that directory to $fpath before `compinit` in ~/.zshrc).
"""


def arg_desc(s: str) -> str:
    """Escape characters that have meaning inside an `_arguments` `[...]`
    description bracket. Shell quoting is handled separately by sh_squote."""
    return (
        s.replace("\\", "\\\\")
         .replace("[", "\\[")
         .replace("]", "\\]")
         .replace(":", "\\:")
    )


def sh_squote(s: str) -> str:
    """Wrap a string in single quotes, escaping embedded single quotes."""
    return "'" + s.replace("'", "'\\''") + "'"


def file_glob(extensions: list[str]) -> str:
    return "*.(" + "|".join(extensions) + ")"


def arg_action(s: spec_mod.Spec, opt: spec_mod.Option) -> str:
    """Return the `:label:action` suffix for a value-taking option, or "" for flags."""
    a = opt.arg
    if a.type == "flag":
        return ""
    if a.type == "enum":
        values = " ".join(s.enum_values(a.values))
        return f":value:({values})"
    if a.type == "file":
        if a.ftype is None:
            return ":file:_files"
        glob = file_glob(s.extensions(a.ftype))
        return f":file:_files -g '{glob}'"
    if a.type == "outfile":
        return ":file:_files"
    labels = {
        "int": "n",
        "real": "x",
        "string": "string",
        "int_list": "n,n,...",
        "real_list": "x,x,...",
    }
    return f":{labels.get(a.type, 'value')}:"


def format_option(s: spec_mod.Spec, opt: spec_mod.Option) -> list[str]:
    """Return one or two `_arguments` specs (one for --long, one for -short)."""
    action = arg_action(s, opt)
    d = arg_desc(opt.desc)
    out = [sh_squote(f"{opt.name}[{d}]{action}")]
    if opt.short:
        out.append(sh_squote(f"{opt.short}[{d}]{action}"))
    return out


def emit_block(label: str, specs: list[str], indent: str = "        ") -> list[str]:
    """Emit a single `_arguments -s ...` call indented for a case branch."""
    if not specs:
        return [f"{indent}_arguments -s"]
    lines = [f"{indent}_arguments -s \\"]
    for i, sp in enumerate(specs):
        suffix = " \\" if i < len(specs) - 1 else ""
        lines.append(f"{indent}    {sp}{suffix}")
    return lines


def main() -> None:
    s = spec_mod.load()
    by_name = {o.name: o for o in s.options}

    out: list[str] = [HEADER]
    out.append("_vsearch() {")
    out.append("    local cmd w")
    out.append('    cmd=""')
    out.append("    # Walk the command line for an action flag")
    out.append('    for w in "${words[@]:1:$CURRENT-1}"; do')
    out.append("        case \"$w\" in")
    out.append("            " + "|".join(s.commands) + ")")
    out.append('                cmd="$w"; break ;;')
    out.append("        esac")
    out.append("    done")
    out.append("")
    out.append('    case "$cmd" in')

    # One branch per command
    for command in s.commands:
        opts = s.options_for_command(command)
        # Also include the command itself first
        cmd_opt = by_name[command]
        specs: list[str] = []
        specs.extend(format_option(s, cmd_opt))
        for o in opts:
            if o.name == command:
                continue
            specs.extend(format_option(s, o))
        out.append(f"        {command})")
        out.extend(emit_block(command, specs))
        out.append("            ;;")

    # No command yet: all action commands (--help / --version included
    # as commands, since vsearch treats them as such). Cross-cutting
    # options like --log / --quiet are surfaced once an action flag is on
    # the line (via command_options).
    pre: list[spec_mod.Option] = [by_name[c] for c in s.commands]
    pre_specs: list[str] = []
    for o in pre:
        pre_specs.extend(format_option(s, o))
    out.append("        *)")
    out.extend(emit_block("(no command)", pre_specs))
    out.append("            ;;")

    out.append("    esac")
    out.append("}")
    out.append("")
    out.append('_vsearch "$@"')
    out.append("")

    sys.stdout.write("\n".join(out))


if __name__ == "__main__":
    main()
