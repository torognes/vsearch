"""Generate the fish completion script for vsearch.

Writes to stdout. Invoked from completion/Makefile.am:
    cd scripts && python3 gen_fish.py > ../completions/vsearch.fish
"""

from __future__ import annotations

import sys

import spec as spec_mod


HEADER = """\
# vsearch completion for fish
# Auto-generated from spec/vsearch.yaml — do not edit by hand.
#
# Install: copy or symlink this file as one of
#   ~/.config/fish/completions/vsearch.fish     (per-user)
#   /usr/share/fish/vendor_completions.d/vsearch.fish  (system-wide)

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

function __vsearch_current_command --description 'Echo the action flag on the line, if any'
    set -l cmds {COMMANDS_SPACE_SEPARATED}
    for word in (commandline -opc)
        if contains -- $word $cmds
            echo $word
            return 0
        end
    end
    return 1
end

function __vsearch_no_command
    not __vsearch_current_command >/dev/null 2>&1
end

function __vsearch_command_is --argument-names cmd
    set -l current (__vsearch_current_command 2>/dev/null)
    or return 1
    test "$current" = "$cmd"
end

"""


def fish_squote(s: str) -> str:
    """Wrap in single quotes, escaping inner single quotes."""
    return "'" + s.replace("\\", "\\\\").replace("'", "\\'") + "'"


def file_arg_flags(s: spec_mod.Spec, ftype: str | None) -> str:
    """Return the trailing flags for a file-taking option.

    Without ftype: -r (require arg, default file completion).
    With ftype:    -rfa '(__fish_complete_suffix .ext1 .ext2 ...)'
                   -r requires arg, -f suppresses default files, -a adds suffix matches.
    """
    if ftype is None:
        return "-r"
    exts = " ".join("." + e for e in s.extensions(ftype))
    return f"-rfa '(__fish_complete_suffix {exts})'"


def option_flags(s: spec_mod.Spec, opt: spec_mod.Option) -> str:
    """Return the per-option fish complete flags excluding -c/-n/-l/-s/-d."""
    a = opt.arg
    if a.type == "flag":
        return "-f"  # no arg, no file completion
    if a.type == "enum":
        vals = " ".join(s.enum_values(a.values))
        return f"-x -a {fish_squote(vals)}"
    if a.type == "file":
        return file_arg_flags(s, a.ftype)
    if a.type == "outfile":
        return "-rF"  # require arg, force file completion (no extension filter)
    # int / real / string / int_list / real_list — require arg, no suggestions
    return "-x"


def complete_line(s: spec_mod.Spec, opt: spec_mod.Option, condition: str | None) -> str:
    parts = ["complete -c vsearch"]
    if condition:
        parts.append(f"-n {fish_squote(condition)}")
    parts.append(f"-l {opt.name.lstrip('-')}")
    if opt.short:
        parts.append(f"-s {opt.short.lstrip('-')}")
    parts.append(option_flags(s, opt))
    if opt.desc:
        parts.append(f"-d {fish_squote(opt.desc)}")
    return " ".join(parts)


def main() -> None:
    s = spec_mod.load()
    by_name = {o.name: o for o in s.options}

    header = HEADER.replace("{COMMANDS_SPACE_SEPARATED}", " ".join(s.commands))
    out: list[str] = [header]

    # ---- Pre-command suggestions (no command on line) ----
    out.append("# ---------------------------------------------------------------------------")
    out.append("# Pre-command suggestions (no command flag on the line yet)")
    out.append("# ---------------------------------------------------------------------------")

    # Each command flag (including --help / --version, which are commands
    # in their own right). Cross-cutting options like --log / --quiet are
    # not commands and are surfaced via the per-command sections below
    # once an action flag is on the line.
    out.append("# Command (action) flags")
    for name in s.commands:
        out.append(complete_line(s, by_name[name], "__vsearch_no_command"))

    # ---- Per-command option suggestions ----
    out.append("")
    out.append("# ---------------------------------------------------------------------------")
    out.append("# Per-command options")
    out.append("# ---------------------------------------------------------------------------")

    for command in s.commands:
        condition = f"__vsearch_command_is {command}"
        out.append("")
        out.append(f"# --- {command} ---")
        # The command flag itself (so fish completes its file arg even after it's typed)
        out.append(complete_line(s, by_name[command], condition))
        for o in s.options_for_command(command):
            if o.name == command:
                continue
            out.append(complete_line(s, o, condition))

    sys.stdout.write("\n".join(out) + "\n")


if __name__ == "__main__":
    main()
