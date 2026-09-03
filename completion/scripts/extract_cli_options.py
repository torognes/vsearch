"""Extract vsearch's command-line option tables from src/cli.cc.

The two tables in cli.cc are the single source of truth for the option
metadata the completion spec has to agree with:

  option_specs   one entry per option: its name, and whether it takes a
                 value (needs_arg)
  valid_options  one row per command, the first entry naming the command
                 and the rest listing the options it accepts, terminated
                 by -1

Writes JSON to stdout:

  {"options": {"--name": needs_arg, ...},
   "commands": {"--command": ["--option", ...], ...}}

Usage: python3 extract_cli_options.py <path to cli.cc>
"""

from __future__ import annotations

import json
import re
import sys


def strip_comments(text: str) -> str:
    """Remove C and C++ comments, so commented-out entries are not read."""
    text = re.sub(r'/\*.*?\*/', '', text, flags=re.S)
    return re.sub(r'//[^\n]*', '', text)


def table_body(text: str, anchor: str, path: str) -> str:
    """The text between anchor's opening '{{' and its closing '},};'."""
    start = text.find(anchor)
    if start < 0:
        sys.exit(f'{path}: no {anchor.rstrip(" =")} table found. The table was '
                 'renamed or removed; this script needs updating.')
    opening = text.find('{{', start)
    closing = text.find('},};', opening)
    if opening < 0 or closing < 0:
        sys.exit(f'{path}: could not delimit the {anchor.rstrip(" =")} table. '
                 'It no longer opens with "{{" and closes with "},};"; this '
                 'script needs updating.')
    return text[opening + 2:closing]


def extract(path: str) -> dict:
    with open(path, encoding='utf-8') as handle:
        text = strip_comments(handle.read())

    options = {
        '--' + name: needs_arg == 'true'
        for name, needs_arg in re.findall(
            r'\{\s*"([^"]+)"\s*,\s*(true|false)\s*\}',
            table_body(text, 'option_specs =', path))
    }
    if not options:
        sys.exit(f'{path}: option_specs parsed as empty.')

    # Rows are terminated by -1 and no option index is negative, so
    # splitting on the literal is safe.
    commands = {}
    for row in table_body(text, 'valid_options =', path).split('-1'):
        names = re.findall(r'option_([A-Za-z0-9_]+)', row)
        if names:
            commands['--' + names[0]] = sorted('--' + name for name in names[1:])
    if not commands:
        sys.exit(f'{path}: valid_options parsed as empty.')

    return {'options': options, 'commands': commands}


def main() -> None:
    if len(sys.argv) != 2:
        sys.exit(f'usage: {sys.argv[0]} <path to cli.cc>')
    json.dump(extract(sys.argv[1]), sys.stdout, indent=1, sort_keys=True)
    sys.stdout.write('\n')


if __name__ == '__main__':
    main()
