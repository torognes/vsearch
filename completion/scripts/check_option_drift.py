"""Check the completion spec against the option tables in src/cli.cc.

The spec is hand-maintained and cli.cc is the authority, so the two drift
apart whenever an option is added, removed or reassigned. A stale spec
does not merely under-suggest: it offers options the binary rejects.

Checked here:

  1. every option in option_specs is in the spec, and vice versa
  2. needs_arg agrees with the spec's argument type
  3. every command row in valid_options is a spec command
  4. the per-command option lists match exactly, both directions
  5. every option the spec hides still exists in cli.cc

Deliberately NOT checked, because cli.cc does not carry it -- these stay
human-maintained, and the manual pages are their reference:

  - the description text
  - the argument type beyond value-or-flag (int, float, string, file,
    outfile, enum) and the file-extension filter
  - the enum value lists
  - which options are hidden (that is an editorial choice)

Usage: python3 check_option_drift.py <cli.json> <spec.yaml>
where cli.json is the output of extract_cli_options.py.
"""

from __future__ import annotations

import json
import sys

import yaml

# vsearch has no short options: getopt_long_only is given an empty short
# string, so --h and --v are one-letter LONG options (which is also why
# -h and -v happen to work). cli.cc therefore lists them in option_specs
# and gives each its own valid_options row, while the spec models them as
# the `short:` field of --help and --version. Map them rather than skip
# them, so a genuine disappearance still shows up.
ALIASES = {'--h': '--help', '--v': '--version'}


def load(cli_path: str, spec_path: str) -> tuple:
    with open(cli_path, encoding='utf-8') as handle:
        cli = json.load(handle)
    with open(spec_path, encoding='utf-8') as handle:
        spec = yaml.safe_load(handle)

    options = {ALIASES.get(name, name): needs_arg
               for name, needs_arg in cli['options'].items()}
    commands = {}
    for command, accepted in cli['commands'].items():
        # An alias row accepts the same options as the row it aliases;
        # merge them so --help is compared once.
        name = ALIASES.get(command, command)
        merged = set(commands.get(name, ())) | {ALIASES.get(o, o) for o in accepted}
        commands[name] = sorted(merged)
    return options, commands, spec


def main() -> int:
    if len(sys.argv) != 3:
        sys.exit(f'usage: {sys.argv[0]} <cli.json> <spec.yaml>')

    cli_options, cli_commands, spec = load(sys.argv[1], sys.argv[2])
    spec_options = {option['name']: option for option in spec['options']}
    spec_commands = set(spec['commands'])
    spec_accepts = spec['command_options']

    problems = []

    def report(label: str, items) -> None:
        items = sorted(items)
        print(f'{label}: {len(items)}')
        for item in items:
            print(f'    {item}')
        if items:
            problems.append(label)

    print('counts')
    print(f'  cli.cc: {len(cli_options)} options, {len(cli_commands)} commands, '
          f'{sum(len(v) for v in cli_commands.values())} command/option pairs')
    print(f'  spec:   {len(spec_options)} options, {len(spec_commands)} commands, '
          f'{sum(len(v) for v in spec_accepts.values())} command/option pairs')
    print()

    print('options')
    report('  in cli.cc but missing from the spec',
           set(cli_options) - set(spec_options))
    report('  in the spec but unknown to cli.cc',
           set(spec_options) - set(cli_options))

    mismatched = []
    for name, needs_arg in cli_options.items():
        option = spec_options.get(name)
        if option is None:
            continue
        spec_type = option.get('arg', {}).get('type')
        if (spec_type == 'flag') == needs_arg:
            mismatched.append(f'{name}: cli.cc needs_arg={needs_arg}, '
                              f'spec type={spec_type}')
    report('  argument requirement disagrees', mismatched)

    report('  hidden by the spec but gone from cli.cc',
           {option['name'] for option in spec['options'] if option.get('hidden')}
           - set(cli_options))
    print()

    print('commands')
    report('  in cli.cc but missing from the spec', set(cli_commands) - spec_commands)
    report('  in the spec but unknown to cli.cc', spec_commands - set(cli_commands))
    print()

    print('per-command option lists')
    for command in sorted(set(cli_commands) & spec_commands):
        accepted = set(cli_commands[command])
        listed = set(spec_accepts.get(command, ()))
        missing, extra = accepted - listed, listed - accepted
        if missing or extra:
            problems.append(command)
            print(f'  {command}')
            for option in sorted(missing):
                print(f'    missing: {option}')
            for option in sorted(extra):
                print(f'    extra:   {option}')
    if not problems:
        print('  all match')

    print()
    if problems:
        print(f'FAIL: {len(problems)} kind(s) of drift. Edit '
              'completion/spec/vsearch.yaml, never completions/, then '
              'regenerate with "make -C completion regenerate-completions".')
        return 1
    print('OK: the spec matches cli.cc')
    return 0


if __name__ == '__main__':
    sys.exit(main())
