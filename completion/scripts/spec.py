"""Load and validate the vsearch completion spec.

Provides a small typed model on top of spec/vsearch.yaml that the
per-shell generators consume.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

import yaml


SPEC_PATH = Path(__file__).resolve().parent.parent / "spec" / "vsearch.yaml"


@dataclass(frozen=True)
class Arg:
    type: str
    ftype: Optional[str] = None
    values: Optional[str] = None


@dataclass(frozen=True)
class Option:
    name: str
    arg: Arg
    short: Optional[str] = None
    command: bool = False
    hidden: bool = False
    desc: str = ""


@dataclass
class Spec:
    enums: dict[str, list[str]]
    file_types: dict[str, list[str]]
    commands: list[str]
    options: list[Option]
    command_options: dict[str, list[str]]

    def __post_init__(self) -> None:
        self._by_name = {o.name: o for o in self.options}

    def get(self, name: str) -> Option:
        return self._by_name[name]

    def visible_options(self) -> list[Option]:
        return [o for o in self.options if not o.hidden]

    def general_options(self) -> list[Option]:
        """Options that appear in *every* command's allowed list.

        These are the options safe to suggest before any command has been
        chosen, together with --help / --version which are not in the
        command_options mapping at all.
        """
        if not self.command_options:
            return []
        common = set.intersection(*(set(v) for v in self.command_options.values()))
        return [o for o in self.visible_options() if o.name in common]

    def options_for_command(self, command: str) -> list[Option]:
        names = self.command_options.get(command, [])
        return [self.get(n) for n in names if not self.get(n).hidden]

    def enum_values(self, key: str) -> list[str]:
        return self.enums[key]

    def extensions(self, ftype: str) -> list[str]:
        return self.file_types[ftype]


def load(path: Path | str | None = None) -> Spec:
    path = Path(path) if path else SPEC_PATH
    raw = yaml.safe_load(path.read_text())

    file_types = {k: list(v["extensions"]) for k, v in raw["file_types"].items()}

    options: list[Option] = []
    for o in raw["options"]:
        arg = Arg(
            type=o["arg"]["type"],
            ftype=o["arg"].get("ftype"),
            values=o["arg"].get("values"),
        )
        options.append(
            Option(
                name=o["name"],
                arg=arg,
                short=o.get("short"),
                command=o.get("command", False),
                hidden=o.get("hidden", False),
                desc=o.get("desc", ""),
            )
        )

    return Spec(
        enums={k: list(v) for k, v in raw["enums"].items()},
        file_types=file_types,
        commands=list(raw["commands"]),
        options=options,
        command_options={k: list(v) for k, v in raw["command_options"].items()},
    )
