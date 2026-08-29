"""Turn a params model into something a template can render.

The Pydantic model stays the single definition of what a parameter is: its type, bounds,
default and help text all come from there, so the form and the validation cannot disagree.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

from pydantic import BaseModel

from .params import TELOTAGS


@dataclass
class Field_:
    name: str
    label: str
    kind: str                      # text | number | select | checkbox
    default: Any = None
    help: str = ""
    required: bool = False
    choices: list[str] = field(default_factory=list)
    minimum: float | None = None
    maximum: float | None = None
    #: Rendered read-only because something else determines it -- currently only chr_arm_ln,
    #: which is read from the registered reference. Letting anyone type a value that disagrees
    #: with the reference silently corrupts q-arm telomere lengths.
    derived_from: str = ""
    advanced: bool = False


#: Everything else is either obvious or safe to leave alone; the paper's thresholds in
#: particular should not be casually retuned, so they go behind a disclosure.
PRIMARY = {"sample_label", "barcode_name", "workers", "threads"}


def describe(model: type[BaseModel]) -> list[Field_]:
    fields: list[Field_] = []
    for name, info in model.model_fields.items():
        annotation = info.annotation
        kind = "number" if annotation in (int, float) else "text"
        choices: list[str] = []

        if name == "barcode_name":
            kind, choices = "select", sorted(TELOTAGS)
        elif annotation is bool:
            kind = "checkbox"

        metadata = {type(m).__name__: m for m in info.metadata}
        fields.append(Field_(
            name=name,
            label=name.replace("_", " ").capitalize(),
            kind=kind,
            default=None if info.is_required() else info.default,
            help=info.description or "",
            required=info.is_required(),
            choices=choices,
            minimum=getattr(metadata.get("Ge"), "ge", None),
            maximum=getattr(metadata.get("Le"), "le", None),
            derived_from="reference" if name == "chr_arm_ln" else "",
            advanced=name not in PRIMARY,
        ))
    return fields


def coerce(model: type[BaseModel], raw: dict[str, str]) -> dict:
    """Pull the model's fields out of form data, dropping blanks so defaults apply."""
    values: dict[str, Any] = {}
    for name, info in model.model_fields.items():
        if name not in raw:
            continue
        value = raw[name]
        if isinstance(value, str):
            value = value.strip()
        if value in ("", None):
            continue
        if info.annotation is bool:
            value = value not in ("0", "false", "off")
        values[name] = value
    return values
