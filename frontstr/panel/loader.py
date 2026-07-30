"""Load panels from YAML files."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import yaml
from pydantic import ValidationError

from frontstr.errors import PanelError
from frontstr.panel.models import Panel


def load_panel(path: Path) -> Panel:
    """Load a panel YAML file into a :class:`Panel` instance.

    Raises:
        PanelError: if the file is missing, malformed, or fails validation.
    """
    if not path.exists():
        raise PanelError(f"Panel YAML not found: {path}")
    try:
        raw: Any = yaml.safe_load(path.read_text())
    except yaml.YAMLError as exc:
        raise PanelError(f"Cannot parse panel YAML {path}: {exc}") from exc
    if not isinstance(raw, dict):
        raise PanelError(f"Panel YAML root must be a mapping; got {type(raw).__name__}")
    try:
        return Panel.model_validate(raw)
    except ValidationError as exc:
        raise PanelError(f"Panel {path} failed validation:\n{exc}") from exc
