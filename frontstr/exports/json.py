"""JSON export of a FRONTStr run.

We re-use :func:`frontstr.report.payload.serialize_run` so the JSON file
contains byte-identical content with the JSON blob embedded in the HTML
report. This is the canonical "machine-readable run record".

Two output modes:

- ``"pretty"`` (default): 2-space indent, sorted keys-by-section preserved
  as written by the serializer. Stable for diffs.
- ``"compact"``: separators=``(",", ":")`` for minimal size (matches the
  HTML embedded payload).
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Literal

from frontstr.errors import FrontstrError

JsonMode = Literal["pretty", "compact"]


def write_run_json(
    payload: dict[str, Any],
    out_path: Path,
    *,
    mode: JsonMode = "pretty",
) -> Path:
    """Write the serialized run payload to ``out_path``.

    Args:
        payload: Output of :func:`frontstr.report.payload.serialize_run`.
        out_path: Destination ``.json`` file.
        mode: ``"pretty"`` (indent=2) or ``"compact"``.

    Returns:
        ``out_path`` after writing.

    Raises:
        FrontstrError: If ``mode`` is unknown or writing fails.
    """
    if mode == "pretty":
        text = json.dumps(payload, indent=2, default=str, ensure_ascii=False)
    elif mode == "compact":
        text = json.dumps(payload, separators=(",", ":"), default=str, ensure_ascii=False)
    else:
        raise FrontstrError(f"Unknown JSON mode: {mode!r}")
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(text + ("\n" if mode == "pretty" else ""), encoding="utf-8")
    return out_path
