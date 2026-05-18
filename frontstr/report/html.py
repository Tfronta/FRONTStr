"""Render the FRONTStr HTML report.

The report is a single self-contained ``.html`` file:

- CSS and JS bundles are read from :mod:`frontstr.report.static` and inlined.
- SVG charts are pre-rendered server-side (no JS dependency for data).
- The full payload is embedded as ``<script type="application/json" id="run-data">``
  for programmatic extraction.
- A SHA-256 of the rendered HTML is computed *after* rendering and substituted
  into a placeholder, so the report contains its own integrity hash.

Public entry point: :func:`build_report`.
"""

from __future__ import annotations

import hashlib
import json
from importlib import resources
from pathlib import Path
from typing import Any

from jinja2 import Environment, FileSystemLoader, select_autoescape

from frontstr.errors import FrontstrError
from frontstr.interp.models import MarkerResult
from frontstr.report.payload import RunContext, serialize_run
from frontstr.report.svg_charts import (
    allele_coverage_svg,
    coverage_bar_svg,
    haplotype_stack_svg,
)

REPORT_HASH_PLACEHOLDER = "@@REPORT_SHA256@@"
_TEMPLATE_NAME = "run_report.html.j2"


def build_report(
    results: list[MarkerResult],
    context: RunContext,
    out_path: Path,
) -> Path:
    """Write a self-contained HTML report to ``out_path`` and return its path.

    Args:
        results: One :class:`MarkerResult` per marker, in panel order.
        context: Run-level metadata.
        out_path: Destination ``.html`` file.

    Raises:
        FrontstrError: If the template cannot be located or rendered.
    """
    payload = serialize_run(results, context)
    rendered = _render(payload)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(rendered, encoding="utf-8")
    _stamp_self_hash(out_path)
    return out_path


def _render(payload: dict[str, Any]) -> str:
    """Render the Jinja template against ``payload``."""
    template_dir = _locate_template_dir()
    env = Environment(
        loader=FileSystemLoader(str(template_dir)),
        autoescape=select_autoescape(("html",)),
        keep_trailing_newline=True,
    )
    try:
        template = env.get_template(_TEMPLATE_NAME)
    except Exception as exc:
        raise FrontstrError(f"Failed to load report template: {exc}") from exc

    css_bundle = _load_static("styles.css")
    app_js_bundle = _load_static("app.js")

    # Attach pre-rendered SVGs to each result so the template can ``| safe`` them.
    enriched_results: list[dict[str, Any]] = []
    for r in payload.get("results", []):
        enriched = dict(r)
        n_groups = len(enriched.get("ngs_panel", {}).get("chart_groups") or [])
        chart_height = min(180, max(120, 88 + 18 * max(n_groups, 1)))
        enriched["svg_allele_coverage"] = allele_coverage_svg(r, height=chart_height)
        enriched["svg_haplotype"] = haplotype_stack_svg(r)
        enriched_results.append(enriched)
    payload = {**payload, "results": enriched_results}

    qc_coverage_svg = coverage_bar_svg(
        payload["qc"]["coverage_table"],
        floor=context_dropout(payload),
    )

    run_data_json = json.dumps(payload, separators=(",", ":"), default=str)
    run_data_json_pretty = json.dumps(payload, indent=2, default=str)

    return template.render(
        payload=payload,
        css_bundle=css_bundle,
        app_js_bundle=app_js_bundle,
        qc_coverage_svg=qc_coverage_svg,
        run_data_json=run_data_json,
        run_data_json_pretty=run_data_json_pretty,
    )


def context_dropout(payload: dict[str, Any]) -> int:
    """Read the dropout floor from the serialized payload (defaults to 30)."""
    return int(payload.get("meta", {}).get("dropout_floor", 30) or 30)


def _stamp_self_hash(out_path: Path) -> None:
    """Replace ``@@REPORT_SHA256@@`` with the SHA-256 of the rendered HTML.

    The hash is computed over the file content with placeholders unchanged, so
    the resulting file's SHA differs from this stamped value — by design (the
    stamp says "this was the hash before stamping", which is what you want to
    pin once the file is final).
    """
    contents = out_path.read_text(encoding="utf-8")
    pre_stamp = contents.encode("utf-8")
    digest = hashlib.sha256(pre_stamp).hexdigest()
    stamped = contents.replace(REPORT_HASH_PLACEHOLDER, digest)
    out_path.write_text(stamped, encoding="utf-8")


def _load_static(name: str) -> str:
    """Read a vendored asset from :mod:`frontstr.report.static`."""
    try:
        return (resources.files("frontstr.report") / "static" / name).read_text(
            encoding="utf-8"
        )
    except (FileNotFoundError, ModuleNotFoundError, AttributeError) as exc:
        # Editable installs sometimes can't resolve via importlib.resources
        # → fall back to the filesystem location of this module.
        fallback = Path(__file__).parent / "static" / name
        if fallback.exists():
            return fallback.read_text(encoding="utf-8")
        raise FrontstrError(f"Static asset {name!r} not found: {exc}") from exc


def _locate_template_dir() -> Path:
    """Return the directory where Jinja templates live."""
    here = Path(__file__).parent / "templates"
    if here.is_dir():
        return here
    try:
        path = resources.files("frontstr.report") / "templates"
        if path.is_dir():
            return Path(str(path))
    except (FileNotFoundError, ModuleNotFoundError):
        pass
    raise FrontstrError("Could not locate FRONTStr report templates directory")
