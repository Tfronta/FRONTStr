"""BED serialization of a :class:`Panel`.

Five columns::

    chrom  start  end  motif[,motif...]  name

Coordinates are 1-based inclusive, written verbatim from the panel YAML (which
already uses that convention, as does HipSTR/LongTR's BED dialect). Markers
without a motif are skipped.

Lives in ``panel`` rather than ``caller``: this is how a panel serializes, not
something an external caller needs. The report embeds it so a reader can see
the exact windows a run used — and paste them straight into samtools or IGV —
rather than taking the marker names on trust.
"""

from __future__ import annotations

from collections import defaultdict
from collections.abc import Iterable
from pathlib import Path

from frontstr.errors import PanelError
from frontstr.panel.models import Panel, System


def write_panel_bed(panel: Panel, out_path: Path) -> Path:
    """Emit the full panel as a single BED file.

    Args:
        panel: Forensic panel to export.
        out_path: Destination ``.bed`` file. Overwritten if it exists.

    Returns:
        ``out_path`` after writing.

    Raises:
        PanelError: If no usable markers were found.
    """
    rows = _format_rows(panel.systems)
    if not rows:
        raise PanelError(f"Panel {panel.name!r} has no markers with motifs; cannot emit BED")
    out_path.write_text("\n".join(rows) + "\n")
    return out_path


def split_panel_by_chromosome(panel: Panel, out_dir: Path) -> dict[str, Path]:
    """Emit one BED per chromosome present in ``panel``.

    Useful when running LongTR in parallel via its ``--chrom`` flag.

    Returns:
        Mapping ``{chromosome: path}``.
    """
    out_dir.mkdir(parents=True, exist_ok=True)
    by_chrom: dict[str, list[str]] = defaultdict(list)
    for s in panel.systems:
        row = _format_row(s)
        if row is not None:
            by_chrom[s.chromosome].append(row)
    paths: dict[str, Path] = {}
    for chrom, rows in by_chrom.items():
        p = out_dir / f"panel.{chrom}.bed"
        p.write_text("\n".join(rows) + "\n")
        paths[chrom] = p
    return paths


def panel_bed_lines(panel: Panel) -> list[str]:
    """The panel's windows as BED lines, for embedding rather than writing.

    Same content :func:`write_panel_bed` puts on disk. The report shows it
    verbatim so a reader can see the exact intervals a run used, and paste them
    into samtools or IGV, instead of taking the marker names on trust.
    """
    return _format_rows(panel.systems)


def _format_rows(systems: Iterable[System]) -> list[str]:
    rows: list[str] = []
    for s in sorted(systems, key=lambda x: (x.chromosome, x.ref_start)):
        line = _format_row(s)
        if line is not None:
            rows.append(line)
    return rows


def _format_row(s: System) -> str | None:
    if not s.motif:
        return None
    return "\t".join((s.chromosome, str(s.ref_start), str(s.ref_end), s.motif, s.name))
