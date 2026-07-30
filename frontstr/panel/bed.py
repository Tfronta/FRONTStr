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

import re
from collections import defaultdict
from collections.abc import Iterable
from pathlib import Path
from typing import Literal

from frontstr.errors import PanelError
from frontstr.panel.models import Panel, System

#: How a BED file's coordinates are meant.
#:
#: ``bed0`` is what BED means everywhere else: 0-based, half-open, so
#: ``chr1 100 110`` is the eleventh through twentieth base. This is what
#: ``samtools view -L`` and UCSC expect.
#:
#: ``panel1`` is 1-based inclusive — the convention the panel YAML uses and
#: what HipSTR/LongTR's region files use. Writing it verbatim is why
#: :func:`write_panel_bed` exists at all.
#:
#: There is no default anywhere in this module. The two differ by one base at
#: the start, which at a window edge is exactly the kind of error that produces
#: a plausible wrong answer instead of a crash.
BedCoords = Literal["bed0", "panel1"]


def write_panel_bed(panel: Panel, out_path: Path, *, coords: BedCoords = "panel1") -> Path:
    """Emit the full panel as a single BED file.

    Args:
        panel: Forensic panel to export.
        out_path: Destination ``.bed`` file. Overwritten if it exists.
        coords: Coordinate convention to write. Defaults to ``panel1`` for the
            LongTR region-file path this was written for; pass ``bed0`` for a
            file any other tool will read correctly.

    Returns:
        ``out_path`` after writing.

    Raises:
        PanelError: If no usable markers were found.
    """
    rows = _format_rows(panel.systems, coords)
    if not rows:
        raise PanelError(f"Panel {panel.name!r} has no markers with motifs; cannot emit BED")
    out_path.write_text("\n".join(rows) + "\n")
    return out_path


def split_panel_by_chromosome(
    panel: Panel, out_dir: Path, *, coords: BedCoords = "panel1"
) -> dict[str, Path]:
    """Emit one BED per chromosome present in ``panel``.

    Useful when running LongTR in parallel via its ``--chrom`` flag, which is
    why this too defaults to ``panel1``.

    Returns:
        Mapping ``{chromosome: path}``.
    """
    out_dir.mkdir(parents=True, exist_ok=True)
    by_chrom: dict[str, list[str]] = defaultdict(list)
    for s in panel.systems:
        row = _format_row(s, coords)
        if row is not None:
            by_chrom[s.chromosome].append(row)
    paths: dict[str, Path] = {}
    for chrom, rows in by_chrom.items():
        p = out_dir / f"panel.{chrom}.bed"
        p.write_text("\n".join(rows) + "\n")
        paths[chrom] = p
    return paths


def panel_bed_lines(panel: Panel, *, coords: BedCoords = "bed0") -> list[str]:
    """The panel's windows as BED lines, for embedding rather than writing.

    Defaults to ``bed0`` — standard BED — because the audience is a person
    pasting the block into samtools or IGV, and every tool but HipSTR reads BED
    that way. :func:`load_panel_from_bed` reads the same convention, so the
    round trip is exact.
    """
    return _format_rows(panel.systems, coords)


#: A 4th column that is nothing but ACGT and commas is a motif, not a name.
_MOTIF_RE = re.compile(r"^[ACGTacgt]+(,[ACGTacgt]+)*$")


def load_panel_from_bed(
    path: Path,
    *,
    coords: BedCoords = "bed0",
    name: str = "bed",
) -> tuple[Panel, list[str]]:
    """Build a :class:`Panel` from a BED file of regions.

    The point is the one HipSTR and LongTR get right and a YAML-only caller does
    not: a user should be able to point the tool at their own regions without
    curating a panel first.

    Accepted columns::

        chrom  start  end  [motif]  [name]

    A 4th column of nothing but ACGT and commas is read as the motif; anything
    else is read as the name. **The motif is required**, and that is not a
    parsing convenience: without it reads cannot be binned by repeat-core
    length, which is the mechanism that stops ONT flank indel errors splitting
    one allele into a dozen clusters. Binning on raw window length instead took
    TH01 from 2 bins to 12. Better to refuse the file than to call from it
    badly.

    Args:
        path: BED file.
        coords: How to read the coordinates — see :data:`BedCoords`. Defaults to
            standard ``bed0``, which is what any BED not written by HipSTR will
            be.
        name: Panel name to record.

    Returns:
        ``(panel, warnings)``. The warnings are the honest part: a BED carries
        no ``period`` and no ``corr_value``, so for any marker STRNaming has no
        reporting range for, the reported number is an uncalibrated repeat count
        rather than a kit allele. Those markers get a ``kit_nomenclature_note``,
        so every call at them raises ``CE_NOMENCLATURE_OFFSET`` — the same
        machinery a curated panel uses to say the same thing.

    Raises:
        PanelError: If the file is unreadable or has no usable rows.
    """
    from frontstr.interp.naming import default_namer

    try:
        text = path.read_text(encoding="utf-8")
    except OSError as exc:
        raise PanelError(f"Cannot read BED {path}: {exc}") from exc

    namer = default_namer()
    systems: list[System] = []
    warnings: list[str] = []
    seen: set[str] = set()

    for lineno, raw in enumerate(text.splitlines(), start=1):
        line = raw.strip()
        if not line or line.startswith(("#", "track", "browser")):
            continue
        fields = line.split("\t") if "\t" in line else line.split()
        if len(fields) < 3:
            raise PanelError(f"{path}:{lineno}: need at least chrom/start/end, got {line!r}")
        chrom, start_s, end_s = fields[0], fields[1], fields[2]
        try:
            start, end = int(start_s), int(end_s)
        except ValueError as exc:
            raise PanelError(f"{path}:{lineno}: non-numeric coordinates in {line!r}") from exc

        motif = ""
        marker = ""
        if len(fields) >= 4:
            if _MOTIF_RE.match(fields[3]):
                motif = fields[3].upper()
                marker = fields[4] if len(fields) >= 5 else ""
            else:
                marker = fields[3]
        if not motif:
            raise PanelError(
                f"{path}:{lineno}: no motif column in {line!r}. FRONTStr bins reads by "
                "repeat-core length, which needs the motif; use "
                "'chrom start end MOTIF name' (e.g. 'chr11 2170987 2171215 AATG TH01')."
            )
        if not marker:
            marker = f"{chrom}_{start}"
        if marker in seen:
            raise PanelError(f"{path}:{lineno}: duplicate marker name {marker!r}")
        seen.add(marker)

        ref_start = start + 1 if coords == "bed0" else start
        if ref_start < 1 or end < ref_start:
            raise PanelError(f"{path}:{lineno}: empty or reversed interval {line!r}")

        has_range = bool(namer and namer.has_range(marker))
        note = None
        if not has_range:
            note = (
                "Region supplied by BED. STRNaming has no reported range for this "
                "marker name, and a BED carries no period or corr_value, so the "
                "number is an uncalibrated repeat count — not a kit CE allele."
            )
            warnings.append(marker)
        systems.append(
            System(
                name=marker,
                chromosome=chrom,
                ref_start=ref_start,
                ref_end=end,
                motif=motif,
                period=-1,
                corr_value=0,
                kit_nomenclature_note=note,
            )
        )

    if not systems:
        raise PanelError(f"{path} has no usable BED rows")
    return Panel(name=name, version="from-bed", systems=systems), warnings


def _format_rows(systems: Iterable[System], coords: BedCoords) -> list[str]:
    rows: list[str] = []
    for s in sorted(systems, key=lambda x: (x.chromosome, x.ref_start)):
        line = _format_row(s, coords)
        if line is not None:
            rows.append(line)
    return rows


def _format_row(s: System, coords: BedCoords) -> str | None:
    if not s.motif:
        return None
    start = s.ref_start - 1 if coords == "bed0" else s.ref_start
    return "\t".join((s.chromosome, str(start), str(s.ref_end), s.motif, s.name))
