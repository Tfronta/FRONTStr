"""Motif-run primitives shared by every layer.

These are pure sequence functions with no forensic opinion, so they sit outside
the Caller / Evidence / Interpretation layering: the Evidence layer needs them
to bin reads by repeat-core length, and the Interpretation layer needs them for
stutter, ISFG compression and catalog lookup. Keeping them here is what lets
Layer 2 use them without importing Layer 3.

Everything is re-exported from :mod:`frontstr.interp.stutter` and
:mod:`frontstr.interp.catalog` so existing import sites keep working.
"""

from __future__ import annotations

from collections.abc import Iterable
from dataclasses import dataclass

#: Non-motif chars between two motif runs above which they belong to different
#: blocks. Matches ``calibrate._GAP_THRESHOLD``: intra-allele spacers (ATG, ta,
#: tccata …) are ≤6 chars, so 12 cleanly separates the repeat core from chance
#: motif hits in the ±100bp flanks of the extraction window.
DEFAULT_GAP_THRESHOLD = 12

_COMPLEMENT = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")


def reverse_complement(seq: str) -> str:
    """Reverse-complement ``seq``, preserving case."""
    return seq.translate(_COMPLEMENT)[::-1]


@dataclass(frozen=True, slots=True)
class MotifRun:
    """One uninterrupted run of a motif inside a sequence."""

    motif: str
    start: int
    n_copies: int

    @property
    def length_bp(self) -> int:
        return self.n_copies * len(self.motif)

    @property
    def end(self) -> int:
        return self.start + self.length_bp


def find_motif_runs(sequence: str, motifs: Iterable[str]) -> list[MotifRun]:
    """Greedy scan for all maximal runs of each motif in ``sequence``.

    Runs of different motifs may not overlap if they share a leading prefix;
    we resolve ties by preferring the longest run starting at each position.
    """
    motif_list = [m for m in motifs if m]
    if not motif_list or not sequence:
        return []
    n = len(sequence)
    runs: list[MotifRun] = []
    i = 0
    while i < n:
        best: MotifRun | None = None
        for m in motif_list:
            k = len(m)
            if i + k > n or sequence[i : i + k] != m:
                continue
            cnt = 1
            j = i + k
            while j + k <= n and sequence[j : j + k] == m:
                cnt += 1
                j += k
            cand = MotifRun(motif=m, start=i, n_copies=cnt)
            if best is None or cand.length_bp > best.length_bp:
                best = cand
        if best is not None:
            runs.append(best)
            i = best.end
        else:
            i += 1
    return runs


def top_two_runs(runs: list[MotifRun]) -> tuple[MotifRun | None, MotifRun | None]:
    """Return (LUS, SLUS) — the two longest runs ordered by copy count then position."""
    if not runs:
        return None, None
    sorted_runs = sorted(runs, key=lambda r: (-r.n_copies, r.start))
    return sorted_runs[0], (sorted_runs[1] if len(sorted_runs) > 1 else None)


def repeat_core_span(
    sequence: str,
    motifs: list[str],
    *,
    gap_threshold: int = DEFAULT_GAP_THRESHOLD,
) -> tuple[int, int] | None:
    """Locate the repeat core inside a flank-inclusive ``sequence``.

    Groups motif runs separated by ≤ ``gap_threshold`` non-motif chars and
    returns the ``(start, end)`` span of the unit-richest group — internal
    spacers preserved, flanks and chance flank motif hits dropped. Returns
    ``None`` when no motif run is present.

    ``sequence`` must already be in canonical (motif) orientation.
    """
    runs = find_motif_runs(sequence, motifs)
    if not runs:
        return None

    groups: list[list[MotifRun]] = [[runs[0]]]
    for run in runs[1:]:
        if run.start - groups[-1][-1].end <= gap_threshold:
            groups[-1].append(run)
        else:
            groups.append([run])

    core = max(groups, key=lambda g: sum(r.length_bp for r in g))
    return core[0].start, core[-1].end


def extract_repeat_core(
    sequence: str,
    motifs: list[str],
    *,
    strand: str = "+",
    gap_threshold: int = DEFAULT_GAP_THRESHOLD,
) -> str:
    """Reduce a flank-inclusive consensus to its repeat-core substring.

    FRONTStr clusters carry the full ±100bp extraction window; catalog entries
    store only the repeat core. Returns the input unchanged if no run is found.
    Reverse-strand markers are canonicalised first.
    """
    if not sequence or not motifs:
        return sequence
    seq = reverse_complement(sequence) if strand == "-" else sequence
    span = repeat_core_span(seq, motifs, gap_threshold=gap_threshold)
    if span is None:
        return seq
    return seq[span[0] : span[1]]


def repeat_core_length(
    sequence: str,
    motifs: list[str],
    *,
    strand: str = "+",
    gap_threshold: int = DEFAULT_GAP_THRESHOLD,
) -> int | None:
    """Length in bp of the repeat core, or ``None`` when no motif run is found.

    This is the quantity the Evidence layer bins reads by: it is the part of
    the window that actually carries allelic information, so an indel error in
    the ±100bp flanks — the large majority of them — no longer splits a real
    allele across two length bins. It stays exact within the array, so a
    microvariant is preserved: ``[AATG]9`` (36 bp) and ``[AATG]6 ATG [AATG]3``
    (39 bp) are both 9 repeat units but land in different bins, which is what
    keeps TH01 9 and 9.3 distinguishable.
    """
    if not sequence or not motifs:
        return None
    seq = reverse_complement(sequence) if strand == "-" else sequence
    span = repeat_core_span(seq, motifs, gap_threshold=gap_threshold)
    if span is None:
        return None
    return span[1] - span[0]
