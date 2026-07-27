"""ISFG bracketed nomenclature compression.

For each cluster consensus, identify runs of the locus motif(s) and emit a
canonical bracketed notation, e.g. ``[AGAT]12 AGAC [AGAT]3`` for D3S1358.

Implementation is a greedy left-to-right scan taking the longest motif run at
each position — the direct reading of the ISFG bracket convention, written in
Python here.

An earlier docstring described this as a port of toaSTR v1's Perl
``Scripts/Library.pm::compressor``. That was wrong and is corrected here: no
toaSTR source exists in this project or in toaSTR's public repository, which
ships only a Docker Compose file and a SQL schema (verified 2026-07). What the
notation itself owes to prior art is cited in the README.

Phase 1.5 of ROADMAP.md.
"""

from __future__ import annotations

import re

_COMPLEMENT = str.maketrans("ACGTacgt", "TGCAtgca")
_BRACKET_RE = re.compile(r"\[([A-Z]+)\](\d+)")


def _rc(seq: str) -> str:
    """Reverse-complement a nucleotide sequence (case-preserving)."""
    return seq.translate(_COMPLEMENT)[::-1]


def _collapse_single_nt_tokens(tokens: list[str]) -> list[str]:
    """Merge consecutive single-base tokens so we do not emit ``C T C`` spacing.

    ``compress_isfg`` used to ``" ".join`` every output piece; runs of
    non-motif nucleotides became one character per token, which is unreadable
    in HTML tables. Collapse those into contiguous strings like ``CTCC``.
    """
    out: list[str] = []
    buf: list[str] = []
    for tok in tokens:
        if len(tok) == 1 and tok in "ACGTacgt":
            buf.append(tok)
        else:
            if buf:
                out.append("".join(buf))
                buf.clear()
            out.append(tok)
    if buf:
        out.append("".join(buf))
    return out


def motif_repeat_summary(sequence: str, motif: str, *, strand: str = "+") -> str:
    """Human-readable repeat counts for multi-motif STRs (e.g. D3S1358).

    Returns a string like ``TCTAx4 + TCTGx1 + TCTAx3 (TR 64 bp)`` listing each
    uninterrupted motif run left-to-right. This is **not** the commercial-kit
    allele number (that needs kit-specific binning); it is the literal repeat
    structure read from the consensus.

    When no motif run is found, falls back to total length only.
    Pass ``strand="-"`` for markers whose canonical ISFG motif is on the
    reverse complement strand (e.g. D5S818, vWA, CSF1PO).
    """
    from frontstr.interp.stutter import find_motif_runs

    motifs = [m for m in motif.split(",") if m]
    if not sequence:
        return ""
    if not motifs:
        return f"{len(sequence)} bp"
    scan_seq = _rc(sequence) if strand == "-" else sequence
    runs = find_motif_runs(scan_seq, motifs)
    if not runs:
        return f"no motif match, {len(sequence)} bp TR"
    parts = [f"{r.motif}x{r.n_copies}" for r in runs]
    return " + ".join(parts) + f" (TR {len(sequence)} bp)"


def compress_isfg(sequence: str, *, motif: str, strand: str = "+") -> str:
    """Compress ``sequence`` to ISFG bracketed notation using ``motif`` (or motifs).

    Args:
        sequence: Raw nucleotide sequence (uppercase, no whitespace).
        motif: Single motif (e.g. ``"AGAT"``) or comma-separated list
            (e.g. ``"TCTA,TCTG"`` for D3S1358).
        strand: ``"+"`` (default) or ``"-"``. For reverse-strand markers the
            consensus is RC'd before scanning so the output is in the
            canonical ISFG orientation.

    Returns:
        Bracketed display string. Returns the input verbatim if no motif matches.
    """
    motifs = [m for m in motif.split(",") if m] if motif else []
    if not motifs or not sequence:
        return sequence
    if strand == "-":
        sequence = _rc(sequence)

    out: list[str] = []
    n = len(sequence)
    i = 0
    while i < n:
        best_motif: str | None = None
        best_count = 0
        best_run_bp = 0
        for m in motifs:
            k = len(m)
            if k == 0 or i + k > n:
                continue
            if sequence[i : i + k] != m:
                continue
            cnt = 1
            j = i + k
            while j + k <= n and sequence[j : j + k] == m:
                j += k
                cnt += 1
            run_bp = cnt * k
            if run_bp > best_run_bp:
                best_motif, best_count, best_run_bp = m, cnt, run_bp
        if best_motif is None:
            out.append(sequence[i].lower())
            i += 1
        else:
            if best_count >= 2:
                out.append(f"[{best_motif}]{best_count}")
            else:
                out.append(best_motif)
            i += best_run_bp
    out = _collapse_single_nt_tokens(out)
    return " ".join(out)


def ce_from_brackets(isfg: str) -> float | None:
    """Compute CE for compound-motif markers by summing repeat units in the ISFG string.

    Used for ``period == -1`` markers (vWA, FGA, D21S11, etc.) where a direct
    length/period calculation is undefined.

    Counting rules (Phillips 2018):
    - ``[MOTIF]n`` blocks contribute ``n`` units each.
    - Bare uppercase motif tokens (single occurrences, not bracketed) contribute 1.
    - Lowercase tokens (uncounted nucleotides) are ignored.

    Returns None for an empty string or when no motif units are found.
    """
    if not isfg:
        return None
    total = 0
    for m in _BRACKET_RE.finditer(isfg):
        total += int(m.group(2))
    # Bare uppercase runs left after removing bracket tokens = single motif copies
    remaining = _BRACKET_RE.sub(" ", isfg)
    for _ in re.finditer(r"[A-Z]{2,}", remaining):
        total += 1
    return float(total) if total > 0 else None


def ce_from_length(length_bp: int, period: int, corr_value: int) -> float | None:
    """Compute forensic CE allele number from raw TR length.

    Uses the ISFG divmod convention: full repeats are the integer part and
    extra bases are the fractional part (e.g. 9 full + 3 extra = 9.3, not 9.75).

    For multi-motif loci where ``period <= 0``, returns None (CE undefined).
    """
    if period <= 0:
        return None
    q, r = divmod(length_bp - corr_value, period)
    return round(q + r / 10, 1)
