"""ISFG bracketed nomenclature compression.

Ported from toaSTR v1 ``Scripts/Library.pm::compressor`` (Perl).

For each cluster consensus, identify runs of the locus motif(s) and emit a
canonical bracketed notation, e.g. ``[AGAT]12 AGAC [AGAT]3`` for D3S1358.

Phase 1.5 of ROADMAP.md. See research-toastr.md §15 and §8 step 4.
"""

from __future__ import annotations


def _collapse_single_nt_tokens(tokens: list[str]) -> list[str]:
    """Merge consecutive single-base tokens so we do not emit ``C T C`` spacing.

    ``compress_isfg`` used to ``" ".join`` every output piece; runs of
    non-motif nucleotides became one character per token, which is unreadable
    in HTML tables. Collapse those into contiguous strings like ``CTCC``.
    """
    out: list[str] = []
    buf: list[str] = []
    for tok in tokens:
        if len(tok) == 1 and tok in "ACGT":
            buf.append(tok)
        else:
            if buf:
                out.append("".join(buf))
                buf.clear()
            out.append(tok)
    if buf:
        out.append("".join(buf))
    return out


def motif_repeat_summary(sequence: str, motif: str) -> str:
    """Human-readable repeat counts for multi-motif STRs (e.g. D3S1358).

    Returns a string like ``TCTAx4 + TCTGx1 + TCTAx3 (TR 64 bp)`` listing each
    uninterrupted motif run left-to-right. This is **not** the commercial-kit
    allele number (that needs kit-specific binning); it is the literal repeat
    structure read from the consensus.

    When no motif run is found, falls back to total length only.
    """
    from frontstr.interp.stutter import find_motif_runs

    motifs = [m for m in motif.split(",") if m]
    if not sequence:
        return ""
    if not motifs:
        return f"{len(sequence)} bp"
    runs = find_motif_runs(sequence, motifs)
    if not runs:
        return f"no motif match, {len(sequence)} bp TR"
    parts = [f"{r.motif}x{r.n_copies}" for r in runs]
    return " + ".join(parts) + f" (TR {len(sequence)} bp)"


def compress_isfg(sequence: str, *, motif: str) -> str:
    """Compress ``sequence`` to ISFG bracketed notation using ``motif`` (or motifs).

    Args:
        sequence: Raw nucleotide sequence (uppercase, no whitespace).
        motif: Single motif (e.g. ``"AGAT"``) or comma-separated list
            (e.g. ``"TCTA,TCTG"`` for D3S1358).

    Returns:
        Bracketed display string. Returns the input verbatim if no motif matches.
    """
    motifs = [m for m in motif.split(",") if m] if motif else []
    if not motifs or not sequence:
        return sequence

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
            out.append(sequence[i])
            i += 1
        else:
            if best_count >= 2:
                out.append(f"[{best_motif}]{best_count}")
            else:
                out.append(best_motif)
            i += best_run_bp
    out = _collapse_single_nt_tokens(out)
    return " ".join(out)


def ce_from_length(length_bp: int, period: int, corr_value: int) -> float | None:
    """Compute forensic CE allele number from raw TR length.

    ``CE = (length - corr_value) / period``

    For multi-motif loci where ``period <= 0``, returns None (CE undefined).
    """
    if period <= 0:
        return None
    return round((length_bp - corr_value) / period, 1)
