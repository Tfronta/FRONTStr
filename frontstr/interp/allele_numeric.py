"""Forensic allele number for reports (CE + LongTR-style offset paths).

When the panel declares a fixed ``period >= 1``, we use the standard
``(length - corr) / period`` CE (same as :func:`frontstr.interp.isfg.ce_from_length`).

For compound markers (``period == -1``) we follow the same *idea* as LongTR:
``Δbp`` versus a reference anchor sequence length, converted with
``allele_bp_step`` (default 4 bp per +1 repeat unit), then added to an optional
``reference_ce`` — the kit-style allele number assigned to the **reference
genome** haplotype in that panel window.

If ``reference_ce`` is unset, we still emit **Δ / step** as a provisional
number (labelled ``delta_only`` in JSON) so the caller always has a numeric
axis; labs should set ``reference_ce`` from STRBase / kit tables for CODIS-scale
output.
"""

from __future__ import annotations

from frontstr.interp.isfg import ce_from_length
from frontstr.panel.models import System


def resolve_ref_anchor_bp(
    system: System,
    *,
    explicit: int | None,
    longtr_ref_len: int | None,
) -> int:
    """Reference TR length (bp) for Δ math: explicit > LongTR REF > panel span."""
    if explicit is not None:
        return explicit
    if longtr_ref_len is not None:
        return longtr_ref_len
    return system.span()


def compute_allele_numeric(
    allele_len_bp: int,
    system: System,
    ref_anchor_bp: int | None,
) -> tuple[float | None, str]:
    """Return ``(forensic_number, source)``.

    ``source`` is one of:
    - ``period_ce`` — standard CE from repeat length and panel period/corr
    - ``reference_offset`` — ``reference_ce + (len - ref) / allele_bp_step``
    - ``delta_only`` — same offset math without an absolute ``reference_ce``
    - ``deletion`` / ``unavailable`` — no meaningful number
    """
    if allele_len_bp <= 0:
        return None, "deletion"
    if ref_anchor_bp is None:
        return None, "unavailable"

    ce = ce_from_length(allele_len_bp, system.period, system.corr_value)
    if ce is not None:
        return float(ce), "period_ce"

    step = float(max(1, system.allele_bp_step))
    bp_delta = float(allele_len_bp - ref_anchor_bp)
    units = bp_delta / step
    if system.reference_ce is not None:
        return round(float(system.reference_ce) + units, 2), "reference_offset"
    return round(units, 2), "delta_only"
