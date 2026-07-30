"""Per-cluster forensic classification.

Implements the two-rule classifier:

1. **Coverage fraction** vs analytical / calling thresholds (handles mixtures
   and low-input ONT data).
2. **Sequence-based expected stutter** from :func:`build_expected_stutter`
   (handles parent ± 1 / ± 2 stutter peaks).

Both rules must agree before a cluster is promoted to ``allele``.
"""

from __future__ import annotations

from frontstr.interp.models import Allele, AlleleStatus


def classify_allele(
    allele: Allele,
    *,
    total_reads: int,
    expected_stutter: dict[str, float],
    analytical_thresh: float,
    calling_thresh: float,
) -> AlleleStatus:
    """Classify a single :class:`Allele` candidate.

    Args:
        allele: The candidate, with raw ``n_reads_total``, ``consensus`` and
            ``is_deletion`` already populated.
        total_reads: Coverage at this locus (sum of all cluster ``n_reads``).
        expected_stutter: Output of
            :func:`frontstr.interp.stutter.build_expected_stutter`.
        analytical_thresh: Below this fraction of total reads, the cluster is
            ``noise`` (e.g. ``0.02`` = 2%).
        calling_thresh: Below this fraction, the cluster is ``artefact``
            (above analytical but below calling) — e.g. ``0.10``.

    Returns:
        The :class:`AlleleStatus` for this cluster.
    """
    if total_reads <= 0:
        return AlleleStatus.NO_DATA
    if allele.is_deletion:
        return AlleleStatus.DELETION

    es = expected_stutter.get(allele.consensus, 0.0)
    if es > 0 and allele.n_reads_total <= es:
        return AlleleStatus.STUTTER

    frac = allele.n_reads_total / total_reads
    if frac < analytical_thresh:
        return AlleleStatus.NOISE
    if frac < calling_thresh:
        return AlleleStatus.ARTEFACT

    return AlleleStatus.ALLELE
