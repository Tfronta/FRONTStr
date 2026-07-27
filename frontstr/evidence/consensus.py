"""Cluster consensus with an explicit, auditable POA backend chain.

Why this module exists
----------------------

The forensic value of FRONTStr is the **sequence** of the called allele: the
ISFG bracket string, the iso-allele catalog match, the microvariant. All three
are computed from ``Cluster.consensus``. If that consensus is a single ONT read
rather than a polished multiple-alignment consensus, it carries that read's
errors — and a mis-called iso-allele is indistinguishable from a correct one in
the output.

Measured on a synthetic 252 bp locus at ONT R10-like error rates (1% sub,
0.3% ins, 0.3% del), consensus vs. the true haplotype:

===========================  ==================  ==================
backend                      3-4 reads in bin    10-16 reads in bin
===========================  ==================  ==================
POA (abPOA / spoa)           0 edits             0 edits
most-common-sequence (mode)  4-6 edits           2-3 edits
===========================  ==================  ==================

The catalog's exact-match window is 2 edits over a ~55 bp repeat core, so the
mode fallback silently degrades iso-allele calling. It is kept only so the
package remains importable without a compiled POA backend — never as a
supported forensic configuration.

Backend chain
-------------

1. ``pyabpoa`` — the ROADMAP's declared backend. Does not build on macOS arm64
   (hardcoded AVX2 x86 intrinsics), so it is tried but not required.
2. ``pyspoa`` — ONT's SPOA bindings. Builds on arm64 and x86.
3. mode — no compiled backend present. Emits
   :class:`ConsensusMethod.MODE`, which the interpretation layer turns into a
   ``CONSENSUS_FALLBACK`` flag on every affected marker.

Alignment mode is **global** (Needleman-Wunsch). Clusters are length-binned
before consensus, so their members are full-window sequences of equal length;
local alignment would be free to trim the flanks and shift the called length,
which directly changes the allele's CE number.

Determinism (chain-of-custody): POA output depends on input order. Cluster
members preserve BAM iteration order, which is coordinate-sorted and stable, so
the same BAM yields the same consensus on every run.
"""

from __future__ import annotations

from collections import Counter
from collections.abc import Callable
from enum import StrEnum
from typing import Any, Protocol


class ConsensusMethod(StrEnum):
    """How a cluster's consensus sequence was derived.

    Serialized into the canonical record so a reviewer can tell a polished
    consensus from an unpolished one. Never repurpose an existing value.
    """

    #: No members — deletion / empty cluster.
    EMPTY = "empty"
    #: Exactly one supporting read; the read verbatim, no polishing possible.
    SINGLE = "single"
    #: Partial-order alignment consensus via ``pyabpoa``.
    POA_ABPOA = "poa_abpoa"
    #: Partial-order alignment consensus via ``pyspoa`` (SPOA).
    POA_SPOA = "poa_spoa"
    #: Most common exact sequence — lossy fallback, no POA backend installed.
    MODE = "mode"

    @property
    def is_poa(self) -> bool:
        """True when a real multiple-alignment consensus was computed."""
        return self in (ConsensusMethod.POA_ABPOA, ConsensusMethod.POA_SPOA)


#: Global alignment (Needleman-Wunsch) in SPOA's ``AlignmentType`` enum.
_SPOA_GLOBAL = 1


class _Backend(Protocol):
    """A POA implementation that turns many sequences into one consensus."""

    method: ConsensusMethod

    def consensus(self, seqs: list[str]) -> str: ...


class _AbpoaBackend:
    method = ConsensusMethod.POA_ABPOA

    def __init__(self, aligner: Any) -> None:
        self._aligner = aligner

    def consensus(self, seqs: list[str]) -> str:
        result = self._aligner.msa(seqs, out_cons=True, out_msa=False)
        cons = getattr(result, "cons_seq", None)
        return str(cons[0]) if cons else ""


class _SpoaBackend:
    method = ConsensusMethod.POA_SPOA

    def __init__(self, poa_fn: Callable[..., tuple[str, list[str]]]) -> None:
        self._poa = poa_fn

    def consensus(self, seqs: list[str]) -> str:
        cons, _msa = self._poa(seqs, genmsa=False, algorithm=_SPOA_GLOBAL)
        return str(cons)


_BACKEND: _Backend | None | bool = False  # False = resolution not yet attempted


def _resolve_backend() -> _Backend | None:
    """Load the first available POA backend, caching the result."""
    global _BACKEND
    if _BACKEND is not False:
        return _BACKEND  # type: ignore[return-value]

    resolved: _Backend | None = None
    try:
        import pyabpoa

        resolved = _AbpoaBackend(pyabpoa.msa_aligner())
    except ImportError:
        try:
            from spoa import poa

            resolved = _SpoaBackend(poa)
        except ImportError:
            resolved = None

    _BACKEND = resolved
    return resolved


def reset_backend_cache() -> None:
    """Forget the resolved backend. For tests that patch the import chain."""
    global _BACKEND
    _BACKEND = False


def poa_backend_name() -> str:
    """Name of the active POA backend, or ``""`` when none is installed.

    Used by ``frontstr doctor`` to surface an unpolished-consensus run before
    it produces results rather than after.
    """
    backend = _resolve_backend()
    return backend.method.value if backend is not None else ""


def build_consensus(seqs: list[str]) -> tuple[str, ConsensusMethod]:
    """Collapse a cluster's member sequences into one consensus.

    Args:
        seqs: Member sequences, in reference orientation. Expected to be of
            equal (or near-equal) length — clusters are length-binned upstream.

    Returns:
        ``(consensus, method)``. ``method`` records how the consensus was
        derived so the report and the canonical JSON can distinguish a polished
        consensus from an unpolished one.
    """
    if not seqs:
        return "", ConsensusMethod.EMPTY
    if len(seqs) == 1:
        return seqs[0], ConsensusMethod.SINGLE

    backend = _resolve_backend()
    if backend is not None:
        try:
            cons = backend.consensus(seqs)
        except Exception:  # a POA failure must degrade, never lose the locus
            cons = ""
        if cons:
            return cons, backend.method

    return Counter(seqs).most_common(1)[0][0], ConsensusMethod.MODE
