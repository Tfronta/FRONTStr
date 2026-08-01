"""Per-locus narrative trace: what the caller did, in order, and why.

The process log (:mod:`frontstr.log`) answers *what was concluded*: one line per
marker with the genotype and the coverage. That is not enough to validate a
call. A reviewer looking at "D3S1358 homozygous 14, 35 reads" cannot see where
those 35 reads came from, what happened to the reads that are missing, how many
candidate alleles were considered, or why the ones that were discarded were
discarded.

This module produces the other thing: a step-by-step account of one locus,
written to be read by a person. It follows the pipeline's actual order —

    fetch → filter → bin → cluster → consensus → name → stutter → classify
    → suppress → call

— and at each step it reports the numbers that step acted on. The design goal
is that someone who does not know the code can follow a genotype back to the
reads, and that someone who does know the code can spot the stage where a
surprising call went wrong.

Structured first, prose second
------------------------------

:class:`LocusTrace` is a plain data record; :func:`render_locus` turns it into
text. Keeping those apart means the same trace can later be serialized into the
audit record or rendered into the HTML report without reformatting strings.

Prior art: HipSTR prints a comparable per-region narrative. The stages here are
FRONTStr's own — there is no haplotype generation or PCR-duplicate removal in a
long-read pileup — but the idea of showing the read funnel with named rejection
reasons is taken from it, and it is the single most useful part.
"""

from __future__ import annotations

import time
from collections.abc import Iterator
from contextlib import contextmanager
from dataclasses import dataclass, field

from frontstr.evidence.pileup import PileupCounts

#: Width of the label column in the rendered output.
_LABEL = 42

#: bp of flank shown on each side of the repeat core. Enough to recognise the
#: anchor sequence and to see a flank variant, short enough that the cores stay
#: aligned in one column.
FLANK_SHOWN = 12

#: Cap on the width the core column is padded to. Right flanks are aligned so
#: a length difference is visible at a glance, but padding every row to the
#: longest core lets one corrupt cluster set the width for the whole locus —
#: a single 200 bp noise core pushed every DYS391 row out to 271 characters.
#: Beyond the cap, long cores simply overflow and only they lose the alignment.
_CORE_ALIGN_CAP = 96


@dataclass(slots=True)
class BinTrace:
    """One repeat-core-length bin and how many reads landed in it."""

    key_bp: int
    n_reads: int
    #: False when the read's repeat core could not be located and the raw
    #: window length was used instead — worth showing, because such a read is
    #: binned on a number that includes flank error.
    from_core: bool = True


@dataclass(slots=True)
class ClusterTrace:
    """One candidate allele, as it looked at the end of the pipeline."""

    index: int
    n_reads: int
    fraction: float
    length_bp: int
    consensus_method: str
    n_hp1: int
    n_hp2: int
    n_untagged: int
    number_label: str
    number_method: str
    naming_status: str
    strnaming_name: str
    isfg: str
    expected_stutter: float
    status: str
    n_reads_absorbed: int
    hp_rescued: bool
    called: bool
    #: The consensus split for display, in **canonical (motif) orientation** so
    #: it reads the same way as the ISFG string. The core is kept whole — it is
    #: what distinguishes one allele from another — while the flanks are cut to
    #: :data:`FLANK_SHOWN` bp, because ~80% of a panel window is flank and
    #: printing all of it would bury the part that matters.
    #: Reads by mapped strand. A real allele is drawn from both strands roughly
    #: in proportion to the locus; a strand-specific basecalling artefact is
    #: not, because ONT reads the two strands through complementary k-mer
    #: contexts and its systematic errors follow that. The statistical test
    #: needs ``strand_bias_min_reads`` before its binomial has any power, so on
    #: the small candidates — exactly where artefacts live — these numbers are
    #: the only evidence there is.
    n_forward: int = 0
    n_reverse: int = 0

    flank_left: str = ""
    core: str = ""
    flank_right: str = ""
    #: False when no motif run was found and the whole consensus is shown raw.
    core_found: bool = True
    #: Phase block (BAM ``PS``) the haplotype-tagged reads came from, and how
    #: many distinct blocks they span. More than one means the ``HP`` counts on
    #: this row are *not* haplotype evidence — see
    #: :mod:`frontstr.interp.haplotype`.
    phase_set: int | None = None
    n_phase_sets: int = 0


@dataclass(slots=True)
class LocusTiming:
    """Wall-clock seconds per pipeline stage, for one locus.

    Not performance tuning — diagnosis. A locus that takes twenty times its
    neighbours is telling you something about the data: a pileup that dominates
    means reads are being fetched and thrown away, and a consensus that
    dominates means POA is aligning many divergent candidates, which is what a
    noisy locus looks like from the inside. Both are visible in the numbers
    before they are visible in the genotype.

    Stages match the order :func:`render_locus` narrates them in, so a slow
    number can be read against the section that produced it. Only populated
    when the run is tracing; untraced runs never construct this.
    """

    #: Fetch reads and apply the read filters.
    pileup: float = 0.0
    #: Group by repeat-core length, then split each bin by pairwise identity.
    clustering: float = 0.0
    #: POA consensus per cluster. Measured inside clustering and subtracted
    #: from it, so the two add up rather than double-counting.
    consensus: float = 0.0
    #: Derive each candidate's allele number and bracket string (STRNaming).
    naming: float = 0.0
    #: Build the expected-stutter table and classify every candidate against it.
    stutter: float = 0.0
    #: Haplotype-aware phantom suppression and catalog annotation.
    suppression: float = 0.0
    #: Apply the calling rules and settle on a genotype.
    calling: float = 0.0

    @property
    def total(self) -> float:
        """Everything above. Stages are disjoint, so this is their sum."""
        return (
            self.pileup
            + self.clustering
            + self.consensus
            + self.naming
            + self.stutter
            + self.suppression
            + self.calling
        )

    def stages(self) -> list[tuple[str, float]]:
        """``(label, seconds)`` in pipeline order, for rendering."""
        return [
            ("Pileup (fetch + filter)", self.pileup),
            ("Clustering (bin + identity)", self.clustering),
            ("Consensus (POA)", self.consensus),
            ("Allele naming", self.naming),
            ("Stutter + classification", self.stutter),
            ("Phantom suppression + catalog", self.suppression),
            ("Genotype call", self.calling),
        ]

    def add(self, other: LocusTiming) -> None:
        """Accumulate ``other`` into this one, for the run total."""
        self.pileup += other.pileup
        self.clustering += other.clustering
        self.consensus += other.consensus
        self.naming += other.naming
        self.stutter += other.stutter
        self.suppression += other.suppression
        self.calling += other.calling


@contextmanager
def timed(timing: LocusTiming | None, stage: str) -> Iterator[None]:
    """Add the block's wall-clock seconds to ``timing.<stage>``.

    A no-op when ``timing`` is ``None``, which is every untraced run, so the
    hot path pays one attribute check rather than two clock reads.
    """
    if timing is None:
        yield
        return
    started = time.perf_counter()
    try:
        yield
    finally:
        setattr(timing, stage, getattr(timing, stage) + time.perf_counter() - started)


@dataclass(slots=True)
class LocusTrace:
    """Everything one locus did, in pipeline order."""

    marker: str
    chrom: str
    start: int  # 1-based inclusive, as the panel writes it
    end: int
    motif: str
    period: int
    strand: str

    min_mapq: int = 20
    flank_anchor: int = 20
    identity_threshold: float = 0.97
    analytical_thresh: float = 0.02
    calling_thresh: float = 0.10

    counts: PileupCounts | None = None
    binned_on_core: bool = True
    bins: list[BinTrace] = field(default_factory=list)
    consensus_backend: str = ""
    #: Fragments folded back together because POA gave them the same consensus.
    #: Without this the reader sees fewer clusters than bins and no reason why.
    merged_on_consensus: int = 0
    clusters: list[ClusterTrace] = field(default_factory=list)

    #: Set for markers that bypass the STR pipeline entirely (amelogenin).
    #: Rendered in place of the funnel/binning/clustering sections.
    note: str = ""
    n_suppressed_phantoms: int = 0
    #: Coverage share of the strongest called allele; ``None`` unless the call
    #: is a heterozygote. See :attr:`MarkerResult.allele_balance`.
    allele_balance: float | None = None
    balanced_ab_max: float = 0.65
    call_rule: str = ""
    tri_type: str = ""
    #: Wall-clock per stage. Set for every traced locus.
    timing: LocusTiming | None = None
    called_labels: list[str] = field(default_factory=list)
    flags: list[tuple[str, str]] = field(default_factory=list)  # (severity, code)

    @property
    def window_bp(self) -> int:
        return self.end - self.start + 1


@dataclass(slots=True)
class RunHeader:
    """What the run was given, before it does anything with it.

    Printed first so a log — especially an unattended benchmark over hundreds of
    samples — states its own inputs. A trace whose provenance has to be
    reconstructed from the invocation is not much use months later.
    """

    inputs: list[str] = field(default_factory=list)
    panel_name: str = ""
    panel_version: str = ""
    n_markers: int = 0
    min_mapq: int = 20
    flank_anchor: int = 20
    identity_threshold: float = 0.97
    analytical_thresh: float = 0.02
    calling_thresh: float = 0.10
    consensus_backend: str = ""
    naming_markers: int = 0
    #: Version and provenance of the stutter model in force. Named in the
    #: header because it decides which candidates become artefacts: a run
    #: whose trace does not say which model it used cannot be reproduced.
    stutter_model: str = ""
    tool_version: str = ""
    #: ``(name, value, default, provenance)`` for every overridden parameter.
    overrides: list[tuple[str, object, object, str]] = field(default_factory=list)


def render_header(h: RunHeader) -> str:
    """Render the run header: inputs, panel, and every threshold in force."""
    lines = [f"── FRONTStr {h.tool_version}".rstrip()]
    kinds: dict[str, int] = {}
    for path in h.inputs:
        kind = "CRAM" if path.lower().endswith(".cram") else "BAM"
        kinds[kind] = kinds.get(kind, 0) + 1
    detected = ", ".join(f"{n} {kind}" for kind, n in sorted(kinds.items())) or "none"
    lines.append(_row("Detected", f"{detected}"))
    for path in h.inputs:
        lines.append(_row(path, "", indent=6).rstrip())
    lines.append(_row("Panel", f"{h.panel_name} {h.panel_version}".strip()))
    lines.append(_row("Markers in panel", h.n_markers))
    lines.append(
        _row("Read filters", f"MAPQ >= {h.min_mapq}, {h.flank_anchor} bp clean flank each side")
    )
    lines.append(
        _row(
            "Thresholds",
            f"analytical {h.analytical_thresh:.0%}, calling {h.calling_thresh:.0%}, "
            f"cluster identity {h.identity_threshold:.2f}",
        )
    )
    if h.consensus_backend:
        lines.append(_row("Consensus backend", h.consensus_backend))
    if h.stutter_model:
        lines.append(_row("Stutter model", h.stutter_model))
    lines.append(
        _row(
            "Allele naming",
            f"STRNaming, offline slice cache, {h.naming_markers} markers"
            if h.naming_markers
            else "legacy CE arithmetic (STRNaming unavailable)",
        )
    )
    if h.overrides:
        lines.append(
            _row(
                "Overridden",
                ", ".join(f"{n}={v} (default {d})" for n, v, d, _ in h.overrides),
            )
        )
        derived = [n for n, _, _, prov in h.overrides if prov == "derived"]
        if derived:
            lines.append(
                _row(
                    "NOT a default run",
                    f"{', '.join(derived)} — measured default(s) overridden",
                )
            )
    lines.append("")
    return "\n".join(lines)


def _row(label: str, value: object, indent: int = 2) -> str:
    """One ``label ....... value`` line, column-aligned where the label allows.

    ``ljust`` alone is not enough: a label longer than the column runs straight
    into its value ("...left flank anchor2"), so a minimum gap is enforced.
    """
    pad = " " * indent
    text = f"{pad}{label}"
    gap = max(_LABEL - len(text), 2) if str(value) else 0
    return f"{text}{' ' * gap}{value}".rstrip()


def _plural(n: int, word: str) -> str:
    return f"{n} {word}{'' if n == 1 else 's'}"


def render_locus(t: LocusTrace) -> str:
    """Render one :class:`LocusTrace` as an indented, readable narrative."""
    lines: list[str] = []
    add = lines.append

    period = "compound" if t.period == -1 else f"period {t.period}"
    add(
        f"── {t.marker}  {t.chrom}:{t.start}-{t.end}  "
        f"({t.window_bp} bp window, motif {t.motif}, {period}"
        + (", minus strand" if t.strand == "-" else "")
        + ")"
    )

    if t.note:
        add(_row(t.note, ""))
        add(_row("Genotype", f"{_genotype(t)}   [{t.call_rule}]{_qc_suffix(t)}"))
        return "\n".join(lines)

    # 1. Read funnel ---------------------------------------------------------
    if t.counts is not None:
        funnel = t.counts
        add(_row("Reads fetched around the window", funnel.fetched))
        # Every reason, including the ones that did not fire. A zero is a
        # positive statement — this was checked and found nothing — and it is
        # what separates an auditable funnel from a plausible one.
        add(_row(f"Rejected ({funnel.n_rejected})", "", indent=2).rstrip())
        for reason, n in funnel.all_reasons():
            add(_row(reason.value, n, indent=6))
        add(_row("Spanning the whole window", f"{funnel.kept}   (total locus coverage)"))
        if funnel.kept == 0:
            add(_row("No usable reads", "locus reported as no_data"))
            add(
                _row(
                    "",
                    f"nothing cleared MAPQ {t.min_mapq} with {t.flank_anchor} bp of clean "
                    "flank each side; --min-mapq and --flank-anchor are the knobs",
                    indent=6,
                ).rstrip()
            )
            return "\n".join(lines)
    add("")

    # 2. Binning -------------------------------------------------------------
    if t.bins:
        basis = (
            "repeat-core length (flank indel errors ignored)"
            if t.binned_on_core
            else "raw window length (no motif configured)"
        )
        add(_row("Step 1: grouped by length", f"{_plural(len(t.bins), 'bin')}, using {basis}"))
        fallbacks = sum(b.n_reads for b in t.bins if not b.from_core)
        for b in sorted(t.bins, key=lambda b: -b.n_reads):
            note = "" if b.from_core else "   (core not locatable, window length used)"
            add(_row(f"{b.key_bp} bp core", f"{_plural(b.n_reads, 'read')}{note}", indent=6))
        if fallbacks:
            add(_row("Reads binned on window length", fallbacks, indent=6))

    # 3. Clustering + consensus ---------------------------------------------
    if t.clusters:
        split = len(t.clusters) - len(t.bins)
        outcome = (
            f"{_plural(len(t.clusters), 'cluster')}, none split"
            if split <= 0
            else f"{_plural(len(t.clusters), 'cluster')}, {split} more than bins"
        )
        add(
            _row(
                "Step 2: split by sequence",
                f"{outcome} (identity below {t.identity_threshold:.2f} separates)",
            )
        )
        if t.consensus_backend:
            add(_row("Step 3: consensus per cluster", t.consensus_backend))
        if t.merged_on_consensus:
            add(
                _row(
                    "Step 4: merged on consensus",
                    f"{_plural(t.merged_on_consensus, 'fragment')} folded back in — "
                    "same allele, split by read-to-read identity, reunited by POA",
                )
            )
        add("")

        # 4. Per-cluster detail ---------------------------------------------
        add(_row("Candidates, strongest first", "reads shown are per-allele coverage"))
        for c in t.clusters:
            mark = "*" if c.called else " "
            hp = f"HP1 {c.n_hp1} / HP2 {c.n_hp2} / untagged {c.n_untagged}"
            if c.n_phase_sets > 1:
                hp += f"  [{c.n_phase_sets} phase blocks — HP not comparable]"
            elif c.phase_set is not None:
                hp += f"  [block {c.phase_set}]"
            add(
                f"    {mark} #{c.index}  {_plural(c.n_reads, 'read'):>10}  "
                f"{c.fraction:5.1%}  {c.length_bp:>4} bp  "
                f"{c.n_forward}+/{c.n_reverse}-  {hp}"
            )
            name = c.strnaming_name or c.isfg or "—"
            if len(name) > 78:
                name = name[:75] + "…"
            add(f"        name      {name}")
            via = (
                "STRNaming"
                if c.number_method == "strnaming"
                else f"{c.number_method} (STRNaming: {c.naming_status})"
            )
            add(f"        number    {c.number_label or '—'}   via {via}")
            if c.consensus_method not in ("poa_spoa", "poa_abpoa"):
                add(f"        consensus {c.consensus_method}   NOT polished by POA")
            why = _explain_status(c, t)
            add(f"        verdict   {c.status}{why}")
            if c.n_reads_absorbed:
                add(f"        absorbed  {c.n_reads_absorbed} read(s) from same-haplotype splits")
        add("")
        lines.extend(_render_sequences(t))

    if t.n_suppressed_phantoms:
        add(
            _row(
                "Same-haplotype phantoms suppressed",
                f"{t.n_suppressed_phantoms} (one allele per haplotype)",
            )
        )

    # 5. The call ------------------------------------------------------------
    add(_row("Genotype", f"{_genotype(t)}   [{t.call_rule}]{_qc_suffix(t)}"))
    if t.counts is not None:
        add(_row("Coverage", _coverage_line(t)))
    if t.allele_balance is not None:
        add(_row("Allele balance", _balance_line(t)))
    if t.tri_type:
        add(_row("Triallelic pattern", t.tri_type))
    for severity, code in t.flags:
        add(_row(f"Flag ({severity})", code, indent=6))

    # 6. What it cost ---------------------------------------------------------
    if t.timing is not None:
        add("")
        for line in _render_timing(t.timing, title="Locus timing"):
            add(line)
    return "\n".join(lines)


def _qc_suffix(t: LocusTrace) -> str:
    """The QC verdict appended to the genotype: the flags themselves, or nothing.

    Deliberately **not** an aggregated ``PASS``. A single green label that
    stands for several checks trains a reviewer to stop reading the individual
    ones, and a label that shows on 95% of loci stops carrying information at
    all. So a clean locus says nothing and a flagged one names what fired.
    """
    if not t.flags:
        return ""
    worst = (
        "error"
        if any(s == "error" for s, _ in t.flags)
        else ("warn" if any(s == "warn" for s, _ in t.flags) else "info")
    )
    return f"   {worst.upper()}: " + ", ".join(code for _, code in t.flags)


def _balance_line(t: LocusTrace) -> str:
    """Allele balance with the scale spelled out, because 0.51 alone is opaque."""
    ab = t.allele_balance
    assert ab is not None
    verdict = "balanced" if ab <= t.balanced_ab_max else "uneven"
    return (
        f"{ab:.2f}  ({verdict}; 0.50 is even, "
        f"balanced up to {t.balanced_ab_max:.2f}, "
        "strongest allele over the called pair)"
    )


def _genotype(t: LocusTrace) -> str:
    """The called genotype with **per-allele** coverage attached to each allele.

    Per-allele read counts are the headline claim of an integer-coverage caller,
    and they were only visible partway up the trace among the rejected
    candidates. A reviewer reading the conclusion should not have to scroll back
    to learn that a "9.3, 7" heterozygote rests on 10 reads and 7 reads.
    """
    called = [c for c in t.clusters if c.called]
    if not called:
        return ", ".join(t.called_labels) if t.called_labels else "none"
    return ", ".join(f"{c.number_label} ({_plural(c.n_reads, 'read')})" for c in called)


def _coverage_line(t: LocusTrace) -> str:
    """Total, and how it splits across called alleles and everything discarded."""
    assert t.counts is not None
    total = t.counts.kept
    called = [c for c in t.clusters if c.called]
    on_alleles = sum(c.n_reads for c in called)
    parts = [f"{total} at the locus", f"{on_alleles} on called allele(s)"]
    if total > on_alleles:
        parts.append(f"{total - on_alleles} on discarded candidates")
    if any(c.n_hp1 or c.n_hp2 for c in called):
        hp = " / ".join(f"{c.number_label}: HP1 {c.n_hp1} HP2 {c.n_hp2}" for c in called)
        parts.append(f"phased {hp}")
    return "; ".join(parts)


def _display_core(c: ClusterTrace) -> str:
    """The core as printed. **A called allele is never truncated.**

    That is the rule that matters: the point of showing bases is to validate the
    genotype, so the sequences behind it print in full however long they are —
    D21S11's core runs past 180 bp and still does.

    Uncalled candidates are cut at the alignment cap. A corrupt ONT read either
    has no locatable core (the whole 250 bp window lands in ``core``) or has one
    that starts at position 0 because a chance motif hit in the flank merged
    into the array. Either way it is noise, printing it whole crowds out the
    candidates that matter, and nobody validates a call by reading a rejected
    read end to end.
    """
    if c.called or len(c.core) <= _CORE_ALIGN_CAP:
        return c.core
    return c.core[: _CORE_ALIGN_CAP - 1] + "…"


def _render_sequences(t: LocusTrace) -> list[str]:
    """The bases themselves: flank … repeat core … flank, cores in one column.

    Aligned deliberately. Read down the core column and a 4 bp length step, an
    interrupted motif or a single substitution is visible without counting
    characters — which is the whole reason a sequence-resolved caller is worth
    having over a length-based one.
    """
    shown = [c for c in t.clusters if c.core or c.flank_left or c.flank_right]
    if not shown:
        return []
    cores = {c.index: _display_core(c) for c in shown}
    left_w = max(len(c.flank_left) for c in shown)
    core_w = min(max(len(x) for x in cores.values()), _CORE_ALIGN_CAP)

    out = [_row("Sequences (flank … repeat core … flank):", "").rstrip()]
    if any(not c.core_found for c in shown):
        out.append(
            _row("rows with no flank: no motif run found, raw consensus", "", indent=6).rstrip()
        )
    if any(len(_display_core(c)) < len(c.core) for c in shown):
        out.append(
            _row(
                f"uncalled candidates cut at {_CORE_ALIGN_CAP} bp (…); called ones never are",
                "",
                indent=6,
            ).rstrip()
        )
    for c in shown:
        mark = "*" if c.called else " "
        row = (
            f"    {mark} #{c.index}  {c.flank_left:>{left_w}}  "
            f"{cores[c.index]:<{core_w}}  {c.flank_right}"
        )
        out.append(row.rstrip())
    out.append("")
    return out


def _explain_status(c: ClusterTrace, t: LocusTrace) -> str:
    """Say *why* a candidate got its status, in the reviewer's terms."""
    if c.status == "noise":
        return f"  (below the {t.analytical_thresh:.0%} analytical threshold)"
    if c.status == "artefact":
        return f"  (above analytical, below the {t.calling_thresh:.0%} calling threshold)"
    if c.status == "stutter":
        return f"  (≤ {c.expected_stutter:.1f} reads expected as stutter from a stronger allele)"
    if c.status == "hp_phantom":
        return "  (same haplotype as a stronger, near-identical candidate)"
    if c.status == "allele" and c.hp_rescued:
        return "  (called on phasing: opposite haplotype, despite the read ratio)"
    if c.status == "allele" and not c.called:
        return "  (real candidate, not in the reported genotype)"
    return ""


def _render_timing(timing: LocusTiming, *, title: str, indent: int = 2) -> list[str]:
    """The stage breakdown, with each stage's share of the total.

    The share is the point. Absolute seconds say a locus was slow; the share
    says *which stage* made it slow, which is the thing you act on. Stages that
    took no measurable time are still listed, for the same reason the read
    funnel lists its zeros.
    """
    total = timing.total
    lines = [_row(title, f"{total:.3f} s", indent=indent).rstrip()]
    for label, seconds in timing.stages():
        share = f"{seconds / total:5.1%}" if total > 0 else "    —"
        lines.append(_row(label, f"{seconds:8.4f} s  {share}", indent=indent + 4))
    return lines


def render_run_summary(traces: list[LocusTrace]) -> str:
    """Closing totals across a run, so the tail of the log is not just the last locus."""
    if not traces:
        return "No loci processed."
    fetched = sum(t.counts.fetched for t in traces if t.counts)
    kept = sum(t.counts.kept for t in traces if t.counts)
    called = sum(1 for t in traces if t.called_labels)
    no_data = [t.marker for t in traces if not t.called_labels]
    fallbacks = [
        (t.marker, c.naming_status)
        for t in traces
        for c in t.clusters
        if c.called and c.number_method != "strnaming"
    ]

    lines = ["", "── Run summary"]
    lines.append(_row("Loci processed", len(traces)))
    lines.append(_row("Loci with a genotype", f"{called}/{len(traces)}"))
    if fetched:
        lines.append(
            _row(
                "Reads fetched, then used",
                f"{fetched} fetched, {kept} spanned ({kept / fetched:.1%})",
            )
        )
    if no_data:
        lines.append(_row("No genotype", ", ".join(no_data)))
    if fallbacks:
        lines.append(_row("Called on the legacy CE path", len(fallbacks)))
        for marker, status in fallbacks:
            lines.append(_row(marker, status, indent=6))
    flags: dict[str, int] = {}
    for t in traces:
        for _severity, code in t.flags:
            flags[code] = flags.get(code, 0) + 1
    if flags:
        lines.append(_row("Flags raised", ""))
        for code, n in sorted(flags.items(), key=lambda kv: (-kv[1], kv[0])):
            lines.append(_row(code, n, indent=6))

    timed_loci = [t for t in traces if t.timing is not None]
    if timed_loci:
        total = LocusTiming()
        for t in timed_loci:
            assert t.timing is not None
            total.add(t.timing)
        lines.append("")
        # "summed over loci", not "whole run": per-run setup — building the
        # STRNaming reference structures, deriving the run-level QC flags — sits
        # outside every locus and is deliberately not counted here. Untimed loci
        # are named rather than quietly dropped from the denominator; AMEL is
        # sex typing, not the STR path, so it has no stages to attribute.
        untimed = len(traces) - len(timed_loci)
        title = "Timing, summed over loci"
        if untimed:
            title += f" ({len(timed_loci)} of {len(traces)}; {untimed} not on the STR path)"
        lines.extend(_render_timing(total, title=title))
        # The slowest locus by name. An aggregate hides the one pathological
        # locus inside twenty-four ordinary ones, and that locus is exactly the
        # one worth opening.
        slowest = max(timed_loci, key=lambda t: t.timing.total if t.timing else 0.0)
        if slowest.timing is not None and slowest.timing.total > 0:
            share = slowest.timing.total / total.total if total.total else 0.0
            lines.append(
                _row(
                    "Slowest locus",
                    f"{slowest.marker}  {slowest.timing.total:.3f} s ({share:.0%} of the run)",
                )
            )
    return "\n".join(lines)
