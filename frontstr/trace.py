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
    flank_left: str = ""
    core: str = ""
    flank_right: str = ""
    #: False when no motif run was found and the whole consensus is shown raw.
    core_found: bool = True


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
    clusters: list[ClusterTrace] = field(default_factory=list)

    #: Set for markers that bypass the STR pipeline entirely (amelogenin).
    #: Rendered in place of the funnel/binning/clustering sections.
    note: str = ""
    n_suppressed_phantoms: int = 0
    call_rule: str = ""
    tri_type: str = ""
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
    tool_version: str = ""


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
    lines.append(
        _row(
            "Allele naming",
            f"STRNaming, offline slice cache, {h.naming_markers} markers"
            if h.naming_markers
            else "legacy CE arithmetic (STRNaming unavailable)",
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
        add(_row("Genotype", f"{_genotype(t)}   [{t.call_rule}]"))
        return "\n".join(lines)

    # 1. Read funnel ---------------------------------------------------------
    if t.counts is not None:
        funnel = t.counts
        add(_row("Reads fetched around the window", funnel.fetched))
        if funnel.n_rejected:
            add(_row(f"Rejected ({funnel.n_rejected})", "", indent=2).rstrip())
            for reason, n in funnel.reasons():
                add(_row(reason.value, n, indent=6))
        add(_row("Spanning the whole window", f"{funnel.kept}   (total locus coverage)"))
        if funnel.kept == 0:
            add(_row("No usable reads", "locus reported as no_data"))
            return "\n".join(lines)
    add("")

    # 2. Binning -------------------------------------------------------------
    if t.bins:
        basis = (
            "repeat-core length (flank indel errors ignored)"
            if t.binned_on_core
            else "raw window length (no motif configured)"
        )
        add(_row("Step 1 — grouped by length", f"{_plural(len(t.bins), 'bin')}, using {basis}"))
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
                "Step 2 — split by sequence",
                f"{outcome} (identity below {t.identity_threshold:.2f} separates)",
            )
        )
        if t.consensus_backend:
            add(_row("Step 3 — consensus per cluster", t.consensus_backend))
        add("")

        # 4. Per-cluster detail ---------------------------------------------
        add(_row("Candidates, strongest first", "reads shown are per-allele coverage"))
        for c in t.clusters:
            mark = "*" if c.called else " "
            hp = f"HP1 {c.n_hp1} / HP2 {c.n_hp2} / untagged {c.n_untagged}"
            add(
                f"    {mark} #{c.index}  {_plural(c.n_reads, 'read'):>10}  "
                f"{c.fraction:5.1%}  {c.length_bp:>4} bp  {hp}"
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
    add(_row("Genotype", f"{_genotype(t)}   [{t.call_rule}]"))
    if t.counts is not None:
        add(_row("Coverage", _coverage_line(t)))
    if t.tri_type:
        add(_row("Triallelic pattern", t.tri_type))
    for severity, code in t.flags:
        add(_row(f"Flag ({severity})", code, indent=6))
    return "\n".join(lines)


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
    return "\n".join(lines)
