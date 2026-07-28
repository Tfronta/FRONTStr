"""Native, sequence-resolved VCF export.

Not an annotation of someone else's VCF. FRONTStr's calls come from its own
evidence layer, so this file exists whether or not LongTR was ever run — which
is the point: it is the interchange format for benchmarking FRONTStr against
other callers, and a format that only exists when a different caller was
involved cannot serve that purpose.

What makes it *sequence-resolved*
---------------------------------

``ALT`` carries the called allele's actual sequence, not a length or a repeat
count. That is the whole long-read advantage: two alleles with the same repeat
count and different internal structure — iso-alleles — are different ALTs here,
where a length-based caller would emit one. The repeat count is still reported,
in ``FORMAT/MC``, but as an annotation of the sequence rather than as the call.

The variant interval is the panel window
----------------------------------------

``POS``/``INFO/END`` span the panel's extraction window, and ``REF`` is the
reference over exactly that interval. This is deliberate: a cluster consensus
*is* the window's sequence as observed in the reads, so REF and every ALT
describe the same interval and are directly comparable. Emitting only the
repeat core would be more compact, but FRONTStr never computes the core's
reference coordinates — only its offset within a read — so the core's POS would
be a guess.

The cost is bulk: REF and ALT run to the window width, typically ~250 bp.

Why a reference FASTA is mandatory
----------------------------------

A VCF whose REF is not the reference sequence is not a VCF. There is no way to
derive REF from a BAM pileup, so rather than emit a plausible-looking file with
a placeholder REF — which would validate, load, and be silently wrong in any
comparison — this function requires the FASTA and fails clearly without it.

Excluded records
----------------

Amelogenin is skipped. It is not a tandem repeat; its "alleles" are the letters
X and Y, which are not nucleotide sequences, and forcing them into REF/ALT
would produce a record that parses but means nothing. Sex typing travels in the
JSON, CSV and XLSX exports instead. A header line records the omission so it is
visible rather than merely absent.
"""

from __future__ import annotations

from datetime import UTC, datetime
from pathlib import Path
from typing import Any

from frontstr.errors import FrontstrError
from frontstr.version import __version__

#: Characters that must not appear raw in a VCF INFO/FORMAT value.
#: Percent-encoded per the VCF 4.3 convention so the original string is
#: recoverable — ISFG notation is space-separated, and losing those spaces
#: would make the field unparseable back into brackets.
_PERCENT_ENCODE = {
    "%": "%25",
    ":": "%3A",
    ";": "%3B",
    "=": "%3D",
    ",": "%2C",
    " ": "%20",
    "\t": "%09",
    "\n": "%0A",
    "\r": "%0D",
}

_MISSING = "."


def _encode(value: str) -> str:
    """Percent-encode a string for use as a VCF INFO/FORMAT value."""
    if not value:
        return _MISSING
    return "".join(_PERCENT_ENCODE.get(c, c) for c in value)


def decode_vcf_value(value: str) -> str:
    """Inverse of :func:`_encode`. Provided so consumers can round-trip."""
    out = value
    for raw, encoded in _PERCENT_ENCODE.items():
        out = out.replace(encoded, raw)
    return out


def _header(
    payload: dict[str, Any],
    reference_fasta: Path,
    contigs: list[tuple[str, int]],
    sample: str,
) -> list[str]:
    meta = payload.get("meta", {})
    audit = payload.get("audit", {})
    lines = [
        "##fileformat=VCFv4.2",
        f"##fileDate={datetime.now(UTC).strftime('%Y%m%d')}",
        f"##source=FRONTStr {__version__}",
        f"##reference=file://{reference_fasta}",
    ]
    lines += [f"##contig=<ID={name},length={length}>" for name, length in contigs]

    lines += [
        f'##frontstrPanel="{meta.get("panel_name", "")} {meta.get("panel_version", "")}"',
        f'##frontstrReferenceBuild="{meta.get("reference_build", "")}"',
        f'##frontstrPoaBackend="{audit.get("poa_backend", "")}"',
        f'##frontstrStutterModel="{audit.get("stutter_model_version", "")} '
        f'({audit.get("stutter_model_protocol", "")})"',
        '##frontstrNote="POS..INFO/END span the panel extraction window; REF and '
        "ALT describe that whole interval, so ALT length is not the repeat length. "
        'Use FORMAT/AL for allele length and FORMAT/MC for the allele number."',
        '##frontstrNote="Amelogenin is not emitted: it is not a tandem repeat and '
        "its X/Y calls are not nucleotide sequences. See the JSON or XLSX export "
        'for sex typing."',
        '##frontstrNote="String FORMAT values are percent-encoded per VCF 4.3 '
        '(%20 space, %3A colon, %3B semicolon, %3D equals, %2C comma, %25 percent)."',
    ]

    lines += [
        '##INFO=<ID=END,Number=1,Type=Integer,Description="End of the panel window">',
        '##INFO=<ID=MARKER,Number=1,Type=String,Description="Forensic marker name">',
        '##INFO=<ID=MOTIF,Number=1,Type=String,Description="Repeat motif(s), '
        'comma-separated on the canonical strand">',
        '##INFO=<ID=PERIOD,Number=1,Type=Integer,Description="Repeat period in bp; '
        '-1 for compound markers whose allele number comes from bracket counting">',
        '##INFO=<ID=STRAND,Number=1,Type=String,Description="Strand of the canonical '
        'ISFG motif relative to the reference">',
        '##INFO=<ID=NOTE,Number=.,Type=String,Description="Informational FRONTStr '
        'flags (severity info) that do not filter the call">',
    ]

    lines += [
        '##FILTER=<ID=PASS,Description="No warning or error flag raised">',
    ]
    for code, description in _FILTER_DESCRIPTIONS.items():
        lines.append(f'##FILTER=<ID={code},Description="{description}">')

    lines += [
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype. Triallelic loci '
        'are emitted as three alleles (e.g. 1/2/3)">',
        '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Reads spanning the whole '
        'window at this locus">',
        '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Integer read depth per '
        "allele, REF first. Counted from FRONTStr's own pileup, not derived from a "
        "caller's length differences\">",
        '##FORMAT=<ID=AL,Number=R,Type=Integer,Description="Allele length in bp, one per '
        'allele in REF-then-ALT order (the window, matching REF/ALT)">',
        '##FORMAT=<ID=MC,Number=R,Type=String,Description="Canonical allele number per '
        "allele in REF-then-ALT order (ISFG designation, CE, or bracket count); "
        '. if not called">',
        '##FORMAT=<ID=MM,Number=R,Type=String,Description="How each MC value was '
        "derived: strnaming | period_ce | reference_offset | bracket_count | delta | "
        'bp_sizing">',
        '##FORMAT=<ID=ISFG,Number=R,Type=String,Description="ISFG bracketed '
        'nomenclature per allele in REF-then-ALT order, percent-encoded">',
        '##FORMAT=<ID=ISO,Number=R,Type=String,Description="Iso-allele suffix per '
        'allele in REF-then-ALT order from the curated catalog; . when unmatched">',
        '##FORMAT=<ID=HP,Number=R,Type=String,Description="Read partition per '
        'allele in REF-then-ALT order as hp1|hp2|untagged">',
        '##FORMAT=<ID=SB,Number=R,Type=String,Description="Strand partition per '
        'allele in REF-then-ALT order as forward|reverse">',
        '##FORMAT=<ID=CM,Number=R,Type=String,Description="Consensus method per '
        'allele in REF-then-ALT order: poa_spoa | poa_abpoa | single | mode">',
    ]

    lines.append(
        "#"
        + "\t".join(
            [
                "CHROM",
                "POS",
                "ID",
                "REF",
                "ALT",
                "QUAL",
                "FILTER",
                "INFO",
                "FORMAT",
                sample,
            ]
        )
    )
    return lines


#: FILTER entries. Only warn/error conditions filter a call; informational
#: flags (iso-allele, triallelic, phantom collapse) go to INFO/NOTE, because a
#: locus is not "bad" for having resolved an iso-allele — that is the tool
#: working.
_FILTER_DESCRIPTIONS = {
    "low_coverage": "Called below the coverage floor; a minor allele may have been missed",
    "dropout": "No allele called at this locus",
    "strand_bias": "A called allele's reads are strand-skewed beyond chance",
    "mixture_suspected": "More than two alleles at a locus without allow_triallelic",
    "longtr_discordant": "The evidence layer and LongTR disagree on the called alleles",
    "consensus_fallback": "No POA backend: consensus is a single unpolished read",
    "ce_nomenclature_offset": (
        "The reported allele number is not the legacy CE-kit designation for this marker"
    ),
}


def _read_reference(reference_fasta: Path, chrom: str, start_1based: int, end: int) -> str:
    """Fetch ``[start_1based, end]`` from the FASTA, in upper case."""
    import pysam

    with pysam.FastaFile(str(reference_fasta)) as fa:
        return fa.fetch(chrom, start_1based - 1, end).upper()


def _genotype(n_called: int, ref_indexes: list[int]) -> str:
    """Build the GT field from allele indexes (0 = REF)."""
    if n_called == 0:
        return "./."
    if n_called == 1:
        return f"{ref_indexes[0]}/{ref_indexes[0]}"
    return "/".join(str(i) for i in ref_indexes)


def _partition(allele: dict[str, Any]) -> str:
    return f"{allele['n_reads_hp1']}|{allele['n_reads_hp2']}|{allele['n_reads_hp_none']}"


def _record(marker: dict[str, Any], reference_fasta: Path) -> tuple[str, int, str] | None:
    """Build one VCF data line. Returns ``(chrom, pos, line)`` for sorting."""
    system = marker["system"]
    chrom = system["chromosome"]
    pos = system["ref_start"]
    end = system["ref_end"]

    try:
        ref_seq = _read_reference(reference_fasta, chrom, pos, end)
    except (KeyError, ValueError) as exc:
        raise FrontstrError(
            f"{marker['marker_name']}: cannot read {chrom}:{pos}-{end} from "
            f"{reference_fasta} ({exc}). Is this the build the panel targets?"
        ) from exc
    if not ref_seq:
        raise FrontstrError(
            f"{marker['marker_name']}: reference {chrom}:{pos}-{end} is empty in {reference_fasta}"
        )

    called = marker.get("alleles_called", [])

    # Map each called allele to an ALT index; a call identical to the reference
    # is allele 0 and must not be repeated in ALT.
    alts: list[str] = []
    indexes: list[int] = []
    for a in called:
        seq = "<DEL>" if a["is_deletion"] else a["consensus"].upper()
        if seq == ref_seq:
            indexes.append(0)
            continue
        if seq in alts:
            indexes.append(alts.index(seq) + 1)
            continue
        alts.append(seq)
        indexes.append(len(alts))

    flags = marker.get("flags", [])
    filters = sorted({f["code"] for f in flags if f["severity"] in ("warn", "error")})
    notes = sorted({f["code"] for f in flags if f["severity"] == "info"})

    info = [
        f"END={end}",
        f"MARKER={_encode(marker['marker_name'])}",
        f"MOTIF={_encode(system['motif'])}",
        f"PERIOD={system['period']}",
        f"STRAND={system.get('strand', '+')}",
    ]
    if notes:
        info.append("NOTE=" + ",".join(notes))

    # Every per-allele field is Number=R: one value per allele in *allele-index*
    # order (REF, then each ALT), with "." where that allele was not called.
    #
    # Emitting them in call order instead — which is by read depth — silently
    # misaligns them against AD whenever one of the calls happens to equal the
    # reference: MC[0] would then describe a different allele than AD[0]. A
    # benchmark that joins those two columns would read every such locus wrong,
    # and nothing about the file would look malformed.
    by_index: list[dict[str, Any] | None] = [None] * (len(alts) + 1)
    depths = [0] * (len(alts) + 1)
    for a, i in zip(called, indexes, strict=True):
        by_index[i] = a
        depths[i] += a["n_reads_total"]

    def per_allele(fn: Any) -> str:
        return ",".join(_MISSING if a is None else (fn(a) or _MISSING) for a in by_index)

    fmt_keys = ["GT", "DP", "AD", "AL", "MC", "MM", "ISFG", "ISO", "HP", "SB", "CM"]
    fmt_values = [
        _genotype(len(called), indexes),
        str(marker["total_reads"]),
        ",".join(str(d) for d in depths),
        per_allele(lambda a: str(a["length_bp"])),
        per_allele(lambda a: a.get("number_label")),
        per_allele(lambda a: a.get("number_method")),
        per_allele(lambda a: _encode(a["isfg"])),
        per_allele(lambda a: (a.get("iso") or {}).get("suffix")),
        per_allele(_partition),
        per_allele(lambda a: f"{a['n_forward']}|{a['n_reverse']}"),
        per_allele(lambda a: a.get("consensus_method")),
    ]

    line = "\t".join(
        [
            chrom,
            str(pos),
            marker["marker_name"],
            ref_seq,
            ",".join(alts) if alts else _MISSING,
            _MISSING,
            ";".join(filters) if filters else "PASS",
            ";".join(info),
            ":".join(fmt_keys),
            ":".join(fmt_values),
        ]
    )
    return chrom, pos, line


def write_run_vcf(
    payload: dict[str, Any],
    out_path: Path,
    *,
    reference_fasta: Path | None,
) -> Path:
    """Write a native, sequence-resolved VCF for one sample.

    Args:
        payload: Output of :func:`frontstr.report.payload.serialize_run`.
        out_path: Destination ``.vcf``.
        reference_fasta: Indexed reference FASTA. Required — see the module
            docstring for why there is no placeholder fallback.

    Returns:
        ``out_path``.

    Raises:
        FrontstrError: If the reference is missing, unindexed, or does not
            contain a marker's interval.
    """
    if reference_fasta is None:
        raise FrontstrError(
            "VCF export needs a reference FASTA (--reference). REF must be the "
            "actual reference sequence; there is no meaningful placeholder, and "
            "emitting one would produce a file that parses but compares wrongly."
        )
    if not reference_fasta.exists():
        raise FrontstrError(f"Reference FASTA not found: {reference_fasta}")

    try:
        import pysam
    except ImportError as exc:  # pragma: no cover - pysam is a hard dependency
        raise FrontstrError("pysam is required for VCF export") from exc

    try:
        with pysam.FastaFile(str(reference_fasta)) as fa:
            contig_lengths = dict(zip(fa.references, fa.lengths, strict=True))
    except (OSError, ValueError) as exc:
        raise FrontstrError(
            f"Cannot open {reference_fasta} — is it indexed (.fai)? ({exc})"
        ) from exc

    results = payload.get("results", [])
    records: list[tuple[str, int, str]] = []
    used_contigs: list[str] = []
    for marker in results:
        system = marker["system"]
        # Amelogenin has no tandem repeat and no nucleotide alleles.
        if system.get("marker_type", "str") != "str":
            continue
        if not marker.get("alleles_called") and marker["total_reads"] == 0:
            # No reads at all: emitting a REF/REF record would assert a
            # homozygous reference call that was never observed.
            continue
        chrom = system["chromosome"]
        if chrom not in contig_lengths:
            raise FrontstrError(
                f"{marker['marker_name']}: contig {chrom} is absent from "
                f"{reference_fasta}. Panel and reference disagree on naming "
                "(chr prefix?)."
            )
        if chrom not in used_contigs:
            used_contigs.append(chrom)
        record = _record(marker, reference_fasta)
        if record is not None:
            records.append(record)

    # VCF must be coordinate-sorted to be indexable; panel order is not.
    contig_order = {name: i for i, name in enumerate(contig_lengths)}
    records.sort(key=lambda r: (contig_order.get(r[0], 1 << 30), r[1]))

    sample = payload.get("meta", {}).get("sample_name") or "SAMPLE"
    lines = _header(
        payload,
        reference_fasta,
        [(c, contig_lengths[c]) for c in sorted(used_contigs, key=contig_order.get)],  # type: ignore[arg-type]
        sample,
    )
    lines += [r[2] for r in records]

    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return out_path
