"""Every knob a run can turn, its default, and where that default came from.

One table, because the parameters were scattered across five modules and only
some of them were reachable from the CLI. A laboratory that cannot see what a
run used cannot reproduce it, and a laboratory that cannot change a threshold
cannot test one.

Provenance is the point
-----------------------

Not all defaults are equal, and treating them as if they were is how an
exploratory profile ends up filed as a production one. Three kinds:

``derived``
    Computed from measured data. ``low_coverage_reads`` = 20 comes from a
    binomial dropout calculation; ``min_reads_third`` = 5 comes from the
    known-bug #6 investigation into ONT basecaller phantoms. Lowering one of
    these is a real change to what the caller will believe, so it **marks the
    run** — see :meth:`RunParameters.derived_overrides`.

``chosen``
    Defensible, but picked rather than measured. The analytical and calling
    thresholds are the honest examples; the documentation has said so from the
    start. Changing them is ordinary tuning.

``convention``
    Inherited from forensic practice rather than from this data.
    ``min_phr_for_het`` = 0.4 is the capillary-electrophoresis heterozygote
    balance rule, which is exactly why phasing is allowed to override it.

The point is not to stop anyone changing a number. It is that six months later
the report says which numbers were changed, and that the ones with measured
backing are called out rather than buried in a list of forty.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Literal

Provenance = Literal["derived", "chosen", "convention"]


@dataclass(frozen=True, slots=True)
class ParamSpec:
    """What a parameter means, defaults to, and why."""

    name: str
    default: Any
    provenance: Provenance
    note: str


#: The canonical table. Order is the order a reader wants: what reads are kept,
#: how they are grouped, what counts as an allele, then the QC policy.
PARAM_SPECS: tuple[ParamSpec, ...] = (
    ParamSpec(
        "min_mapq",
        20,
        "chosen",
        "Reads below this MAPQ are dropped. Also what excludes X/Y paralogue "
        "mismappings at the sex markers.",
    ),
    ParamSpec(
        "flank_anchor",
        20,
        "chosen",
        "bp of cleanly aligned flank required on each side of the window.",
    ),
    ParamSpec(
        "identity_threshold",
        0.97,
        "chosen",
        "Pairwise identity to join a cluster seed. Known to be miscalibrated: "
        "it was set as if comparing a read to a consensus, but the comparison "
        "is raw read to raw read, where two ONT reads of one allele diverge "
        "2-4%. Worth experimenting with.",
    ),
    ParamSpec(
        "len_tolerance_bp",
        0,
        "derived",
        "Merges adjacent length bins. Must stay 0: with repeat-core binning the "
        "flank errors it was meant to absorb are already gone, and any tolerance "
        "merges genuine microvariants — TH01 9 and 9.3 are 3 bp apart.",
    ),
    ParamSpec(
        "analytical_thresh",
        0.02,
        "chosen",
        "Below this fraction of locus coverage a cluster is noise.",
    ),
    ParamSpec(
        "calling_thresh",
        0.10,
        "chosen",
        "Below this fraction (but above analytical) a cluster is an artefact. "
        "Not data-derived; a chosen default.",
    ),
    ParamSpec(
        "min_phr_for_het",
        0.4,
        "convention",
        "Minor allele as a fraction of the major before a heterozygote is "
        "called on read counts. Inherited from capillary electrophoresis, which "
        "is why phasing evidence is allowed to override it.",
    ),
    ParamSpec(
        "min_reads_third",
        5,
        "derived",
        "Absolute read floor a 3rd candidate must clear before it can promote a "
        "locus to triallelic or raise a mixture flag. From the known-bug #6 "
        "work: 2-4 read basecaller phantoms sit at 5-10% allele fraction at ONT "
        "coverage. Lowering it re-admits those phantoms.",
    ),
    ParamSpec(
        "low_coverage_reads",
        20,
        "derived",
        "Read floor below which a called locus is flagged, measured against the "
        "reads supporting the genotype rather than every read spanning the "
        "window. Re-derived by binomial against ONT depths: at the most "
        "unbalanced heterozygote the caller accepts, the minor allele is missed "
        "~5.7% of the time from N=20 upward, ~9.7% from 17 and ~4.6% from 25. "
        "20 is the knee — below it the risk climbs fast, above it barely falls.",
    ),
    ParamSpec(
        "balanced_ab_max",
        0.65,
        "chosen",
        "Largest allele balance a heterozygote may have and still count as "
        "balanced. 0.50 is even; 0.714 is where min_phr_for_het stops calling "
        "hets at all, so this band warns inside the callable range.",
    ),
)

_BY_NAME = {spec.name: spec for spec in PARAM_SPECS}


@dataclass(slots=True)
class RunParameters:
    """The values a run actually used, against the table above."""

    values: dict[str, Any] = field(default_factory=dict)

    @classmethod
    def defaults(cls) -> RunParameters:
        return cls({spec.name: spec.default for spec in PARAM_SPECS})

    @classmethod
    def of(cls, **kwargs: Any) -> RunParameters:
        """Build from explicit values, filling anything unset with its default."""
        unknown = set(kwargs) - set(_BY_NAME)
        if unknown:
            raise KeyError(f"unknown parameter(s): {sorted(unknown)}")
        merged = {spec.name: spec.default for spec in PARAM_SPECS}
        merged.update({k: v for k, v in kwargs.items() if v is not None})
        return cls(merged)

    def __getitem__(self, name: str) -> Any:
        return self.values[name]

    def spec(self, name: str) -> ParamSpec:
        return _BY_NAME[name]

    def is_default(self, name: str) -> bool:
        return bool(self.values[name] == _BY_NAME[name].default)

    def overrides(self) -> list[ParamSpec]:
        """Specs whose value differs from the default, in table order."""
        return [s for s in PARAM_SPECS if self.values.get(s.name) != s.default]

    def derived_overrides(self) -> list[ParamSpec]:
        """The overrides that change a *measured* default. These mark the run."""
        return [s for s in self.overrides() if s.provenance == "derived"]

    def as_dict(self) -> dict[str, Any]:
        """Flat mapping for the payload, audit record and report."""
        return dict(self.values)

    def as_audit_rows(self) -> list[dict[str, Any]]:
        """One row per parameter: value, default, provenance, and why."""
        return [
            {
                "name": s.name,
                "value": self.values.get(s.name),
                "default": s.default,
                "is_default": self.is_default(s.name),
                "provenance": s.provenance,
                "note": s.note,
            }
            for s in PARAM_SPECS
        ]


def render_echo(params: RunParameters, *, width: int = 26) -> str:
    """The block printed at the start of a run.

    Shows every parameter, not only the overridden ones: the values that decide
    a call are usually the ones nobody typed, so a run that lists only its
    command line reads as if it were barely configured.
    """
    lines = ["Parameters in force"]
    for spec in PARAM_SPECS:
        value = params[spec.name]
        if params.is_default(spec.name):
            lines.append(f"  {spec.name.ljust(width)}{value}")
        else:
            lines.append(
                f"  {spec.name.ljust(width)}{value}"
                f"   CHANGED (default {spec.default}, {spec.provenance})"
            )
    derived = params.derived_overrides()
    if derived:
        lines.append("")
        lines.append(
            f"  {len(derived)} measured default(s) overridden — this run is not "
            "comparable with a default one:"
        )
        for spec in derived:
            lines.append(f"    {spec.name}: {spec.note}")
    return "\n".join(lines)
