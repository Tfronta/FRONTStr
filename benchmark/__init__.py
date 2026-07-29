"""Development benchmark: score FRONTStr against other callers on known samples.

**This is not part of FRONTStr and never runs for a user of the caller.**

FRONTStr is for anyone with their own ONT samples, where there is no truth
table, no Illumina profile and no second caller to compare against. There, the
caller has to stand on the evidence it shows for its own calls — which is what
``frontstr interpret --trace`` is for.

What lives here is the opposite situation: a handful of public 1000 Genomes
samples that happen to have Illumina and other-caller genotypes published, used
during development to answer "is the caller right". That is benchmarking, and
keeping it under its own name is the point — a comparison against longTR or
STRspy must never look like a step in calling a sample.

``pyproject.toml`` ships only ``frontstr``, so nothing here is installed with
the package and nothing in the caller imports it.
"""
