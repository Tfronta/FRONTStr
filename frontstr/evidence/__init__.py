"""Sequence-pileup evidence layer (Layer 2).

This is the unique forensic value of FRONTStr: instead of approximating
per-allele coverage from VCF fields, we re-read the BAM, extract the TR
subsequence per read, cluster by sequence (with ONT-aware tolerance), and
emit integer counts.

Phase 1.4 of ROADMAP.md.
"""
