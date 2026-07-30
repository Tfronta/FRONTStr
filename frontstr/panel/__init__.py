"""Panel and marker data model.

A *panel* is an ordered collection of *markers* (TR loci) with forensic
metadata: motif(s), corr_value, CODIS name, allow_triallelic, etc.

"""

from frontstr.panel.models import Panel, System

__all__ = ["Panel", "System"]
