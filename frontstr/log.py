"""Structured process logging.

Layer 1 of the audit trail (see :mod:`frontstr.audit`): a JSONL record of what
the pipeline did, one event per line.

Deliberately free of any FRONTStr import. The interpretation layer logs, and
the audit record describes the interpretation layer's output — putting both in
one module makes the domain depend on the audit trail and the audit trail
depend on the domain.
"""

from __future__ import annotations

import logging
import sys
from pathlib import Path
from typing import Any

import structlog
from structlog.typing import Processor

#: Default filename for the per-run process log, written into the output dir.
PROCESS_LOG_NAME = "frontstr.log.jsonl"


#: ``sort_keys`` matters beyond tidiness: a JSONL log with a stable key order
#: diffs cleanly between runs, which is how you find what changed.
_JSONL_PROCESSORS: list[Processor] = [
    structlog.contextvars.merge_contextvars,
    structlog.processors.add_log_level,
    structlog.processors.TimeStamper(fmt="iso", utc=True),
    structlog.processors.StackInfoRenderer(),
    structlog.processors.format_exc_info,
    structlog.processors.JSONRenderer(sort_keys=True),
]


def _apply(level: int) -> None:
    structlog.configure(
        processors=_JSONL_PROCESSORS,
        wrapper_class=structlog.make_filtering_bound_logger(level),
        logger_factory=structlog.stdlib.LoggerFactory(),
        cache_logger_on_first_use=False,
    )


def configure_logging(
    log_path: Path | None = None, *, level: int = logging.INFO, console: bool = False
) -> None:
    """Set up structured logging for a run.

    Args:
        log_path: JSONL destination. ``None`` leaves the log without a file
            sink, which is what library consumers and tests want.
        level: Standard library level; ``DEBUG`` adds per-marker events.
        console: Also emit lines on stderr. Off by default so the CLI's own
            Rich output stays the primary channel and log lines do not
            interleave with progress tables.

    Safe to call more than once; the last call wins.
    """
    handlers: list[logging.Handler] = []
    if log_path is not None:
        log_path.parent.mkdir(parents=True, exist_ok=True)
        handlers.append(logging.FileHandler(log_path, mode="w", encoding="utf-8"))
    if console or not handlers:
        handlers.append(logging.StreamHandler(sys.stderr))

    logging.basicConfig(level=level, format="%(message)s", handlers=handlers, force=True)
    _apply(level)


#: Importing FRONTStr must not start printing. A library that configures
#: logging on the application's behalf is a library that writes to someone
#: else's stdout — so the default drops everything below WARNING, and
#: :func:`configure_logging` is what opens the tap.
_apply(logging.WARNING)


def get_logger(name: str = "frontstr") -> Any:
    """Return a bound structlog logger."""
    return structlog.get_logger(name)
