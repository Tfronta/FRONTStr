"""Structured process logging.

Layer 1 of the audit trail (see :mod:`frontstr.audit`): a JSONL record of what
the pipeline did, one event per line.

Deliberately free of any FRONTStr import. The interpretation layer logs, and
the audit record describes the interpretation layer's output — putting both in
one module makes the domain depend on the audit trail and the audit trail
depend on the domain.

Two sinks, two renderings
-------------------------

The file sink and the terminal sink want opposite things. The audit file wants
JSONL with sorted keys, because that is what diffs cleanly between runs and
what a downstream tool can parse. A person watching a run wants aligned
key=value pairs they can read at a glance.

Rather than pick one and make the other unpleasant, the shared processor chain
stops short of rendering and hands off to :class:`structlog.stdlib.
ProcessorFormatter`, so each handler renders the same event its own way.
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


#: Everything both sinks share. Rendering is deliberately *not* here — it is
#: per-handler, see the module docstring.
_SHARED_PROCESSORS: list[Processor] = [
    structlog.contextvars.merge_contextvars,
    structlog.processors.add_log_level,
    structlog.processors.TimeStamper(fmt="iso", utc=True),
    structlog.processors.StackInfoRenderer(),
    structlog.processors.format_exc_info,
]

#: ``sort_keys`` matters beyond tidiness: a JSONL log with a stable key order
#: diffs cleanly between runs, which is how you find what changed.
_JSON_RENDERER: Processor = structlog.processors.JSONRenderer(sort_keys=True)


def _apply(level: int) -> None:
    structlog.configure(
        processors=[*_SHARED_PROCESSORS, structlog.stdlib.ProcessorFormatter.wrap_for_formatter],
        wrapper_class=structlog.make_filtering_bound_logger(level),
        logger_factory=structlog.stdlib.LoggerFactory(),
        cache_logger_on_first_use=False,
    )


def _formatter(renderer: Processor) -> logging.Formatter:
    return structlog.stdlib.ProcessorFormatter(processor=renderer, foreign_pre_chain=[])


def configure_logging(
    log_path: Path | None = None, *, level: int = logging.INFO, console: bool = False
) -> None:
    """Set up structured logging for a run.

    Args:
        log_path: JSONL destination. ``None`` leaves the log without a file
            sink, which is what library consumers and tests want.
        level: Standard library level; ``DEBUG`` adds per-marker events.
        console: Also emit lines on stderr, rendered for reading rather than
            for parsing. Off by default so the CLI's own Rich output stays the
            primary channel and log lines do not interleave with result tables.
            **stderr, not stdout** — piping the table somewhere must not pick
            up the log.

    Safe to call more than once; the last call wins.
    """
    handlers: list[logging.Handler] = []
    if log_path is not None:
        log_path.parent.mkdir(parents=True, exist_ok=True)
        file_handler = logging.FileHandler(log_path, mode="w", encoding="utf-8")
        file_handler.setFormatter(_formatter(_JSON_RENDERER))
        handlers.append(file_handler)
    if console or not handlers:
        stream = logging.StreamHandler(sys.stderr)
        # Colour only when stderr is a terminal: a redirected log full of ANSI
        # escapes is worse than no colour at all.
        stream.setFormatter(
            _formatter(structlog.dev.ConsoleRenderer(colors=sys.stderr.isatty()))
            if console
            else _formatter(_JSON_RENDERER)
        )
        handlers.append(stream)

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
