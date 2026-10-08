#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Terminal presentation helpers: banner, section rules and coloured logging.

Colours are used only when the target stream is a TTY and not disabled via the
``NO_COLOR`` environment variable. On Windows, ANSI/VT mode is enabled for the
console when possible. Every helper degrades to plain text gracefully, so the
output stays readable when redirected to a file.
"""

import logging
import os
import sys

_WIDTH = 72

_RESET = "\x1b[0m"
_DIM = "\x1b[2m"
_CYAN = "\x1b[36m"
_YELLOW = "\x1b[33m"
_RED = "\x1b[31m"
_BOLD = "\x1b[1m"


def _enable_windows_vt(stream) -> None:
    """Best-effort enabling of ANSI escape sequences on the Windows console."""
    try:
        import ctypes

        kernel32 = ctypes.windll.kernel32
        handle = kernel32.GetStdHandle(-11 if stream is sys.stdout else -12)
        mode = ctypes.c_uint32()
        if kernel32.GetConsoleMode(handle, ctypes.byref(mode)):
            kernel32.SetConsoleMode(handle, mode.value | 0x0004)  # ENABLE_VIRTUAL_TERMINAL_PROCESSING
    except Exception:
        pass


def ansi_enabled(stream) -> bool:
    """True if ANSI colours should be used for the given stream."""
    if os.environ.get("NO_COLOR"):
        return False
    try:
        if not stream.isatty():
            return False
    except Exception:
        return False
    if os.name == "nt":
        _enable_windows_vt(stream)
    return True


class ColorFormatter(logging.Formatter):
    """Log formatter with a dim timestamp and a coloured, padded level name."""

    _LEVEL_COLORS = {
        logging.DEBUG: _DIM,
        logging.INFO: _CYAN,
        logging.WARNING: _YELLOW,
        logging.ERROR: _RED,
        logging.CRITICAL: _BOLD + _RED,
    }

    def __init__(self, use_color: bool = True):
        super().__init__(fmt="%(asctime)s %(levelname)-8s %(message)s", datefmt="%H:%M:%S")
        self._use_color = use_color

    def format(self, record: logging.LogRecord) -> str:
        ts = self.formatTime(record, self.datefmt)
        level = f"{record.levelname:<8}"
        message = record.getMessage()
        if self._use_color:
            ts = f"{_DIM}{ts}{_RESET}"
            level = f"{self._LEVEL_COLORS.get(record.levelno, '')}{level}{_RESET}"
        line = f"{ts} {level} {message}"
        if record.exc_info:
            line += "\n" + self.formatException(record.exc_info)
        return line


def section(title: str, width: int = _WIDTH) -> str:
    """A log-friendly section rule, e.g. ``── Job mini ──────────────``."""
    dashes = max(0, width - len(title) - 4)
    return f"── {title} " + "─" * dashes


def _clip(text: str, width: int) -> str:
    text = str(text)
    return text if len(text) <= width else text[: max(0, width - 1)] + "…"


def banner(title: str, subtitle: str, rows) -> str:
    """A double-line box with a title, a subtitle and key/value rows."""
    content_w = _WIDTH - 4
    value_w = content_w - 13  # room for a 12-char key + separating space
    lines = [
        f"╔{'═' * (_WIDTH - 2)}╗",
        f"║ {_clip(title, content_w):<{content_w}} ║",
        f"║ {_clip(subtitle, content_w):<{content_w}} ║",
        f"╠{'═' * (_WIDTH - 2)}╣",
    ]
    for key, value in rows:
        lines.append(f"║ {_clip(key, 12):<12} {_clip(value, value_w):<{value_w}} ║")
    lines.append(f"╚{'═' * (_WIDTH - 2)}╝")
    return "\n".join(lines)
