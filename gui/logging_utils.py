from __future__ import annotations

import logging
from pathlib import Path


_LOGGERS = {}


def get_gui_logger(log_path: Path) -> logging.Logger:
    key = str(log_path)
    if key in _LOGGERS:
        return _LOGGERS[key]

    logger = logging.getLogger(f"merlino.gui.{log_path.name}")
    logger.setLevel(logging.INFO)
    logger.propagate = False

    if not any(isinstance(h, logging.FileHandler) and h.baseFilename == str(log_path) for h in logger.handlers):
        log_path.parent.mkdir(parents=True, exist_ok=True)
        handler = logging.FileHandler(log_path, encoding="utf-8")
        handler.setLevel(logging.INFO)
        formatter = logging.Formatter("%(asctime)s %(message)s", datefmt="%Y-%m-%d %H:%M:%S")
        handler.setFormatter(formatter)
        logger.addHandler(handler)

    _LOGGERS[key] = logger
    return logger
