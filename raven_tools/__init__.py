"""
API reference documentation for the 'raven_tools' package.
"""
__all__ = ['processing', 'raven_diag', 'model', 'config', 'log']

import logging
import sys

import config
from raven_tools import model

logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)

handler = logging.StreamHandler(sys.stdout)
frm = logging.Formatter("{levelname}: {message} ({filename}/{funcName}/{lineno})",
                        style="{")
handler.setFormatter(frm)
logger.addHandler(handler)
logger.debug("Logging from _raven_tools_ to console started")
