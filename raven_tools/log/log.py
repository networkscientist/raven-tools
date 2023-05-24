import logging
import sys

logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)

handler = logging.StreamHandler(sys.stdout)
frm = logging.Formatter("{levelname}: {message} ({filename}/{funcName}/{lineno})",
                        style="{")
handler.setFormatter(frm)
logger.addHandler(handler)
