import logging

import log

print("This is your log.py file speaking...")
logger = logging.getLogger(__name__)
logger.debug("Logging from _log_ to console started")
