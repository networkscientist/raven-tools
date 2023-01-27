import logging

from model import raven_model

print("This is package _model_ speaking...")
logger = logging.getLogger(__name__)
logger.debug("Logging from _model_ to console started")
