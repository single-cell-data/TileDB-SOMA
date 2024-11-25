"""Experimental SOMA features

This module is for experimental features. Support for these features
may be dropped.
"""

from .ingest import VisiumPaths, from_visium

__all__ = ["from_visium", "VisiumPaths"]
