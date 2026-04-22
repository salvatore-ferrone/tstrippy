"""Helpers implemented in pure Python under ``tstrippy.code``."""

from importlib import import_module

__all__ = ["sampling", "bfe", "orbits"]


def __getattr__(name):
	if name in __all__:
		return import_module(f".{name}", __name__)
	raise AttributeError(f"module {__name__!r} has no attribute {name!r}")


def __dir__():
	return sorted(set(globals()) | set(__all__))