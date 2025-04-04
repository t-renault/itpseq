"""itpseq API entry point"""

import importlib
import sys

__version__ = '0.0.1a18'

_LAZY_MODULES = {
    '.core': ['DataSet', 'Sample', 'Replicate'],
    '.parsing': ['parse', 'parse_all'],
}

_ORIGINS = {}
for module, items in _LAZY_MODULES.items():
    for item in items:
        _ORIGINS[item] = module

__all__ = list(_ORIGINS)


def __dir__():
    return list(sys.modules.get(__name__, {}).__dict__) + __all__


def __getattr__(name):
    if name in _ORIGINS:
        module = importlib.import_module(_ORIGINS[name], package=__package__)
        return getattr(module, name) if name in dir(module) else module
    raise AttributeError(f"module '{__name__}' has no attribute '{name}'")
