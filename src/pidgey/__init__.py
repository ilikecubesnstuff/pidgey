"""Python interface for dynamics in galaxies, supporting multiple backends."""

__all__ = [
    "get_backend_from",
    # Backends
    "AgamaBackend",
    "GalaBackend",
    "GalpyBackend",
]

import inspect as _inspect
from collections.abc import Sequence as _Sequence

from .agama_backend import AgamaBackend
from .gala_backend import GalaBackend
from .galpy_backend import GalpyBackend


def get_backend_from(obj):
    """
    Get the backend associated with an object's module of origin.

    Args:
        obj (object): An object whose module is used to determine the appropriate backend.

    Returns:
        Backend: The backend associated with the object's module.
    """
    if isinstance(obj, _Sequence) and len(obj) > 0:
        obj = obj[0]
    module = _inspect.getmodule(obj) or _inspect.getmodule(type(obj))
    if module is None:
        return None

    pkg, *_ = module.__name__.partition(".")
    if pkg == "agama":
        return AgamaBackend()
    elif pkg == "gala":
        return GalaBackend()
    elif pkg == "galpy":
        return GalpyBackend()
    return None
