import inspect as _inspect
from collections.abc import Sequence as _Sequence

from .agama_backend import AgamaBackend
from .gala_backend import GalaBackend
from .galpy_backend import GalpyBackend


def get_backend_from(obj):
    if isinstance(obj, _Sequence) and len(obj) > 0:
        obj = obj[0]
    module = _inspect.getmodule(obj)
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
