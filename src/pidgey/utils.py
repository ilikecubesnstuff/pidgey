import inspect
from collections.abc import Collection, Iterator, Sequence
from itertools import islice
from math import ceil

import astropy.units as u
import numpy as np


def make_quantity(obj, unit: u.Unit = u.dimensionless_unscaled):
    if isinstance(obj, u.Quantity):
        return obj
    return obj * unit


def make_collection(obj, cls: type = list):
    if not issubclass(cls, Collection):
        raise ValueError("Class provided must subclass collections.abc.Collection")

    # mypy doesn't let me instantiate cls with arguments?
    # look into this later
    if not isinstance(obj, Collection):
        obj = [obj]
    elif isinstance(obj, np.ndarray) and obj.size == 1 and obj.ndim == 0:
        obj = [obj]
    elif len(obj) == 0:
        obj = [obj]

    if isinstance(obj, u.Quantity):
        return make_collection(obj.value)
    return obj


def make_sequence(obj, cls: type = list):
    if not issubclass(obj, Sequence):
        raise ValueError("Class provided must subclass collections.abc.Sequence")
    if not isinstance(obj, Sequence):
        obj = cls([obj])
    return obj
