# pidgey

[![PyPI - Version](https://img.shields.io/pypi/v/pidgey)](https://pypi.org/project/pidgey/)
[![PyPI - Python Version](https://img.shields.io/pypi/pyversions/pidgey)](https://pypi.org/project/pidgey/)
[![tests](https://github.com/ilikecubesnstuff/pidgey/actions/workflows/tests.yml/badge.svg)](https://github.com/ilikecubesnstuff/pidgey/actions/workflows/tests.yml)
[![pdm-managed](https://img.shields.io/badge/pdm-managed-blueviolet)](https://pdm.fming.dev)
[![pre-commit](https://img.shields.io/badge/pre--commit-enabled-brightgreen?logo=pre-commit)](https://github.com/pre-commit/pre-commit)
[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)
[![Imports: isort](https://img.shields.io/badge/%20imports-isort-%231674b1?style=flat&labelColor=ef8336)](https://pycqa.github.io/isort/)

A python interface for dynamics of galaxies, using existing galactic dynamics libraries in Python ([`galpy`](https://github.com/jobovy/galpy), [`gala`](https://github.com/adrn/gala) and [`agama`](https://github.com/GalacticDynamics-Oxford/Agama)) as backends.

This is currently used in the [commensurability](https://github.com/ilikecubesnstuff/commensurability/) library.

## Installation

Install this package via `pip`:

```
python -m pip install pidgey
```

## Usage

Get the points from an orbit integration stored in an `astropy.coordinates.SkyCoord` object by calling the `compute_orbit` method of a backend. The backend can generally be obtained from the potential used for the orbit integration with the `get_backend_from` function.

```py
import astropy.coordinates as c
import astropy.units as u

from pidgey import get_backend_from

def calculate_orbit(ic: c.SkyCoord,
                    potential: Any
                    ) -> c.SkyCoord:
    backend = get_backend_from(potential)
    orbit = backend.compute_orbit(ic, potential, dt=0.01 * u.Gyr, steps=1000)
    return orbit
```
