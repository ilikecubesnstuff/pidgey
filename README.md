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

## Table of Contents

- [Installation](#installation)
- [Usage](#usage)
    - [With `agama`](#with-agama)
    - [With `gala`](#with-gala)
    - [With `galpy`](#with-galpy)

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

### With `agama`

```py
import agama

bar_par = dict(
    type="Ferrers",
    mass=1e9,
    scaleRadius=1.0,
    axisRatioY=0.5,
    axisratioz=0.4,
    cutoffStrength=2.0,
    patternSpeed=30,
)
disk_par = dict(type="Disk", mass=5e10, scaleRadius=3, scaleHeight=0.4)
bulge_par = dict(type="Sersic", mass=1e10, scaleRadius=1, axisRatioZ=0.6)
halo_par = dict(type="NFW", mass=1e12, scaleRadius=20, axisRatioZ=0.8)
potgal = agama.Potential(disk_par, bulge_par, halo_par, bar_par)

ic = c.SkyCoord(
    x=6 * u.kpc,
    y=0 * u.kpc,
    z=2 * u.kpc,
    v_x=0 * u.km / u.s,
    v_y=150 * u.km / u.s,
    v_z=0 * u.km / u.s,
    frame="galactocentric",
    representation_type="cartesian",
)
orbit = calculate_orbit(ic, potential)
```

### With `gala`

```py
import gala.potential as gp
from gala.units import galactic

disk = gp.MiyamotoNagaiPotential(m=6e10 * u.Msun, a=3.5 * u.kpc, b=280 * u.pc, units=galactic)
halo = gp.NFWPotential(m=6e11 * u.Msun, r_s=20.0 * u.kpc, units=galactic)
bar = gp.LongMuraliBarPotential(
    m=1e10 * u.Msun,
    a=4 * u.kpc,
    b=0.8 * u.kpc,
    c=0.25 * u.kpc,
    alpha=25 * u.degree,
    units=galactic,
)
pot = gp.CCompositePotential()
pot["disk"] = disk
pot["halo"] = halo
pot["bar"] = bar

frame = gp.ConstantRotatingFrame(Omega=[0, 0, 30] * u.km / u.s / u.kpc, units=galactic)
hamiltonian = gp.Hamiltonian(pot, frame=frame)
# this hamiltonian is pidgey's "potential"
# since it has an orbit integration method

ic = c.SkyCoord(
    x=6 * u.kpc,
    y=0 * u.kpc,
    z=2 * u.kpc,
    v_x=0 * u.km / u.s,
    v_y=150 * u.km / u.s,
    v_z=0 * u.km / u.s,
    frame="galactocentric",
    representation_type="cartesian",
)
orbit = calculate_orbit(ic, hamiltonian)
```


### With `galpy`

```py
import galpy.potential as gp

omega = 30 * u.km / u.s / u.kpc
halo = gp.NFWPotential(conc=10, mvir=1)
disc = gp.MiyamotoNagaiPotential(amp=5e10 * u.solMass, a=3 * u.kpc, b=0.1 * u.kpc)
bar = gp.SoftenedNeedleBarPotential(
    amp=1e9 * u.solMass, a=1.5 * u.kpc, b=0 * u.kpc, c=0.5 * u.kpc, omegab=omega
)
potential = [halo, disc, bar]

ic = c.SkyCoord(
    x=6 * u.kpc,
    y=0 * u.kpc,
    z=2 * u.kpc,
    v_x=0 * u.km / u.s,
    v_y=150 * u.km / u.s,
    v_z=0 * u.km / u.s,
    frame="galactocentric",
    representation_type="cartesian",
)
orbit = calculate_orbit(ic, potential)
```
