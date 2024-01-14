import astropy.coordinates as coord
import astropy.units as u
import galpy
import numpy as np
import pytest
from astropy.coordinates import CartesianDifferential

from pidgey import GalpyBackend, get_backend_from


@pytest.fixture
def backend():
    return GalpyBackend()


@pytest.fixture
def potential():
    return galpy.potential.MWPotential2014


def test_detection(potential):
    backend = get_backend_from(potential)
    assert isinstance(backend, GalpyBackend)


def test_single_orbit(potential, backend):
    STEPS = 50
    c = coord.SkyCoord(
        ra=20.0 * u.deg,
        dec=30.0 * u.deg,
        distance=2.0 * u.kpc,
        pm_ra_cosdec=-10.0 * u.mas / u.yr,
        pm_dec=20.0 * u.mas / u.yr,
        radial_velocity=50.0 * u.km / u.s,
        galcen_distance=8.0 * u.kpc,
        z_sun=15.0 * u.pc,
        galcen_v_sun=CartesianDifferential([10.0, 235.0, 7.0] * u.km / u.s),
    )
    orbit = backend.compute_orbit(c, potential, dt=0.1 * u.Gyr, steps=STEPS)
    assert isinstance(orbit, coord.representation.CartesianRepresentation)
    assert len(orbit) == STEPS


def test_multiple_orbits(potential, backend):
    ORBITS = 10
    STEPS = 50

    xs = np.ones(ORBITS) * u.kpc
    ys = np.ones(ORBITS) * u.kpc
    zs = np.ones(ORBITS) * u.kpc
    v_xs = np.zeros(ORBITS) * u.km / u.s
    v_ys = np.linspace(100, 300, ORBITS) * u.km / u.s
    v_zs = np.zeros(ORBITS) * u.km / u.s
    c = coord.SkyCoord(
        x=xs,
        y=ys,
        z=zs,
        v_x=v_xs,
        v_y=v_ys,
        v_z=v_zs,
        frame="galactocentric",
        representation_type="cartesian",
    )

    orbit = backend.compute_orbit(c, potential, dt=0.1 * u.Gyr, steps=STEPS)
    assert isinstance(orbit, coord.representation.CartesianRepresentation)
    assert len(orbit) == ORBITS
    assert len(orbit[0]) == STEPS
