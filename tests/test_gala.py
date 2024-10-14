import astropy.coordinates as coord
import astropy.units as u
import numpy as np
import pytest
from astropy.coordinates import CartesianDifferential
from gala import integrate as gi
from gala import potential as gp
from gala import units as gu

from pidgey import GalaBackend, get_backend_from


@pytest.fixture
def backend():
    return GalaBackend()


@pytest.fixture
def potential():
    return gp.NFWPotential.from_circular_velocity(
        v_c=200 * u.km / u.s, r_s=10.0 * u.kpc, units=gu.galactic
    )


def test_detection(potential):
    backend = get_backend_from(potential)
    assert isinstance(backend, GalaBackend)


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
    orbit = backend.get_points()
    assert isinstance(orbit, coord.representation.CartesianRepresentation)
    assert len(orbit) == STEPS


def test_single_orbit_with_alternate_integrator(potential, backend):
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
    orbit = backend.compute_orbit(
        c,
        potential,
        dt=0.1 * u.Gyr,
        steps=STEPS,
        integrator=gi.LeapfrogIntegrator,
    )
    assert isinstance(orbit, coord.representation.CartesianRepresentation)
    assert len(orbit) == STEPS
    orbit = backend.get_points()
    assert isinstance(orbit, coord.representation.CartesianRepresentation)
    assert len(orbit) == STEPS


def test_single_orbit_with_pattern_speed(potential, backend):
    # NOTE: potential is not rotating, only tests that the pattern speed is accepted by the backend
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
    orbit = backend.compute_orbit(
        c,
        potential,
        dt=0.1 * u.Gyr,
        steps=STEPS,
        pattern_speed=50 * u.km / u.s / u.kpc,
    )
    assert isinstance(orbit, coord.representation.CartesianRepresentation)
    assert len(orbit) == STEPS
    orbit = backend.get_points()
    assert isinstance(orbit, coord.representation.CartesianRepresentation)
    assert len(orbit) == STEPS


@pytest.mark.parametrize("orbits", [1, 10])
def test_multiple_orbits(potential, backend, orbits):
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
    print(orbit.T)
    assert isinstance(orbit, coord.representation.CartesianRepresentation)
    assert len(orbit) == ORBITS
    assert len(orbit[0]) == STEPS
    orbit = backend.get_points()
    assert isinstance(orbit, coord.representation.CartesianRepresentation)
    assert len(orbit) == ORBITS
    assert len(orbit[0]) == STEPS
