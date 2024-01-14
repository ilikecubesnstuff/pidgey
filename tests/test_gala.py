import astropy.coordinates as coord
import astropy.units as u
import gala
import numpy as np
import pytest
from astropy.coordinates import CartesianDifferential

from pidgey import GalaBackend, get_backend_from


@pytest.fixture
def backend():
    return GalaBackend()


@pytest.fixture
def potential():
    return gala.potential.NFWPotential.from_circular_velocity(
        v_c=200 * u.km / u.s, r_s=10.0 * u.kpc, units=gala.units.galactic
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
    orbit = backend.compute_orbit(c, 0.1 * u.Gyr, STEPS, potential)
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

    orbit = backend.compute_orbit(c, 0.1 * u.Gyr, STEPS, potential)
    print(orbit.T)
    assert isinstance(orbit, coord.representation.CartesianRepresentation)
    assert len(orbit) == ORBITS
    assert len(orbit[0]) == STEPS


# import numpy as np
# import astropy.units as u
# import astropy.coordinates as coord
# from astropy.coordinates import CartesianDifferential

# c= coord.SkyCoord(ra=20.*u.deg,dec=30.*u.deg,distance=2.*u.kpc,
#                 pm_ra_cosdec=-10.*u.mas/u.yr,pm_dec=20.*u.mas/u.yr,
#                 radial_velocity=50.*u.km/u.s,
#                 galcen_distance=8.*u.kpc,z_sun=15.*u.pc,
#                 galcen_v_sun=CartesianDifferential([10.0,235.,7.]*u.km/u.s))
# c.representation_type = 'cartesian'
# print(c)

# # c = coord.representation.CartesianRepresentation([1, 2, 3], [1, 2, 3], [1, 2, 3])

# # N = (2, 2)
# # xs = np.ones(N) * u.kpc
# # ys = np.ones(N) * u.kpc
# # zs = np.ones(N) * u.kpc
# # v_xs = np.zeros(N) * u.km/u.s
# # v_ys = np.linspace(100, 300, np.prod(N)).reshape(*N) * u.km/u.s
# # v_zs = np.zeros(N) * u.km/u.s
# # c = coord.SkyCoord(x=xs, y=ys, z=zs, v_x=v_xs, v_y=v_ys, v_z=v_zs, frame='galactocentric', representation_type='cartesian')
# # print(c, c.shape)


# import astropy.units as u
# import matplotlib.pyplot as plt
# import numpy as np
# import gala.integrate as gi
# import gala.dynamics as gd
# import gala.potential as gp
# from gala.units import galactic

# pot = gp.NFWPotential.from_circular_velocity(v_c=200*u.km/u.s,
#                                              r_s=10.*u.kpc,
#                                              units=galactic)

# pos = [c.x.value, c.y.value, c.z.value] * c.x.unit
# vel = [c.v_x.value, c.v_y.value, c.v_z.value] * c.v_x.unit

# ics = gd.PhaseSpacePosition(pos, vel)
# orbit = gp.Hamiltonian(pot).integrate_orbit(ics, dt=2., n_steps=10)
# # print(orbit)
# # print(type(orbit))

# # print(orbit.data.x)


# from pydgin.gala_backend import GalaBackend

# backend = GalaBackend()
# o = backend.compute_orbit(c, 0.1 * u.Gyr, 10, pot)
# print(o)

# print(len(o))
