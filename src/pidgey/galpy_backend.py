import astropy.coordinates as coord
import astropy.units as u
import numpy as np
from galpy.util.conversion import get_physical

from .base import Backend


class GalpyBackend(Backend):
    def __imports__():
        from galpy import orbit

    def ORBIT_TYPE(self):
        return self.orbit.Orbit

    def _compute_orbit(
        self, skycoord, pot, dt, steps, pattern_speed=0 * u.km / u.s / u.kpc
    ):
        # unit consistency issues... look into this later
        pot_units = get_physical(pot)
        orbit = self.orbit.Orbit(skycoord, **pot_units)

        t = np.arange(steps) * dt
        orbit.integrate(t, pot)
        return orbit

    def _extract_points(self, orbit, pattern_speed=0 * u.km / u.s / u.kpc):
        skycoord, pot, dt, steps = self._args
        t = np.arange(steps) * dt
        phi_offset = t * pattern_speed

        R = orbit.R(t)
        phi = orbit.phi(t) + phi_offset.to(u.dimensionless_unscaled)

        x = R * np.cos(phi.value)
        y = R * np.sin(phi.value)
        z = orbit.z(t)

        return coord.representation.CartesianRepresentation(x, y, z)
