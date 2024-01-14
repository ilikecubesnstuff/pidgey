import astropy.coordinates as coord
import numpy as np

from .base import Backend


class GalpyBackend(Backend):
    def __imports__():
        from galpy import orbit

    def ORBIT_TYPE(self):
        return self.orbit.Orbit

    def _compute_orbit(self, skycoord, dt, steps, pot):
        orbit = self.orbit.Orbit(skycoord)
        t = np.arange(steps) * dt
        orbit.integrate(t, pot)
        return orbit

    def _extract_points(self, orbit):
        skycoord, dt, steps, pot = self._args
        t = np.arange(steps) * dt
        return coord.representation.CartesianRepresentation(
            orbit.x(t), orbit.y(t), orbit.z(t)
        )
