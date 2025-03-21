import astropy.coordinates as coord
import astropy.units as u
import numpy as np

from .base import Backend


class GalpyBackend(Backend):
    def __imports__():
        from galpy import orbit
        from galpy.util import conversion

    @property
    def ORBIT_TYPE(self):
        return self.orbit.Orbit

    def _compute_orbit(
        self,
        skycoord,
        pot,
        dt,
        steps,
        pattern_speed=0 * u.km / u.s / u.kpc,
        **integration_kwargs,
    ):
        rovo = self.conversion.get_physical(pot)
        ro = rovo["ro"]
        vo = rovo["vo"]
        z_sun = skycoord.z_sun
        skycoord = coord.SkyCoord(
            x=skycoord.galactocentric.x,
            y=skycoord.galactocentric.y,
            z=skycoord.galactocentric.z,
            v_x=skycoord.galactocentric.v_x,
            v_y=skycoord.galactocentric.v_y,
            v_z=skycoord.galactocentric.v_z,
            frame="galactocentric",
            representation_type="cartesian",
            galcen_distance=np.sqrt(ro**2 * u.kpc**2 + z_sun**2),
            galcen_v_sun=[vo, 0.0, 0.0] * u.km / u.s,
            z_sun=z_sun,
        )
        orbit = self.orbit.Orbit(skycoord, **rovo)

        t = np.arange(steps) * dt
        orbit.integrate(t, pot, **integration_kwargs)
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
