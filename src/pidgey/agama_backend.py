import astropy.coordinates as coord
import astropy.units as u
import numpy as np

from .base import Backend


class AgamaBackend(Backend):
    def __imports__():
        import agama

    def ORBIT_TYPE(self):
        return tuple

    def _compute_orbit(self, skycoord, pot, dt, steps):
        pos = [getattr(skycoord, attr).to(u.kpc).value for attr in ("x", "y", "z")]
        vel = [
            getattr(skycoord, attr).to(u.km / u.s).value
            for attr in ("v_x", "v_y", "v_z")
        ]
        posvel = np.array([*pos, *vel]).T
        orbit = self.agama.orbit(
            potential=pot, ic=posvel, time=(dt * steps).to(u.Gyr).value, trajsize=steps
        )
        return orbit

    def _extract_points(self, orbit):
        skycoord, pot, dt, steps = self._args
        if not skycoord.shape:
            orbit = [orbit]
        xs = []
        ys = []
        zs = []
        for times, posvel in orbit:
            x, y, z, *_ = posvel.T
            xs.append(x)
            ys.append(y)
            zs.append(z)
        if not skycoord.shape:
            (xs,) = xs
            (ys,) = ys
            (zs,) = zs
        return coord.representation.CartesianRepresentation(xs, ys, zs)
