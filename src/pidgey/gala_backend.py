import astropy.coordinates as coord
import astropy.units as u
import numpy as np

from .base import Backend


class GalaBackend(Backend):
    def __imports__():
        from gala import dynamics, potential, units

    @property
    def ORBIT_TYPE(self):
        return self.dynamics.orbit.Orbit

    def _compute_orbit(
        self,
        skycoord,
        pot,
        dt,
        steps,
        pattern_speed=0 * u.km / u.s / u.kpc,
        integrator=None,
        integrator_kwargs=None,
        cython_if_possible=True,
        store_all=True,
        **time_spec,
    ):
        pos = [skycoord.x.value, skycoord.y.value, skycoord.z.value] * skycoord.x.unit
        vel = [
            skycoord.v_x.value,
            skycoord.v_y.value,
            skycoord.v_z.value,
        ] * skycoord.v_x.unit
        ics = self.dynamics.PhaseSpacePosition(pos, vel)

        if integrator_kwargs is None:
            integrator_kwargs = {}
        orbit = pot.integrate_orbit(
            ics,
            integrator,
            integrator_kwargs,
            cython_if_possible,
            store_all,
            dt=dt,
            n_steps=steps - 1,
        )

        # In the case of a single element collection of initial conditions,
        # the orbit is returned as a 1D array. We need to reshape it to the
        # expected shape.
        expected_shape = (steps, *skycoord.shape)
        if orbit.shape != expected_shape:
            orbit = orbit.reshape(expected_shape)
        return orbit

    def _extract_points(self, orbit, pattern_speed=0 * u.km / u.s / u.kpc):
        skycoord, pot, dt, steps = self._args
        _, *shape = orbit.shape

        t = np.arange(steps) * dt
        t = (
            np.array([t for _ in range(int(np.prod(shape)))]).reshape(orbit.shape)
            * dt.unit
        )
        phi_offset = t * pattern_speed

        orbit = orbit.represent_as(coord.CylindricalRepresentation)
        rho = orbit.data.rho
        phi = orbit.data.phi + phi_offset.to(u.dimensionless_unscaled)

        x = rho * np.cos(phi.value)
        y = rho * np.sin(phi.value)
        z = orbit.data.z

        x = np.moveaxis(x, 0, -1)
        y = np.moveaxis(y, 0, -1)
        z = np.moveaxis(z, 0, -1)
        return coord.representation.CartesianRepresentation(x, y, z)
