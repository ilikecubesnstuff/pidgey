import astropy.units as u

from .base import Backend


class GalaBackend(Backend):
    def __imports__():
        from gala import dynamics, potential, units

    def ORBIT_TYPE(self):
        return self.dynamics.orbit.Orbit

    def _compute_orbit(
        self, skycoord, pot, dt, steps, pattern_speed=0 * u.km / u.s / u.kpc
    ):
        pos = [skycoord.x.value, skycoord.y.value, skycoord.z.value] * skycoord.x.unit
        vel = [
            skycoord.v_x.value,
            skycoord.v_y.value,
            skycoord.v_z.value,
        ] * skycoord.v_x.unit
        ics = self.dynamics.PhaseSpacePosition(pos, vel)

        frame = self.potential.ConstantRotatingFrame(
            Omega=[0, 0, pattern_speed.value] * pattern_speed.unit,
            units=self.units.galactic,
        )
        orbit = self.potential.Hamiltonian(pot, frame=frame).integrate_orbit(
            ics, dt=dt, n_steps=steps - 1
        )
        return orbit

    def _extract_points(self, orbit, pattern_speed=0 * u.km / u.s / u.kpc):
        frame = self.potential.ConstantRotatingFrame(
            Omega=[0, 0, pattern_speed.value] * pattern_speed.unit,
            units=self.units.galactic,
        )
        orbit = orbit.to_frame(frame)
        return orbit.data.T
