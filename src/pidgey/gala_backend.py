from .base import Backend


class GalaBackend(Backend):
    def __imports__():
        from gala import dynamics, potential

    def ORBIT_TYPE(self):
        return self.dynamics.orbit.Orbit

    def _compute_orbit(self, skycoord, dt, steps, pot):
        pos = [skycoord.x.value, skycoord.y.value, skycoord.z.value] * skycoord.x.unit
        vel = [
            skycoord.v_x.value,
            skycoord.v_y.value,
            skycoord.v_z.value,
        ] * skycoord.v_x.unit
        ics = self.dynamics.PhaseSpacePosition(pos, vel)
        orbit = self.potential.Hamiltonian(pot).integrate_orbit(
            ics, dt=dt, n_steps=steps - 1
        )
        return orbit

    def _extract_points(self, orbit):
        return orbit.data.T
