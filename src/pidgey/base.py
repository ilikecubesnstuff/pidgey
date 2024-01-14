from abc import abstractmethod, abstractproperty

import astropy.coordinates as coord
import astropy.units as u
from iext import ExtendImports


class Backend(ExtendImports):
    def __imports__():
        ...

    def __init__(self):
        self._args = None
        self._result = None

    @abstractproperty
    def ORBIT_TYPE(self):
        pass

    @abstractmethod
    def _compute_orbit(self, skycoord, pot, dt, steps):
        pass

    def compute_orbit(self, skycoord, pot, dt, steps):
        if not isinstance(skycoord, coord.SkyCoord):
            raise TypeError(
                "coord must be passed in as a astropy.coordinates.SkyCoord object."
                f"{skycoord} is not a astropy.coordinates.SkyCoord object."
            )
        skycoord.representation_type = "cartesian"
        if not (
            hasattr(skycoord, "v_x")
            and hasattr(skycoord, "v_y")
            and hasattr(skycoord, "v_z")
        ):
            raise TypeError(
                "coord must have 3-D velocity information to compute an orbit."
            )
        if not isinstance(dt, u.Quantity):
            raise TypeError(
                "dt must be passed in as a astropy.units.Quantity object."
                f"{dt} is not a astropy.coordinates.SkyCoord object."
            )
        self._args = skycoord, pot, dt, steps
        self._result = self._compute_orbit(skycoord, pot, dt, steps)
        return self._extract_points(self._result)

    @abstractmethod
    def _extract_points(self, orbit):
        pass

    def get_points(self):
        if not isinstance(self._result, self.ORBIT_TYPE):
            raise TypeError("No computed orbit to retrive points from.")
        return self._extract_points(self._result)
