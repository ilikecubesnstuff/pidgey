from abc import abstractmethod

import astropy.coordinates as coord
import astropy.units as u
from iext import ExtendImports


class Backend(ExtendImports):
    """
    Base class for orbit integration backends.

    This class defines the interface for orbit integration through
    two methods: `compute_orbit` and `get_points`. Two abstract
    methods, `_compute_orbit` and `_extract_points`, and an abstract
    attribute, `ORBIT_TYPE`, must be implemented to define a new backend.

    Attributes:
        ORBIT_TYPE (type): The type of the orbit object returned by the backend.
    """

    def __imports__():
        ...

    def __init__(self):
        self._args = None
        self._result = None

    @property
    @abstractmethod
    def ORBIT_TYPE(self):
        """
        The type of the orbit object returned by the backend.
        """
        pass

    @abstractmethod
    def _compute_orbit(
        self,
        skycoord,
        pot,
        dt,
        steps,
        pattern_speed=0 * u.km / u.s / u.kpc,
        **integration_kwargs,
    ):
        pass

    def compute_orbit(
        self,
        skycoord,
        pot,
        dt,
        steps,
        pattern_speed=0 * u.km / u.s / u.kpc,
        **integration_kwargs,
    ):
        """
        Compute an orbit for a given initial condition (through a SkyCoord object),
        through a given potential.

        Args:
            skycoord (astropy.coordinates.SkyCoord): The initial conditions of the orbit.
            pot (object): The potential to compute the orbit in.
            dt (astropy.units.Quantity): The time step for the integration.
            steps (int): The number of steps to integrate over.
            pattern_speed (astropy.units.Quantity): The pattern speed of the potential.
            **integration_kwargs: Additional arguments for integration routines.

        Returns:
            astropy.coordinates.SkyCoord: The points from the orbit integration.
        """
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

        with u.add_enabled_equivalencies(u.dimensionless_angles()):
            self._result = self._compute_orbit(
                skycoord, pot, dt, steps, pattern_speed, **integration_kwargs
            )
            return self._extract_points(self._result, pattern_speed)

    @abstractmethod
    def _extract_points(self, orbit, pattern_speed=0 * u.km / u.s / u.kpc):
        pass

    def get_points(self, pattern_speed=0 * u.km / u.s / u.kpc):
        """
        Get the points from the orbit integration result in a SkyCoord object.

        Args:
            pattern_speed (astropy.units.Quantity): The pattern speed of the potential.

        Returns:
            astropy.coordinates.SkyCoord: The points from the orbit integration.
        """
        if not isinstance(self._result, self.ORBIT_TYPE):
            raise TypeError("No computed orbit to retrive points from.")

        with u.add_enabled_equivalencies(u.dimensionless_angles()):
            return self._extract_points(self._result, pattern_speed)
