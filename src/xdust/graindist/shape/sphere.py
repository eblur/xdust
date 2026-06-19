import numpy as np
import astropy.units as u

from . import Shape

class Sphere(Shape):
    """
    Attributes
    ----------
    shape : string : Describes the shape ('sphere')
    """
    def __init__(self):
        Shape.__init__(self)
        self.shape = 'Sphere'

    def vol(self, a):
        """
        Return the grain's volume in units of cm^3

        Parameters
        ----------
        a : astropy.units.Quantity or float
            Grain radius; plain floats are assumed to be in microns.

        Returns
        -------
        numpy.ndarray
            Grain volume :math:`\\frac{4}{3}\\pi a^3` [cm^3].
        """
        if isinstance(a, u.Quantity):
            a_cm = a.to('cm').value
        else:
            a_cm = (a * u.micron).to('cm').value
        return (4.0/3.0) * np.pi * np.power(a_cm, 3)  # cm^3

    def cgeo(self, a):
        """
        Return the geometric cross-section of a spherical particle, in units of cm^2

        Parameters
        ----------
        a : astropy.units.Quantity or float
            Grain radius; plain floats are assumed to be in microns.

        Returns
        -------
        numpy.ndarray
            Geometric cross-section :math:`\\pi a^2` [cm^2].
        """
        if isinstance(a, u.Quantity):
            a_cm = a.to('cm').value
        else:
            a_cm = (a * u.micron).to('cm').value
        return np.pi * np.power(a_cm, 2)  # cm^2
