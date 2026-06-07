import numpy as np
import astropy.units as u

from .. import shape

__all__ = ['Grain']

# Some default values
AMICRON = 1.0 * u.micron
RHO     = 3.0  # g cm^-3 (average grain material density)

SHAPE    = shape.Sphere()

#-------------------------------

class Grain(object):
    """
    A single grain size distribution

    Attributes
    ----------
    a : astropy.units.Quantity
        Grain radius.

    dtype : str
        ``'Grain'``.
    """
    def __init__(self, rad=AMICRON):
        """
        Parameters
        ----------
        rad : astropy.units.Quantity or float
            Grain radius; plain floats are assumed to be in microns.
        """
        self.dtype = 'Grain'
        assert np.size(rad) == 1
        if isinstance(rad, u.Quantity):
            self.a = np.array([rad.value]) * rad.unit
        else:
            self.a = np.array([rad]) * u.micron

    def ndens(self, md, rho=RHO, shape=SHAPE):
        """
        Calculate number density of dust grains, given a dust mass column

        Parameters
        ----------
        md : float
            Mass column density [g cm^-2].

        rho : float
            Grain material density [g cm^-3].

        shape : xdust.graindist.shape object
            Grain shape (default: ``Sphere``).

        Returns
        -------
        numpy.ndarray
            Column density of grains [cm^-2].
        """
        gvol = shape.vol(self.a) # cm^3
        return md / (gvol * rho)  # cm^-2

    def mdens(self, md, rho=RHO, shape=SHAPE):
        """
        Calculate dust mass distribution -- in this case, no calculation needed
        because it is a single grain size.

        Why do I have this function here? Mostly to preserve interoperability. 
        Other distributions return a continuous function that can be integrated.
        """
        return md
