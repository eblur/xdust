"""
-------
API
-------

*Sizedist* must contain the following attributes

`dtype` : a string description

`a` : an array

And must contain the following two methods:

`ndens` (md, rho, shape) returns number density of dust grains [e.g. cm^-2 um^-1]

`mdens` (md, rho, shape) returns mass density of dust grains [e.g. g cm^-2 um^-1]
"""

class Sizedist(object):
    """
    Abstract class for grain size distributions

    Attributes
    ----------
    dtype : str
        String description of the size distribution.

    a : astropy.units.Quantity
        Grain radius grid.
    """
    def __init__(self):
        self.dtype = 'Abstract Sizedist'
        self.a = None

    def ndens(self, md, rho, shape):
        """
        Calculate number density of dust grains, given a dust mass column

        Parameters
        ----------
        md : float
            Mass column density [g cm^-2].

        rho : float
            Grain material density [g cm^-3].

        shape : xdust.graindist.shape object
            Grain shape (default: ``Sphere``, the only supported option)

        Returns
        -------
        numpy.ndarray
            Column density of grains [cm^-2 um^-1].
        """
        raise NotImplementedError('ndens method not implemented for Sizedist')   

    def mdens(self, md, rho, shape):
        """
        Calculate mass density function for the dust grains, given a total dust mass column

        Parameters
        ----------
        md : float
            Mass column density [g cm^-2].

        rho : float
            Grain material density [g cm^-3].

        shape : xdust.graindist.shape object
            Grain shape (default: ``Sphere``, the only supported option)

        Returns
        -------
        numpy.ndarray
            Mass column distribution [g cm^-2 um^-1].
        """
        raise NotImplementedError('mdens method not implemented for Sizedist') 

from .grain import Grain
from .powerlaw import Powerlaw
from .exp_cutoff import ExpCutoff
from .astrodust import Astrodust
from .weingartner import WD01

from .. import shape
