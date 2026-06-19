import numpy as np
import astropy.units as u
from scipy.integrate import trapezoid as trapz
from .. import shape
from . import Sizedist

__all__ = ['Powerlaw']

# Some default values
RHO      = 3.0     # g cm^-3 (average grain material density)

NA       = 24     # default number for grain size dist resolution
PDIST    = 3.5     # default slope for power law distribution

# min and max grain radii for MRN distribution
AMIN     = 0.005   # micron
AMAX     = 0.3     # micron

SHAPE    = shape.Sphere()

#------------------------------------

class Powerlaw(Sizedist):
    """
    A power law grain size distribution

    Parameters
    ----------
    amin : astropy.units.Quantity or float
        Minimum grain radius; plain floats are assumed to be in microns.

    amax : astropy.units.Quantity or float
        Maximum grain radius; plain floats are assumed to be in microns.

    p : float
        Power law slope for :math:`dn/da \\propto a^{-p}`.

    na : int
        Number of grain size grid points.

    log : bool
        If ``True`` (default), use log-spaced grain size grid; otherwise, use a linear grid.
    """
    def __init__(self, amin=AMIN, amax=AMAX, p=PDIST, na=NA, log=True):
        Sizedist.__init__(self)
        # Set the name of this size disribution
        self.dtype = 'Powerlaw'

        # Put amin and amax into units of micron
        if isinstance(amin, u.Quantity):
            amin_um = amin.to('micron').value
        else:
            amin_um = amin
        if isinstance(amax, u.Quantity):
            amax_um = amax.to('micron').value
        else:
            amax_um = amax
        
        # Set up the grid of grain sizes
        if log:
            self.a = np.logspace(np.log10(amin_um), np.log10(amax_um), na) * u.micron
        else:
            self.a = np.linspace(amin_um, amax_um, na) * u.micron
        
        # Power-law slope to use
        self.p    = p

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
            Column density of grains [cm^-2 um^-1].
        """
        a_um = self.a.to('micron').value

        # power law slope component
        adep  = np.power(a_um, -self.p)   # um^-p
        
        # get the mass dependence, units of g um^-p
        mgra  = shape.vol(self.a) * rho     # g (mass of each grain)
        dmda  = adep * mgra                 # g um^-p
        
        # Integrate over dmda and use that with total mass to get the 
        # correct constant for the entire function
        const = md / trapz(dmda, a_um)  # cm^-2 um^p-1
        
        # Final units are number column density per grain size unit (default:micron)
        return const * adep  # cm^-2 um^-1

    def mdens(self, md, rho=RHO, shape=SHAPE):
        """
        Calculate mass density function for the dust grains, given a total dust mass column

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
            Mass column distribution [g cm^-2 um^-1].
        """
        nd = self.ndens(md, rho, shape)  # dn/da [cm^-2 um^-1]
        mg = shape.vol(self.a) * rho     # grain mass for each radius [g]
        return nd * mg  # g cm^-2 um^-1
