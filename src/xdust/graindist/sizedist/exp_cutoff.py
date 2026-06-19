import numpy as np
import astropy.units as u
from scipy.integrate import trapezoid as trapz
from . import Sizedist
from .. import shape

__all__ = ['ExpCutoff']

# Some default values
RHO      = 3.0     # g cm^-3 (average grain material density)

NA       = 24      # default number for grain size dist resolution
PDIST    = 3.5     # default slope for power law distribution

# min and max grain radii for MRN distribution
AMIN     = 0.005   # micron
ACUT     = 0.1     # micron
NFOLD    = 5       # Number of e-foldings (a/amax) to cover past the amax point

SHAPE    = shape.Sphere()

#------------------------------------

class ExpCutoff(object):
    """
    Power law grain size distribution with an exponential cut-off at the large end

    Parameters
    ----------
    amin : astropy.units.Quantity or float
        Minimum grain radius; plain floats are assumed to be in microns.

    acut : astropy.units.Quantity or float
        Exponential cut-off grain radius; plain floats are assumed to be in microns.

    p : float
        Power law slope for :math:`dn/da \\propto a^{-p}`.

    na : int
        Number of grain size grid points.

    log : bool
        If ``True`` (default), use log-spaced grain size grid; otherwise, use a linear grid.

    nfold : int
        Number of e-foldings past ``acut`` to extend the grid.
    """
    def __init__(self, amin=AMIN, acut=ACUT, p=PDIST, na=NA, log=True, nfold=NFOLD):
        Sizedist.__init__(self)
        self.dtype = 'ExpCutoff'

        # Put amin and acut into units of micron
        if isinstance(amin, u.Quantity):
            amin_um = amin.to('micron').value
        else:
            amin_um = amin
        if isinstance(acut, u.Quantity):
            acut_um = acut.to('micron').value
        else:
            acut_um = acut

        # Set up the grid of grain sizes
        if log:
            self.a = np.logspace(np.log10(amin_um), np.log10(acut_um * nfold), na) * u.micron
        else:
            self.a = np.linspace(amin_um, acut_um * nfold, na) * u.micron

        # Log the relevant params
        self.p    = p
        self.acut = acut_um * u.micron

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
        acut_um = self.acut.to('micron').value

        # power law slope component
        adep  = np.power(a_um, -self.p) * np.exp(-a_um/acut_um)   # um^-p

        # get the mass dependence, units of g um^-p
        mgra  = shape.vol(self.a) * rho  # g (mass of each grain)
        dmda  = adep * mgra              # g um^-p

        # integrate to get the correct scaling constant
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
