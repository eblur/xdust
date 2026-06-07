import numpy as np
import astropy.units as u
import astropy.constants as c
from scipy.integrate import trapezoid as trapz
from scipy.special import erf
from .. import shape

__all__ = ['WD01']

# Some default values
RHO_C    = 2.24    # g cm^-3 (graphite grain material density), p.298
RHO_S    = 3.5     # g cm^-3 (silicate grain material density), p.298

NA       = 100     # default number for grain size dist resolution

# min and max grain radii
AMIN     = 0.00035   # micron
AMAX     = 0.5    # micron, relevant for silicates
# Silicate grain models mostly max out at 0.5 micron
# Carbonaceous grain models can go out to 10 micron for RV=5.5, 1 micron for RV=3.1

SHAPE    = shape.Sphere()

# I will take cases only where 10^5 bc = 0 because with current instruments
# we don't have a lot of constraints on small graphitic grains
class _WD01_params(object):
    def __init__(self, table1):
        self.RV       = table1['RV-case'] # Note the RV value and "Case" for solving
        # Case A - assumes A(I)/NH is independent of RV (makes it so that less grain volume can account for high RV)
        # Case B - assumes grain volume is held fixed to value for RV=3.1
        # Case B is more physical (it is unlikely that grains are removed from solid phase in denser regions, where RV is higher)
        self.bc       = table1['bc']
        self.alpha_g  = table1['alpha_g']
        self.beta_g   = table1['beta_g']
        self.a_tg     = table1['a_tg'] # micron
        self.a_cg     = table1['a_cg'] # micron
        self.C_g      = table1['C_g']
        self.alpha_s  = table1['alpha_s']
        self.beta_s   = table1['beta_s']
        self.a_ts     = table1['a_ts'] # micron
        self.C_s      = table1['C_s']
        self.a_cs     = 0.1 # micron, see Table 1 of the paper
    
    def silicate_function(self, agrid):
        """
        Returns shape of silicate grain size distribution function [um^-1]
        agrid : numpy ndarray of grain sizes in micron
        """
        F = 0.0 # default, if not set properly, the final result will be all zeros
        if self.beta_s >= 0:
            F = 1 + self.beta_s * (agrid / self.a_ts)
        elif self.beta_s < 0:
            F = 1 / (1 - self.beta_s * (agrid / self.a_ts))

        result = (self.C_s / agrid) * (agrid / self.a_ts)**self.alpha_s * F

        ii = (agrid > self.a_ts)
        result[ii] *= np.exp(-((agrid[ii] - self.a_ts) / self.a_cs)**3)
        return result
    
    def graphite_function(self, agrid):
        """
        Returns shape of graphite grain size distribution function [um^-1]
        agrid : numpy ndarray of grain sizes in micron
        """
        m_c = 12.011 * c.u.to('g').value
        B1 = (3 / (2 * np.pi)**1.5)
        B1 *= np.exp(-4.5 * 0.4**2) / (2.24 * 3.5**3 * u.Angstrom.to('cm')**3)
        B1 *= 0.75 * self.bc * m_c
        B1 /= 1 + erf(3 * 0.4 / np.sqrt(2) + np.log(1.0)/(0.4 * np.sqrt(2)))
        B2 = (3 / (2 * np.pi)**1.5)
        B2 *= np.exp(-4.5 * 0.4**2) / (2.24 * 30.**3 * u.Angstrom.to('cm')**3)
        B2 *= 0.25 * self.bc * m_c
        B2 /= 1 + erf(3 * 0.4 / np.sqrt(2) + np.log(30./3.5)/(0.4 * np.sqrt(2)))
        D = (B1 / agrid) * np.exp(-0.5 * (np.log(agrid / 0.00035) / 0.4)**2)
        D += (B2 / agrid) * np.exp(-0.5 * (np.log(agrid / 0.00300) / 0.4)**2)

        F = 0.0 # default, if not set properly, the final result will be all zeros
        if self.beta_g >= 0:
            F = 1 + self.beta_g * (agrid / self.a_tg)
        elif self.beta_g < 0:
            F = 1 / (1 - self.beta_g * (agrid / self.a_tg))
        
        result = F * (self.C_g / agrid) * (agrid / self.a_tg)**self.alpha_g

        ii = (agrid > self.a_tg)
        result[ii] *= np.exp(-((agrid[ii] - self.a_tg) / self.a_cg)**3)

        result += D
        return result

# bc = 0 cases
# MW_RV3p1_params = _WD01_params({
#     'RV-case': '3.1-A',
#     'bc'     : 0.0,
#     'alpha_g': -2.25,
#     'beta_g' : -0.0648,
#     'a_tg'   : 0.00745,
#     'a_cg'   : 0.606,
#     'C_g'    : 9.94e-11,
#     'alpha_s': -1.48,
#     'beta_s' : -9.34,
#     'a_ts'   : 0.172,
#     'C_s'    : 1.02e-12
# })

# MW_RV4p0_params = _WD01_params({
#     'RV-case': '4.0-B',
#     'bc'     : 0.0,
#     'alpha_g': -2.62,
#     'beta_g' : -0.0144,
#     'a_tg'   : 0.0187,
#     'a_cg'   : 5.74,
#     'C_g'    : 6.46e-12,
#     'alpha_s': -2.01,
#     'beta_s' : 0.894,
#     'a_ts'   : 0.198,
#     'C_s'    : 4.95e-14
# })

# MW_RV5p5_params = _WD01_params({
#     'RV-case': '5.5-B',
#     'bc'     : 0.0,
#     'alpha_g': -2.80,
#     'beta_g' : 0.0356,
#     'a_tg'   : 0.0203,
#     'a_cg'   : 3.43,
#     'C_g'    : 2.74e-12,
#     'alpha_s': -1.09,
#     'beta_s' : -0.370,
#     'a_ts'   : 0.218,
#     'C_s'    : 1.17e-13
# })


# distribution favored by WD01
MW_RV3p1_params = _WD01_params({
    'RV-case': '3.1-A',
    'bc'     : 6.0e-5,
    'alpha_g': -1.54,
    'beta_g' : -0.165,
    'a_tg'   : 0.0107,
    'a_cg'   : 0.428,
    'C_g'    : 9.99e-12,
    'alpha_s': -2.21,
    'beta_s' : 0.300,
    'a_ts'   : 0.164,
    'C_s'    : 1.00e-13
})

# distribution favored by WD01
MW_RV4p0_params = _WD01_params({
    'RV-case': '4.0-B',
    'bc'     : 4.0e-5,
    'alpha_g': -1.96,
    'beta_g' : -0.813,
    'a_tg'   : 0.0693,
    'a_cg'   : 3.48,
    'C_g'    : 2.95e-13,
    'alpha_s': -2.11,
    'beta_s' : 2.10,
    'a_ts'   : 0.198,
    'C_s'    : 3.13e-14
})

MW_RV5p5_params = _WD01_params({
    'RV-case': '5.5-B',
    'bc'     : 3.0e-5,
    'alpha_g': -1.90,
    'beta_g' : -0.0517,
    'a_tg'   : 0.0120,
    'a_cg'   : 7.28,
    'C_g'    : 2.86e-12,
    'alpha_s': -1.13,
    'beta_s' : -0.109,
    'a_ts'   : 0.211,
    'C_s'    : 1.04e-13
})



#------------------------------------

class WD01(object):
    """
    Weingartner & Draine (2001) grain size distribution
    """
    def __init__(self, galaxy='MW', RV=3.1, form='silicate', amin=AMIN, amax=AMAX, na=NA, log=True):
        """
        Parameters
        ----------
        galaxy : str
            Name of the galaxy. Default ``'MW'``.

        RV : float
            Total-to-selective extinction ratio. Default ``3.1``.

        form : str
            Dust grain composition; ``'silicate'`` or ``'graphite'``.

        amin : astropy.units.Quantity or float
            Minimum grain radius; plain floats are assumed to be in microns.

        amax : astropy.units.Quantity or float
            Maximum grain radius; plain floats are assumed to be in microns.

        na : int
            Number of grain size bins. Default ``100``.

        log : bool
            If ``True``, use log-spaced grain size grid. Default ``True``.
        """
        self.dtype = f'WD01-{galaxy}-{RV}-{form}'

        # Put amin and acut into units of micron
        if isinstance(amin, u.Quantity):
            amin_um = amin.to('micron').value
        else:
            amin_um = amin
        if isinstance(amax, u.Quantity):
            amax_um = amax.to('micron').value
        else:
            amax_um = amax

        # Set up the grid of grain sizes
        self.a = None
        if log:
            self.a = np.logspace(np.log10(amin_um), np.log10(amax_um), na) * u.micron
        else:
            self.a = np.linspace(amin_um, amax_um, na) * u.micron

        # Set up the relevant parameters
        self.params = None
        if galaxy == 'MW':
            if RV == 3.1:
                self.params = MW_RV3p1_params
            elif RV == 4.0:
                self.params = MW_RV4p0_params
            elif RV == 5.5:
                self.params = MW_RV5p5_params
            else:
                raise ValueError(f"RV={RV} not implemented for galaxy={galaxy}")
        
        self.form = None
        if form.lower() in ['silicate', 'silicates']:
            self.form = 'silicate'
        elif form.lower() in ['carbon', 'carbonaceous', 'graphite', 'graphitic']:
            self.form = 'graphite'
        else:
            raise ValueError(f"form={form} not recognized; choose from 'silicate' or 'graphite'")
        
        self.rho = None
        if self.form == 'silicate':
            self.rho = RHO_S
        elif self.form == 'graphite':
            self.rho = RHO_C

    def ndens(self, md, rho=None, shape=SHAPE):
        """
        Calculate number density of dust grains, given a dust mass column

        Parameters
        ----------
        md : float
            Mass column density [g cm^-2].

        rho : float, optional
            Grain material density [g cm^-3]. If provided and different from
            ``self.rho``, overrides the default density.

        shape : xdust.graindist.shape object
            Grain shape (default: ``Sphere``).

        Returns
        -------
        numpy.ndarray
            Column density of grains [cm^-2 um^-1].
        """
        a_um = self.a.to('micron').value

        adep = None # unit will be um^-1
        if self.form == 'silicate':
            adep = self.params.silicate_function(a_um)
        elif self.form == 'graphite':
            adep = self.params.graphite_function(a_um)

        if (rho is not None) and (rho != self.rho):
            print(f"Overriding default grain density of {self.rho} g cm^-3 with {rho} g cm^-3")
            self.rho = rho

        # get the mass dependence, units of g um^-1
        mgra  = shape.vol(self.a) * self.rho  # g (mass of each grain)
        dmda  = adep * mgra

        # integrate to get the correct scaling constant
        const = md / trapz(dmda, a_um)  # g cm^-2 / g = cm^-2

        # Final units are number column density per grain size unit
        return const * adep  # cm^-2 um^-1

    def mdens(self, md, rho=None, shape=SHAPE):
        """
        Calculate mass density function for the dust grains, given a total dust mass column

        Parameters
        ----------
        md : float
            Mass column density [g cm^-2].

        rho : float, optional
            Grain material density [g cm^-3]. If provided and different from
            ``self.rho``, overrides the default density.

        shape : xdust.graindist.shape object
            Grain shape (default: ``Sphere``).

        Returns
        -------
        numpy.ndarray
            Mass column distribution [g cm^-2 um^-1].
        """
        if (rho is not None) and (rho != self.rho):
            print(f"Overriding default grain density of {self.rho} g cm^-3 with {rho} g cm^-3")
            self.rho = rho

        nd = self.ndens(md, self.rho, shape)  # dn/da [cm^-2 um^-1]
        mg = shape.vol(self.a) * self.rho     # grain mass for each radius [g]
        return nd * mg  # g cm^-2 um^-1



def WD01_validation_plot(galaxy='MW', RV=3.1, amax=AMAX):
    """
    Make a plot to validate the WD01 implementation
    """
    import matplotlib.pyplot as plt

    agrid = np.logspace(np.log10(0.00035), np.log10(amax), 200) # micron
    # agrid = np.logspace(np.log10(0.00035), np.log10(amax), 20) # for debugging

    if galaxy == 'MW' and RV == 3.1:
        params = MW_RV3p1_params
    elif galaxy == 'MW' and RV == 4.0:
        params = MW_RV4p0_params
    elif galaxy == 'MW' and RV == 5.5:
        params = MW_RV5p5_params
    else:
        raise ValueError(f"RV={RV} not implemented for galaxy={galaxy}")
    
    dnda_sil = params.silicate_function(agrid)
    dnda_gra = params.graphite_function(agrid)

    plt.figure()
    plt.loglog(agrid, agrid**4 * u.micron.to('cm')**3 * 1e29 * dnda_gra, label=f'{galaxy} RV={RV} graphite')
    plt.loglog(agrid, agrid**4 * u.micron.to('cm')**3 * 1e29 * dnda_sil, label=f'{galaxy} RV={RV} silicate')
    plt.xlabel('Grain radius [micron]')
    plt.ylabel('$10^{29} n_H^{-1} a^4 dn/da$ [cm$^{3}$]')
    plt.title(f'Weingartner & Draine (2001) distribution')
    plt.loglog()
    plt.legend()
