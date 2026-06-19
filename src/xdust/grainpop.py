import numpy as np
from scipy.integrate import trapezoid as trapz
import astropy.units as u

from . import graindist
from . import scatteringmodel

__all__ = ['SingleGrainPop','GrainPop','make_MRN','make_MRN_RGDrude']

MD_DEFAULT    = 1.e-4  # g cm^-2
AMIN, AMAX, P = 0.005, 0.3, 3.5  # um, um, unitless
RHO_AVG       = 3.0  # g cm^-3

# Make this a subclass of GrainDist at some point
class SingleGrainPop(graindist.GrainDist):
    """
    SingleGrainPop unifies the grain size distribution (its parent class) 
    with extinction model calculations. It provides attributes and convenience
    functions for easy access to frequently needed information. Only one dust 
    composition can be modeled at a time.

    Parameters
    ----------

    dtype : string or custom :ref:`sizedist` 
        Defines the grain radius distribution, string options: 'Grain', 'Powerlaw', 'ExpCutoff'

    cmtype : string or custom :ref:`composition` 
        Defines the optical constants and material density of the compound, string options: 'Drude', 'Silicate', 'Graphite'

    stype : string or custom :ref:`scatmodels`
        Defines what extinction model calculator to use, string options: 'RG' (Rayleigh-Gans), 'Mie' (Mie scattering)
    
    shape : string
        Only 'Sphere' is supported, some day might be used to define a custom shape

    md : float 
        Dust mass column [g cm^-2]
    
    scatm_from_file : string
        If provided, will read in a scattering model from a FITS file, and use the scattering model parameters and scattering properties from the file.
    
    **kwargs
        Extra inputs passed to ``GrainDist.__init__``.

    Note
    ----
    If an input for ``scatm_from_file`` is provided, then the ``stype`` input will be ignored.

    Attributes
    ----------
    
    lam : astropy.units.Quantity 
        wavelength or energy used for the extinction computation

    tau_sca : numpy.ndarray float 
        scattering optical depth for this grain population

    tau_abs : numpy.ndarray float 
        absorption optical depth for this grain population

    tau_ext : numpy.ndarray float 
        extinction (scattering + absorption) optical depth for this grain population

    diff : astropy.units.Quantity 
        [cm^2 ster^-1] differential scattering cross-section as a function of wavelength/energy, grain size, and angle (NE x NA x NTH)

    int_diff : astropy.units.Quantity 
        [ster^-1] differential cross-section integrated over grain size distribution, :math:`d\tau/d\Omega` (NE x NTH)

    See Also
    --------
    :ref:`GrainDist`
    :ref:`composition`
    :ref:`scatmodels`

    """
    def __init__(self, dtype, cmtype, stype, shape='Sphere', md=MD_DEFAULT, scatm_from_file=None, **kwargs):
        graindist.GrainDist.__init__(self, dtype, cmtype, shape=shape, md=md, **kwargs)

        self.lam      = None  # NE
        self.tau_sca  = None  # NE
        self.tau_abs  = None  # NE
        self.tau_ext  = None  # NE
        self.diff     = None  # NE x NA x NTH [cm^2 ster^-1]
        self.int_diff = None  # NE x NTH [ster^-1], differential xsect integrated over grain size

        # Handling scattering model FITS input, if requested
        if scatm_from_file is not None:
            self.scatm = scatteringmodel.ScatteringModel(from_file=scatm_from_file)
            assert isinstance(stype, str)
            self.scatm.stype = stype
            self.lam = self.scatm.pars['lam']
            self._calculate_tau()
        # Otherwise choose from existing (or custom) scattering calculators
        elif isinstance(stype, str):
            self._assign_scatm_from_string(stype)
        else:
            self.scatm = stype

    def _assign_scatm_from_string(self, stype):
        assert stype in ['RG', 'Mie']
        if stype == 'RG':
            self.scatm = scatteringmodel.RGscattering()
        if stype == 'Mie':
            self.scatm = scatteringmodel.Mie()

    # Run scattering model calculation, then compute optical depths
    def calculate_ext(self, lam, theta=0.0, **kwargs):
        """
        Calculate the extinction model.

        Parameters
        ----------

        lam : astropy.units.Quantity or numpy.ndarray
            Wavelength or energy values for calculating the cross-sections;
            if no units specified, defaults to keV

        theta : astropy.units.Quantity or numpy.ndarray or float
            Scattering angles for computing the differential scattering cross-section;
            if no units specified, defaults to radian.

        **kwargs
            Passed to ``self.scatm.calculate``.
        """
        self.scatm.calculate(lam, self.a, self.comp, theta=theta, **kwargs)
        self.lam      = self.scatm.pars['lam']
        self._calculate_tau()

    # Compute optical depths only
    def _calculate_tau(self):
        NE, NA, NTH = np.shape(self.scatm.diff)
        # Recall cgeo is cm^2 and ndens is cm^-2 um^-1
        # In single size grain case
        if len(self.a) == 1:
            self.tau_ext = self.ndens * self.scatm.qext[:,0] * self.cgeo
            self.tau_sca = self.ndens * self.scatm.qsca[:,0] * self.cgeo
            self.tau_abs = self.ndens * self.scatm.qabs[:,0] * self.cgeo
        # Otherwise, integrate over grain size (axis=1)
        else:
            geo_fac = self.ndens * self.cgeo  # array of length NA, unit is um^-1
            geo_2d  = np.repeat(geo_fac.reshape(1, NA), NE, axis=0)  # NE x NA
            a_um = self.a.to('micron').value
            self.tau_ext = trapz(geo_2d * self.scatm.qext, a_um, axis=1)
            self.tau_sca = trapz(geo_2d * self.scatm.qsca, a_um, axis=1)
            self.tau_abs = trapz(geo_2d * self.scatm.qabs, a_um, axis=1)

        # Recall that scatm.diff is diffrential scattering efficiency [ster^-1]
        # diff shape is NE x NA x NTH
        area_2d = np.repeat(self.cgeo.reshape(1, NA), NE, axis=0) # cm^2
        area_3d = np.repeat(area_2d.reshape(NE, NA, 1), NTH, axis=2)
        self.diff = self.scatm.diff * area_3d * u.Unit('cm^2 rad^-2') # NE x NA x NTH, [cm^2 ster^-1]

        # If a single grain size, operate in 2D (shape: NE x NTH)
        if np.size(self.a) == 1:
            int_diff = np.sum(self.scatm.diff * self.cgeo[0] * self.ndens[0], axis=1)
        # Otherwise, integrate differential scattering cross-section over NA
        else:
            a_um = self.a.to('micron').value
            agrid        = np.repeat(
                np.repeat(a_um.reshape(1, NA, 1), NE, axis=0),
                NTH, axis=2)
            ndgrid       = np.repeat(
                np.repeat(self.ndens.reshape(1, NA, 1), NE, axis=0),
                NTH, axis=2)
            int_diff = trapz(self.scatm.diff * area_3d * ndgrid, agrid, axis=1)

        self.int_diff = int_diff * u.Unit('rad^-2')  # NE x NTH, [ster^-1]

    # Plot information about the grain size distribution
    def plot_sdist(self, ax, **kwargs):
        """
        Plot information about the grain size distribution (calls GrainDist.plot)

        Parameters
        ----------

        ax : matplotlib.pyplot.axes object
        """
        self.plot(ax, **kwargs)

    # Plot the extinction properties
    def plot_ext(self, ax, keyword, unit=None, **kwargs):
        """
        Plot  the extinction properties of the grain population.

        Parameters
        ----------

        ax : matplotlib.pyplot.axes object

        keyword : str
            Extinction value(s) to plot. Options: ``'ext'``, ``'sca'``, ``'abs'``, ``'all'``.

        unit : str, optional
            Unit for the x-axis, parsable by ``astropy.units``.
            If ``None``, defaults to the unit of ``self.lam``.

        **kwargs
            Passed to ``ax.legend()``.
        """
        assert keyword in ['ext','sca','abs','all']
        try:
            assert self.lam is not None
        except:
            print("Need to run calculate_ext")
            pass

        # Handle units on the wavelength/energy scale
        xval = self.lam.value
        xunit = self.lam.unit.to_string()
        if unit is not None:
            xval = self.lam.to(unit, equivalencies=u.spectral()).value
            xunit = unit

        if keyword == 'ext':
            ax.plot(xval, self.tau_ext, **kwargs)
            ax.set_xlabel(xunit)
            ax.set_ylabel(r"$\tau_{ext}$")
        if keyword == 'sca':
            ax.plot(xval, self.tau_sca, **kwargs)
            ax.set_xlabel(xunit)
            ax.set_ylabel(r"$\tau_{sca}$")
        if keyword == 'abs':
            ax.plot(xval, self.tau_abs, **kwargs)
            ax.set_xlabel(xunit)
            ax.set_ylabel(r"$\tau_{abs}$")
        if keyword == 'all':
            ax.plot(xval, self.tau_ext, 'k-', lw=2, label='Extinction')
            ax.plot(xval, self.tau_sca, 'r--', label='Scattering')
            ax.plot(xval, self.tau_abs, 'r:', label='Absorption')
            ax.set_xlabel(xunit)
            ax.set_ylabel(r"$\tau$")
            ax.legend(**kwargs)

    # Printing information
    def info(self):
        """
        Print information about this grain population
        """
        print("Size distribution: %s" % self.size.dtype)
        print("Extinction calculated with: %s" % self.scatm.stype)
        print("Grain composition: %s" % self.comp.cmtype)
        print("rho = %.2f g cm^-3, M_d = %.2e g cm^-2" % (self.rho, self.md))

    # Write an extinction table
    def write_extinction_table(self, outfile, **kwargs):
        """
        Write the contents of the extinction calculation to a FITS file. 
        (Runs ScatteringModel.write_table)

        Parameters
        ----------

        outfile : str
            Name of file to write.

        **kwargs
            Passed to ``self.scatm.write_table``.
        """
        self.scatm.write_table(outfile, **kwargs)
        return

class GrainPop(object):
    """
    A collection of SingleGrainPop objects. Each SingleGrainPop represents 
    a different composition (or perhaps different geometries of grain size distribution). 
    Each SingleGrainPop object has a label string, and can be accessed in the same was as 
    key-value pairs are accessed in a Python dictionary.

    Parameters
    ----------
    
    gpoplist : Python list of :ref:`SingleGrainPop` objects
        Will be paired, in order, with the values in ``keys``

    keys : Python list of strings 
        Labels for each :ref:`SingleGrainPop` object.
        (Default: None, will make the keys a list of integers starting with 0)

    description : string 
        A description for this collection

    Attributes
    ----------

    keys : Python list of strings 
        Labels for each :ref:`SingleGrainPop` object. 

    gpoplist : Python list of :ref:`SingleGrainPop` objects 
        Paired, in order, with the values in ``keys``

    description : string 
        A description for this collection

    lam : astropy.units.Quantity
        Wavelength or energy used for the extinction computation
    
    md : float
        Total dust mass column [g cm^-2], sum of the `md` values from each SingleGrainPop in the collection    
    
    tau_ext : numpy.ndarray
        Total extinction optical depth as a function of wavelength / energy
    
    tau_sca : numpy.ndarray
        Total scattering optical depth as a function of wavelength / energy
    
    tau_abs : numpy.ndarray 
        Total absorption optical depth as a function of wavelength / energy
        
    See Also
    --------
    :ref:`GrainDist`
    :ref:`composition`
    :ref:`scatmodels`

    """
    def __init__(self, gpoplist, keys=None, description='Custom_GrainPopDict'):
        assert isinstance(gpoplist, list)
        if keys is None:
            self.keys = list(range(len(gpoplist)))
        else:
            self.keys = keys
        self.description = description
        self.gpoplist    = gpoplist

        # Assign the key as a description on each SingleGrainPop, helps with testing
        for k in self.keys:
            i = self.keys.index(k)
            self.gpoplist[i].description = str(self.keys[i])

        # No wavelength/energy grid assigned until the calculation is done
        self.lam = None
        
    def calculate_ext(self, lam, **kwargs):
        """
        Calculate the extinction model.

        Parameters
        ----------
        lam : astropy.units.Quantity or numpy.ndarray
            Wavelength or energy values for calculating the cross-sections;
            if no units specified, defaults to keV.

        Other Parameters
        ----------------
        theta : astropy.units.Quantity or numpy.ndarray or float
            Scattering angles for computing the differential scattering cross-section;
            if no units specified, defaults to radian.

        **kwargs
            Passed to each ``SingleGrainPop.calculate_ext``.
        """
        # Assing units if an Astropy Quantity is not input
        input_lam = lam
        if not isinstance(lam, u.Quantity):
            input_lam = lam * u.keV
        
        # Run extinction calculation on each SingleGrainPop in the list
        for gp in self.gpoplist:
            gp.calculate_ext(input_lam, **kwargs)

        # If everything went fine, store inthe input wavlength/energy grid
        self.lam = input_lam

    # Makes this object work like a dictionary
    def __getitem__(self, key):
        if isinstance(key, int):
            return self.gpoplist[key]
        else:
            assert key in self.keys
            k = self.keys.index(key)
            return self.gpoplist[k]

    @property
    def md(self):
        """
        Return the sum of the `md` values from each SingleGrainPop
        """
        result = 0.0
        for gp in self.gpoplist:
            result += gp.md
        return result

    @property
    def tau_ext(self):
        """
        Return the sum of the `tau_ext` values from each SingleGrainPop
        """
        result = 0.0
        if self.lam is None:
            print("ERROR: Extinction properties need to be calculated")
        else:
            for gp in self.gpoplist:
                result += gp.tau_ext
        return result

    @property
    def tau_sca(self):
        """
        Return the sum of the `tau_sca` values from each SingleGrainPop
        """
        result = 0.0
        if self.lam is None:
            print("ERROR: Extinction properties need to be calculated")
        else:
            for gp in self.gpoplist:
                result += gp.tau_sca
        return result

    @property
    def tau_abs(self):
        """
        Return the sum of the `tau_abs` values from each SingleGrainPop
        """
        result = 0.0
        if self.lam is None:
            print("ERROR: Extinction properties need to be calculated")
        else:
            for gp in self.gpoplist:
                result += gp.tau_abs
        return result

    def plot_ext(self, ax, keyword, unit=None, **kwargs):
        """
        Plot the sum of the extinction properties across all grain populations

        Parameters
        ----------

        ax : matplotlib.pyplot.axes object

        keyword : str
            Extinction value(s) to plot. Options: ``'ext'``, ``'sca'``, ``'abs'``, ``'all'``.

        unit : str, optional
            Unit for the x-axis, parsable by ``astropy.units``.
            If ``None``, defaults to the unit of ``self.lam``.

        **kwargs
            Passed to ``ax.legend()``.
        """
        assert keyword in ['all','ext','sca','abs']

        # Handle units on the wavelength/energy scale
        xval = self.lam.value
        xunit = self.lam.unit.to_string()
        if unit is not None:
            xval = self.lam.to(unit, equivalencies=u.spectral()).value
            xunit = unit

        if keyword == 'ext':
            ax.plot(xval, self.tau_ext, **kwargs)
            ax.set_xlabel(xunit)
            ax.set_ylabel(r"$\tau_{ext}$")
        if keyword == 'sca':
            ax.plot(xval, self.tau_sca, **kwargs)
            ax.set_xlabel(xunit)
            ax.set_ylabel(r"$\tau_{sca}$")
        if keyword == 'abs':
            ax.plot(xval, self.tau_abs, **kwargs)
            ax.set_xlabel(xunit)
            ax.set_ylabel(r"$\tau_{abs}$")
        if keyword == 'all':
            ax.plot(xval, self.tau_ext, 'k-', lw=2, label='Extinction')
            ax.plot(xval, self.tau_sca, 'r--', label='Scattering')
            ax.plot(xval, self.tau_abs, 'r:', label='Absorption')
            ax.set_xlabel(xunit)
            ax.set_ylabel(r"$\tau$")
            ax.set_title(self.description)
            ax.legend(**kwargs)

    def info(self, key=None):
        """
        Print information about this collection of grain populations

        Parameters
        ----------

        key : int or string
            The key-value for the SingleGrainPopulation of interest. 
            If None, systematically prints information for every SingleGrainPop in the collection.
        """
        if key is None:
            print("General information for %s dust grain population" % self.description)
            for gp in self.gpoplist:
                print("---")
                gp.info()
        else:
            assert key in self.keys
            self[key].info()


#---------- Basic helper functions for fast production of GrainPop objects

def make_MRN(amin=AMIN, amax=AMAX, p=P, md=MD_DEFAULT, fsil=0.6, **kwargs):
    """
    Builds a :ref:`GrainPop` object describing an MRN dust grain size distribution, 
    which is a mixture of silicate and graphite grains. Applies the 1/3 parallel, 
    2/3 perpendicular assumption of graphite grain orientations.

    Parameters
    ----------
    
    amin : astropy.units.Quantity or float
        Minimum grain radius; if a float, micron units assumed

    amax : astropy.units.Quantity or float
        Maximum grain radius; if a float, micron units assumed

    p : float 
        Powerlaw slope for grain size distribution

    md : float
        Dust mass column [g cm^-2]

    fsil : float 
        Fraction of dust mass in silicate grains

    **kwargs
        Passed to ``graindist.sizedist.Powerlaw()``.

    Returns
    -------
    GrainPop
        A :ref:`GrainPop` with keys ``'sil'``, ``'gra_para'``, and ``'gra_perp'`` for the
        silicate, graphite parallel, and graphite perpendicular grain populations, respectively.
    """
    assert isinstance(fsil, float)
    assert fsil >= 0.0 and fsil <= 1.0
    md_sil  = fsil * md
    # Graphite grain assumption: 1/3 parallel and 2/3 perpendicular
    md_gra_para = (1.0 - fsil) * md * (1.0/3.0)
    md_gra_perp = (1.0 - fsil) * md * (2.0/3.0)

    pl_sil  = graindist.sizedist.Powerlaw(amin=amin, amax=amax, p=p, **kwargs)
    pl_gra  = graindist.sizedist.Powerlaw(amin=amin, amax=amax, p=p, **kwargs)

    sil    = graindist.composition.CmSilicate()
    gra_ll = graindist.composition.CmGraphite(orient='para')
    gra_T  = graindist.composition.CmGraphite(orient='perp')

    mrn_sil = SingleGrainPop(pl_sil, sil, 'Mie', md=md_sil)
    mrn_gra_para = SingleGrainPop(pl_gra, gra_ll, 'Mie', md=md_gra_para)
    mrn_gra_perp = SingleGrainPop(pl_gra, gra_T, 'Mie', md=md_gra_perp)

    gplist = [mrn_sil, mrn_gra_para, mrn_gra_perp]
    keys   = ['sil','gra_para','gra_perp']
    return GrainPop(gplist, keys=keys, description='MRN')

def make_MRN_RGDrude(amin=AMIN, amax=AMAX, p=P, rho=RHO_AVG, md=MD_DEFAULT, **kwargs):
    """
    Builds a :ref:`SingleGrainPop` object describing an MRN dust grain size distribution, but uses 
    Rayleigh-Gans scattering with the  Drude approximation, which approximates 
    the dust grain as a sphere of free electrons.

    Parameters
    ----------

    amin : astropy.units.Quantity or float
        Minimum grain radius; if a float, micron units assumed

    amax : astropy.units.Quantity or float
        Maximum grain radius; if a float, micron units assumed

    p : float
        Powerlaw slope for grain size distribution

    rho : float
        Density of the dust grain material [g cm^-3]

    md : float
        Dust mass column [g cm^-2]

    **kwargs
        Passed to ``graindist.sizedist.Powerlaw()``.

    Returns
    -------
    SingleGrainPop
        A :ref:`SingleGrainPop` using the ``'RG'`` scattering model and ``'Drude'`` composition.
    """
    pl      = graindist.sizedist.Powerlaw(amin=amin, amax=amax, p=p, **kwargs)
    dru     = graindist.composition.CmDrude(rho=rho)
    mrn_rgd = SingleGrainPop(pl, dru, 'RG', md=md)
    return mrn_rgd
