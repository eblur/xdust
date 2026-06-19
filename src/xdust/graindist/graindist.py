import numpy as np
from . import sizedist
from . import composition
from . import shape as sh

MD_DEFAULT = 1.e-4  # g cm^-2
RHO        = 3.0    # g cm^-3
AMAX       = 0.3  # um

ALLOWED_SIZES  = ['Grain','Powerlaw','ExpCutoff','Astrodust','WD01']
ALLOWED_COMPS  = ['Drude', 'Silicate', 'Graphite']
SHAPES = {'Sphere':sh.Sphere()}

__all__ = ['GrainDist']

class GrainDist(object):
    """
    Graindist ties together the size distribution (which has no set column density) 
    and composition (which contains optical constants and density) and anchors them
    with dust mass density and a series of convenience functions to get at some of the
    frequently needed information.

    Parameters
    ----------
    dtype : str or xdust.graindist.sizedist object
        Grain radius distribution. String options: ``'Grain'``, ``'Powerlaw'``,
        ``'ExpCutoff'``, ``'Astrodust'``, ``'WD01'``.

    cmtype : str or xdust.graindist.composition object
        Optical constants and compound density. String options: ``'Drude'``,
        ``'Silicate'``, ``'Graphite'``.

    shape : str
        Grain shape. ``'Sphere'`` is the only supported option.

    md : float
        Dust mass column [g cm^-2].

    rho : float, optional
        If provided, passed to the ``rho`` keyword of the composition object.

    **kwargs
        Extra inputs passed to ``sizedist.__init__``.
    
    Attributes
    ----------
    size : :ref:`sizedist` object
        Grain size distribution.

    comp : :ref:`composition` object
        Optical constants and material density.

    shape : :ref:`shape` object
        Grain shape object.

    md : float
        Dust mass column density [g cm^-2].
    """
    def __init__(self, dtype, cmtype, shape='Sphere', md=MD_DEFAULT, rho=None, **kwargs):
        self.md = md

        if isinstance(dtype, str):
            self._assign_sizedist_from_string(dtype, **kwargs)
        else:
            self.size = dtype

        if isinstance(cmtype, str):
            self._assign_comp_from_string(cmtype, rho)
        else:
            self.comp = cmtype

        if isinstance(shape, str):
            self._assign_shape_from_string(shape)
        else:
            self.shape = shape

    @property
    def a(self):
        """
        Grain radius grid

        Returns
        -------
        astropy.units.Quantity
        """
        return self.size.a

    @property
    def ndens(self):
        """
        Number column density of dust grains as a function of radius [cm^-2 um^-1], using the dust mass column density in ``self.md`` and the grain material density in ``self.comp.rho``.
        
        Returns
        -------
        numpy.ndarray
        """
        return self.size.ndens(self.md, rho=self.comp.rho, shape=self.shape)

    @property
    def mdens(self):
        """
        Mass column density of dust grains as a function of radius [g cm^-2 um^-1], using the dust mass column density in ``self.md`` and the grain material density in ``self.comp.rho``.

        Returns
        -------
        numpy.ndarray
        """
        mg = self.shape.vol(self.a) * self.comp.rho  # mass of each dust grain [g]
        return self.ndens * mg

    @property
    def rho(self):
        """
        Grain material density [g cm^-3]

        Returns
        -------
        float
        """
        return self.comp.rho

    @property
    def cgeo(self):
        """
        Geometric cross section of the grains [float, cm^2].

        Returns
        -------
        numpy.ndarray    
        """
        return self.shape.cgeo(self.a)

    @property
    def vol(self):
        """
        Volume of the grains [float, cm^3].
        
        Returns
        -------
        numpy.ndarray
        """
        return self.shape.vol(self.a)

    def plot(self, ax, **kwargs):
        """
        Plots the grain size distribution on a given matplotlib axis        

        Parameters
        ----------
        ax : matplotlib.axes.Axes
            Axes on which to plot the grain size distribution
        """
        ax.plot(self.a.to('micron').value, self.ndens * np.power(self.a.to('micron').value, 4), **kwargs)
        ax.set_xlabel("Radius (micron)")
        ax.set_ylabel("$(dn/da) a^4$ (cm$^{-2}$ um$^{3}$)")
        ax.set_xscale('log')
        ax.set_yscale('log')
        return

    def _assign_sizedist_from_string(self, dtype, **kwargs):
        assert dtype in ALLOWED_SIZES
        if dtype == 'Grain':
            self.size = sizedist.Grain(**kwargs)
        if dtype == 'Powerlaw':
            self.size = sizedist.Powerlaw(**kwargs)
        if dtype == 'ExpCutoff':
            self.size = sizedist.ExpCutoff(**kwargs)
        if dtype == 'Astrodust':
            self.size = sizedist.Astrodust(**kwargs)
        if dtype == 'WD01':
            self.size = sizedist.WD01(**kwargs)
        return

    def _assign_comp_from_string(self, cmtype, rho):
        assert cmtype in ALLOWED_COMPS
        if cmtype == 'Drude':
            if rho is not None: self.comp = composition.CmDrude(rho=rho)
            else: self.comp = composition.CmDrude()
        if cmtype == 'Silicate':
            if rho is not None: self.comp = composition.CmSilicate(rho=rho)
            else: self.comp = composition.CmSilicate()
        if cmtype == 'Graphite':
            if rho is not None: self.comp = composition.CmGraphite(rho=rho)
            else: self.comp = composition.CmGraphite()
        return

    def _assign_shape_from_string(self, shape):
        assert shape in SHAPES.keys()
        self.shape = SHAPES[shape]
        return
