
xdust.scatteringmodel
=====================

Abstract Class: *ScatteringModel*
---------------------------------

A dust scattering model must contain the following:

**Attributes:**

- ``stype`` -- a string describing the model

- ``cite`` -- a citation for the model--

- ``qsca`` -- unitless scattering efficiency

- ``qext`` -- unit-less extinction efficiency

- ``qabs`` -- unit-less absorption efficiency

- ``diff`` -- differential scattering efficiency, ster^-1

- ``pars`` -- dict, stores the parameters used to run `calculate`

**Methods**

- ``calculate(lam, a, cm, mem_lim=8.0, theta=0.0)`` 
    calculates the extinction 
    efficiencies and differential scattering cross section for the given 
    wavelength/energy grid, grain size, and composition.  
    Stores the results in the attributes `qsca`, `qext`, `qabs`, and `diff`. 
    The parameters used to run `calculate` are stored in the attribute `pars` 
    as a dictionary.

           ``lam`` -- astropy Quantity [wavelength or energy grid]

           ``a``   -- float [grain size, micron] or astropy Quantity 

           ``cm `` -- xdust.graindist.composition object

           ``mem_lim`` -- float [memory limit for scattering calculation, GB, default 8.0]

           ``theta`` -- float or np.array [radians] or astropy Quantity, angles to calculate differential scattering [default 0.0]

- ``write_table(outfile)`` writes a FITS table of efficiency values to the specified filename   
            
            ``outfile`` -- string [filename for writing a FITS table of efficiency values]

- ``read_from_table(infile)`` reads efficiency values from a FITS table and stores them in the attributes `qsca`, `qext`, `qabs`, and `diff`. The parameters used to run the scattering model are stored in the attribute `pars` as a dictionary. 

            ``infile`` -- string [filename for loading efficiency values from FITS file]


RGscattering
------------
.. autoclass:: xdust.scatteringmodel.RGscattering

Mie
---
.. autoclass:: xdust.scatteringmodel.Mie

GGADT
-----
.. autoclass:: xdust.scatteringmodel.GGADT
