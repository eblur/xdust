from .rgscat import RGscattering
from .miescat import Mie
from .pah import PAH
from .scatteringmodel import ScatteringModel
from .make_ggadt import make_fits
from .make_ggadt_astrodust import make_fits_astrodust
from .ggadt import GGADT

"""
--------------------------------------------------------------
    API for Abstract Class 'ScatteringModel'
--------------------------------------------------------------
 A dust scattering model should contain the following attributes

calculate( lam : astropy Quantity [wavelength or energy grid]
           a   : float [grain size, micron] or astropy Quantity
           cm  : xdust.graindist.composition cm object (abstract class)
           mem_lim = : float [memory limit for scattering calculation, GB, default 1.0]
           theta = : float or np.array [radians] or astropy Quantity, angles to calculate differential scattering [default 0.0]
           **kwargs
           )

qsca : np.array, scattering efficiency [unitless]
qabs : np.array, absorption efficiency [unitless]
qext : np.array, extinction efficiency [unitless]
diff : np.array, differential scattering efficiency [ster^-1]
pars : dict, stores the parameters used to run `calculate`

write_table( outfile : string [filename for writing a FITS table of efficiency values] )
read_from_table( infile : string [filename for loading efficiency values from FITS file])
"""
