
from .scatteringmodel import ScatteringModel

class GGADT(ScatteringModel):
    """
    Scattering model loaded from pre-computed GGADT output files.

    Attributes inherited from :ref:`scatmodels`.
    See `Hoffman & Draine (2016) <https://ui.adsabs.harvard.edu/abs/2016ApJ...817..139H>`_.
    """
    def __init__(self, from_file):
      """
      Parameters
      ----------
      from_file : str
          Name of the FITS file containing pre-computed GGADT data.
      """

      self.read_from_table(from_file)
      self.stype = 'GGADT'
      self.citation = 'https://ui.adsabs.harvard.edu/abs/2016ApJ...817..139H/abstract'

