
class Shape:
    """
    Abstract class for grain shape.
    """
    def __init__(self):
        self.shape = 'Undefined'
    
    def vol(self, a):
        """
        Volume of the grain as a function of grain size.

        Parameters
        ----------
        a : astropy.units.Quantity or float
            Grain size; plain floats are assumed to be in microns.

        Returns
        -------
        numpy.ndarray
            Volume of the grain [cm^-3]
        """
        raise NotImplementedError("Subclasses must implement this method.")

    def cgeo(self, a):
        """
        Geometric cross-section of the grain as a function of grain size.

        Parameters
        ----------
        a : astropy.units.Quantity or float
            Grain size; plain floats are assumed to be in microns.

        Returns
        -------
        float or numpy.ndarray
            Geometric cross-section of the grain [cm^2]
        """
        raise NotImplementedError("Subclasses must implement this method.")


from .sphere import Sphere
