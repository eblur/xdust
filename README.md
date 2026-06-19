# Xdust

(
    *Previously eblur/newdust* [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.7500048.svg)](https://doi.org/10.5281/zenodo.7500048) 
)

This package calculates extinction curves and small-angle scattering halos from a user-defined dust grain size distribution. It calculates scattering and absorption from first principles (optical constants of the material, Mie or Rayleigh-Gans scattering).

**Updated 2026.06** to by installable via PyPI and documented with Sphinx

**Documentation:** https://eblur.github.io/xdust

## Installation

### Typical installation

```
pip install xdust
```

### For developers

If you need to stay up to date with development versions of xdust, use:

```
git clone https://github.com/eblur/xdust.git

cd xdust

pip install -e .
```

## To invoke:

```
import xdust
```

## How to use:

See the jupyter notebooks in **examples/** 
for examples of setting up grain populations and modeling scattering halos from Galactic dust.

### For simulating cosmological halos

(e.g. [Corrales & Paerels, 2012](http://adsabs.harvard.edu/abs/2012ApJ...751...93C) and 
[Corrales 2015](http://adsabs.harvard.edu/abs/2015ApJ...805...23C))
see the `cosmhalo` branch. Some of the cosmhalo tests do not pass. Use with caution.

