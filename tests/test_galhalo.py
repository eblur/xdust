import pytest
import numpy as np
from scipy.integrate import trapezoid as trapz
import astropy.units as u
from astropy.io import fits
from astropy.table import Table

from xdust.halos import *
from xdust import grainpop
from . import percent_diff

# I tested a few values of NTH to reach the 5% integration threshold in some tests
NE, NTH = 10, 80
EVALS   = np.logspace(-1, 1, NE)   # keV
THVALS  = np.logspace(-1, 4, NTH)  # arcsec

GPOP    = grainpop.make_MRN_RGDrude()

E0, A0, TH0 = 1.0, 0.3, 10.0  # keV, um, arcsec
FABS = 1.0 * np.power(EVALS, -2.0) * np.exp(-0.1 * np.power(EVALS, -2.5))

# Unit Default Input with Astropy Unit
AMIN, AMAX, P, RHO = 0.005*u.micron, 0.5*u.micron, 3.5, 3*u.Unit('g cm^-3') # micron, micron, unitless, g cm^-3
MD = 1.e-6*u.Unit('g cm^-2') # g cm^-2

# Unitless Default Input for testing unit functionality
AMIN_UL, AMAX_UL, P_UL, RHO_UL = 0.005, 0.5, 3.5, 3 # micron, micron, unitless, g cm^-3
MD_UL = 1.e-6 # g cm^-2

# Default Input with Different Unit
AMIN_DU, AMAX_DU, P_DU, RHO_DU = AMIN.to(u.cm), AMAX.to(u.cm), P, RHO.to(u.Unit('kg m^-3')) # micron, micron, unitless, g cm^-3
MD_DU = MD.to(u.Unit('kg m^-2')) # g cm^-2

def test_Halo_dimensions():
    test = Halo(EVALS, THVALS)
    assert len(test.lam) == NE
    assert test.lam.unit == 'keV'
    assert len(test.theta) == NTH

# ---- Test Galactic Halo stuff ---- #

# Tests with default values
UNI_HALO = galhalo.UniformGalHalo(EVALS, THVALS)
SCR_HALO = galhalo.ScreenGalHalo(EVALS, THVALS)
UNI_HALO_CP15 = galhalo.UniformGalHaloCP15(EVALS, THVALS)
SCR_HALO_CP15 = galhalo.ScreenGalHaloCP15(EVALS,THVALS)
# Test for different unit
UNI_HALO_CP15_UL = galhalo.UniformGalHaloCP15(EVALS,THVALS)
SCR_HALO_CP15_UL = galhalo.ScreenGalHaloCP15(EVALS,THVALS)
UNI_HALO_CP15_DU = galhalo.UniformGalHaloCP15(EVALS,THVALS)
SCR_HALO_CP15_DU = galhalo.ScreenGalHaloCP15(EVALS,THVALS)

# Test that calculations run
def test_galhalo_uniform():
    UNI_HALO.calculate(GPOP)
    assert UNI_HALO.norm_int.unit == 'arcsec^-2'
    # Test CP15 version
    UNI_HALO_CP15.calculate(MD)
    assert UNI_HALO_CP15.norm_int.unit == 'arcsec^-2'
    UNI_HALO_CP15_UL.calculate(MD_UL,amin=AMIN_UL, amax=AMAX_UL, p=P_UL, rho = RHO_UL)
    assert UNI_HALO_CP15_UL.norm_int.unit == 'arcsec^-2'
    UNI_HALO_CP15_DU.calculate(MD_DU,amin=AMIN_DU, amax=AMAX_DU, p=P_DU, rho = RHO_DU)
    assert UNI_HALO_CP15_DU.norm_int.unit == 'arcsec^-2'
    # Check consistency of the results from different unit choices
    assert all(percent_diff(np.concatenate(UNI_HALO_CP15.norm_int.value), np.concatenate(UNI_HALO_CP15_UL.norm_int.value)) <= 0.05)
    assert all(percent_diff(np.concatenate(UNI_HALO_CP15.norm_int.value), np.concatenate(UNI_HALO_CP15_DU.norm_int.value)) <= 0.05)

@pytest.mark.parametrize('x', [1.0, 0.5])
def test_galhalo_screen(x):
    SCR_HALO.calculate(GPOP, x)
    assert SCR_HALO.norm_int.unit == 'arcsec^-2'
    
    # Observed angle should be equal to scattering angle when x = 1,
    # so halo should match differential scattering cross section integrated over dust grain size distributions
    if x == 1.0:
        # have to convert int_diff to arcsec^-2
        test = np.abs(SCR_HALO.norm_int.value - GPOP.int_diff.to('arcsec^-2').value)
        assert np.all(test < 0.01)
    # Test CP15 version
    SCR_HALO_CP15.calculate(MD, x=x)
    assert SCR_HALO_CP15.norm_int.unit == 'arcsec^-2'
    SCR_HALO_CP15_UL.calculate(MD_UL,amin=AMIN_UL, amax=AMAX_UL, p=P_UL, rho = RHO_UL, x=x)
    assert SCR_HALO_CP15_UL.norm_int.unit == 'arcsec^-2'
    SCR_HALO_CP15_DU.calculate(MD_DU,amin=AMIN_DU, amax=AMAX_DU, p=P_DU, rho = RHO_DU, x=x)
    assert SCR_HALO_CP15_DU.norm_int.unit == 'arcsec^-2'
    # Check consistency of the results from different unit choices
    assert all(percent_diff(np.concatenate(SCR_HALO_CP15.norm_int.value), np.concatenate(SCR_HALO_CP15_UL.norm_int.value)) <= 0.05)
    assert all(percent_diff(np.concatenate(SCR_HALO_CP15.norm_int.value), np.concatenate(SCR_HALO_CP15_DU.norm_int.value)) <= 0.05)

@pytest.mark.parametrize('test', [UNI_HALO, SCR_HALO, UNI_HALO_CP15, SCR_HALO_CP15])
def test_halos_general(test):
    # Test basic shape properties of outputs
    assert np.shape(test.norm_int) == (NE, NTH)
    assert np.shape(test.taux) == (NE,)

    # Test that the scattering halos integrate to total scattering optical depth
    alph_grid = np.repeat(THVALS.reshape(1, NTH), NE, axis=0)
    int_halo  = trapz(test.norm_int.value * 2.0 * np.pi * alph_grid, THVALS, axis=1)
    assert all(percent_diff(int_halo, test.taux) <= 0.05)

    #--- Test that the flux and intensity functions work
    # Check that calculation runs
    test.calculate_intensity(FABS)
    assert test.intensity is not None
    assert np.size(test.fabs) == NE
    assert np.size(test.fext) == NE
    assert np.size(test.fhalo) == NE
    assert test.intensity.unit == 'arcsec^-2'
    # Integrated intensity profile should add up to total halo flux
    # Calculation uses optically thin approximation, so it adds up to FABS * tau_sca
    tot_halo = trapz(test.intensity * 2.0 * np.pi * THVALS * u.arcsec, THVALS * u.arcsec)
    assert percent_diff(tot_halo, np.sum(test.fabs * test.taux)) <= 0.05

    # The fhalo function calculations halo flux without optically thin approximation
    # So fext = fabs - fhalo
    assert all(percent_diff(test.fext, test.fabs - test.fhalo) <= 0.01)

    # Test that distinction between extincted flux and absorbed flux is handled properly
    FEXT = FABS * np.exp(-test.taux)
    assert all(percent_diff(test.fext, FEXT) <= 0.01)
    fh1  = test.fhalo
    pfa1 = test.percent_fabs
    pfe1 = test.percent_fext
    test.calculate_intensity(FEXT, ftype='ext')
    assert all(percent_diff(test.fext, FEXT) <= 0.01)
    assert all(percent_diff(fh1, test.fhalo) <= 0.01)
    assert percent_diff(pfa1, test.percent_fabs) <= 0.01
    assert percent_diff(pfe1, test.percent_fext) <= 0.01

@pytest.mark.parametrize('test', [UNI_HALO, SCR_HALO, UNI_HALO_CP15, SCR_HALO_CP15])
def test_halo_slices(test):
    # Test the halo slice functions
    # Split the energy values in half
    EMID = (EVALS[-1] - EVALS[0]) / 2.0
    i1 = (EVALS < EMID)
    i2 = (EVALS >= EMID)
    N1, N2 = np.size(EVALS[i1]), np.size(EVALS[i2])
    h1 = test[0:EMID]
    h2 = test[EMID:]
    # Test that the slices return Halo objects
    assert isinstance(h1, Halo)
    assert isinstance(h2, Halo)
    # Test the sizes of things
    assert np.size(h1.lam) == N1
    assert np.size(h2.lam) == N2
    assert np.shape(h1.norm_int) == (N1, NTH)
    assert np.shape(h2.norm_int) == (N2, NTH)
    # Test the values of the halo slices
    assert np.all(percent_diff(np.append(h1.taux,h2.taux), test.taux) == 0.0)
    assert np.all(percent_diff(np.append(h1.lam,h2.lam), test.lam) == 0.0)

@pytest.mark.parametrize('test', [UNI_HALO, SCR_HALO, UNI_HALO_CP15, SCR_HALO_CP15])
def test_halo_index(test):
    i = 0
    tt = test[0]
    assert isinstance(tt, Halo)
    assert np.size(tt.lam) == 1
    assert tt.lam == test.lam[i]
    assert tt.taux == test.taux[i]
    assert np.shape(tt.norm_int) == (1, NTH)
    assert np.shape(tt.intensity) == (NTH,)

@pytest.mark.parametrize('test', [UNI_HALO, SCR_HALO, UNI_HALO_CP15, SCR_HALO_CP15])
def test_halo_io(test):
    test.write('test_halo_io.fits')
    new_halo = Halo(from_file='test_halo_io.fits')
    assert np.all(percent_diff(new_halo.lam,test.lam) < 0.01)
    assert np.all(percent_diff(new_halo.theta,test.theta) < 0.01)
    assert np.all(percent_diff(new_halo.taux,test.taux) < 0.01)
    assert np.all(percent_diff(new_halo.norm_int.flatten(),test.norm_int.flatten()) < 0.01)
    assert new_halo.lam.unit == 'keV'

def create_arf(tmp_path, outfile, evals, effarea) -> str:
    """
    Create a dummy ARF file for testing
    """
    t = Table()
    t['ENERG_LO'] = evals[:-1]
    t['ENERG_HI'] = evals[1:]
    t['SPECRESP'] = np.full(len(evals) - 1, effarea)

    hdu_table = fits.BinTableHDU(t.as_array(), name='SPECRESP')
    hdu_primary = fits.PrimaryHDU()
    hdul = fits.HDUList([hdu_primary, hdu_table])
    hdul.writeto(tmp_path / outfile, overwrite=True)
    return str(tmp_path / outfile)

@pytest.fixture
def halo():
    new_halo = galhalo.UniformGalHalo(EVALS, THVALS)
    new_halo.calculate(GPOP)
    return new_halo

@pytest.fixture
def screen_halo():
    new_halo = galhalo.ScreenGalHalo(EVALS, THVALS)
    new_halo.calculate(GPOP, x=0.5)
    new_halo.calculate_intensity(FABS)
    return new_halo

@pytest.fixture
def arf(tmp_path):
    arf = create_arf(tmp_path, 'test_arf.fits', EVALS, 100.0)
    return arf

def test_fake_image(halo, arf):
    image = halo.fake_image(arf, src_flux=FABS, exposure=1e4, pix_scale=1.0, num_pix=[8, 16])
    assert image.shape == (16, 8)
    assert np.all(image >= 0.0)
    assert np.all(image == np.floor(image))
    assert image.sum() > 0

def test_fake_image_lmin_lmax(halo, arf):
    restricted = halo.fake_image(arf, src_flux=FABS, exposure=1e4, pix_scale=1.0, num_pix=[8, 16], lmin=0.5, lmax=5.0)
    unrestricted = halo.fake_image(arf, src_flux=FABS, exposure=1e4, pix_scale=1.0, num_pix=[8, 16])
    assert 0 < restricted.sum() < unrestricted.sum()

def test_time_delay():
    result1 = galhalo.time_delay(100.0, 0.5, 8.0)
    assert result1 == pytest.approx(96770, rel=1e-3)

    result2 = galhalo.time_delay(200.0, 0.5, 8.0)
    assert result2 == pytest.approx(4 * result1, rel=1e-3)

def test_variable_profile(screen_halo):
    time = np.linspace(0, 100, 10)
    lc = np.ones_like(time)
    dist = 8.0
    var_profile1 = screen_halo.variable_profile(time, lc, dist=dist)
    assert var_profile1.shape == (len(EVALS), len(THVALS))
    assert np.sum(var_profile1) > 0

    var_profile2 = screen_halo.variable_profile(time, 2 * lc, dist=dist)
    assert var_profile2.shape == (len(EVALS), len(THVALS))
    assert np.sum(var_profile2) > 0
    # 2 * lc and assert the result is ~2×

def test_fake_variable_image(screen_halo, arf):
    time = np.linspace(0, 100, 10)
    lc = np.ones_like(time)
    image = screen_halo.fake_variable_image(time, lc=lc, arf=arf, exposure=1e4, num_pix=[8, 16])
    assert image.shape == (16, 8)
    assert np.all(image >= 0.0)
    assert np.all(image == np.floor(image))
    assert image.sum() > 0
