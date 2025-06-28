
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from astropy.io import fits
from astropy.wcs import WCS

def between(a, x, b):
    return (a <= x) & (x < b)


def cosd(theta):
    """Return cosine of theta in degrees"""
    return np.cos(np.deg2rad(theta))


def sind(theta):
    """Return sine of theta in degrees"""
    return np.sin(np.deg2rad(theta))


def tand(x):
    """Return tangent of theta in degrees"""
    return sind(x) / cosd(x)


def xieta2radec(xi, eta, A, D):
    """Tangent-plane projection around A, D (with all params in degrees) to ra, dec"""
    DEGREE = np.pi / 180.
    scalar = np.isscalar(xi)
    if (not scalar):
        assert (len(xi) == len(eta))

    xirad = np.atleast_1d(xi*DEGREE)
    etarad = np.atleast_1d(eta*DEGREE)

    rho = np.hypot(xirad, etarad)
    cang = np.arctan(rho)

    alpha = xirad*0.
    delta = xirad*0.

    rholimit = 1.e-8
    set = np.where(rho >= rholimit)
    if (set[0].size > 0):
        alpha1 = A + (1./DEGREE) * np.arctan((xirad[set]*np.sin(cang[set]))
           / (rho[set]*cosd(D)*np.cos(cang[set]) - etarad[set]*sind(D)*np.sin(cang[set])))
        delta1 = (1./DEGREE) * np.arcsin(np.cos(cang[set])*sind(D)
                          + etarad[set]*np.sin(cang[set])*cosd(D)/rho[set])
        alpha[set] = alpha1
        delta[set] = delta1
    set = np.where(rho < rholimit)
    if (set[0].size > 0):
        alpha1 = A + (1./DEGREE) * np.arctan(xirad[set]
                     / (cosd(D) - etarad[set]*sind(D)))
        delta1 = (1./DEGREE) * np.arcsin(sind(D) + etarad[set]*cosd(D))
        alpha[set] = alpha1
        delta[set] = delta1

    if (scalar):
        alpha = np.asscalar(alpha)
        delta = np.asscalar(delta)

    return alpha, delta


def is_sorted(a):
    for i in range(a.size-1):
        if a[i+1] < a[i] :
           return False
    return True


class FitsImage(object):
    """Utility class for handling FITS images"""
    def __init__(self, filename, unit=0):
        from IPython.utils import io
        with io.capture_output() as captured:
            with fits.open(filename) as hdul:
                hdul.info()
                self.header = hdul[unit].header
                self.data = hdul[unit].data
        self.wcs = WCS(self.header, naxis=2)  # without naxis=2, barfs on color images
        if len(self.data.shape) == 3:
            self.data = np.moveaxis(self.data, 0, -1)

            
def getmasktable(filename):
    """Read mask table fits file and return the most useful fields.
    Returns astropy table for slit/target properties, and dict for mask-level information.
    Based on very limited experience, this routine works, but check it carefully!"""
    hdul = fits.open(filename)
    assert is_sorted(hdul[1].data['ObjectId'])
    assert is_sorted(hdul[4].data['dSlitId'])
    assert (hdul[5].data['dSlitId'] == hdul[4].data['dSlitId']).all()
    assert (hdul[7].data['dSlitId'] == hdul[4].data['dSlitId']).all()
    ind = np.searchsorted(hdul[1].data['ObjectId'], hdul[5].data['ObjectId'])
    # slit table
    t = Table()
    for key in ['OBJECT', 'RA_OBJ', 'DEC_OBJ']:
        t[key] = hdul[1].data[key][ind]
    for key in ['dSlitId', 'SlitName', 'slitRA', 'slitDec', 'slitTyp', 'slitLen', 'slitLPA', 'slitWid']:
        t[key] = hdul[4].data[key]
    for key in ['TopDist', 'BotDist']:
        t[key] = hdul[5].data[key]
    for key in ['slitX1', 'slitY1', 'slitX2', 'slitY2', 'slitX3', 'slitY3', 'slitX4', 'slitY4']:
        t[key] = hdul[7].data[key]
    # general mask info
    info = dict()
    for key in ['DesName', 'DesAuth', 'DesDate', 'ProjName', 'INSTRUME', 'RA_PNT', 'DEC_PNT', 'PA_PNT']:
        info[key] = hdul[3].data[key][0]
    return t, info