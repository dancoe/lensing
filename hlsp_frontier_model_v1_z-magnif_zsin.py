# http://archive.stsci.edu/pub/hlsp/frontier/abell2744/models/
# hlsp_frontier_model_v1_z-magnif.py

#################################
# Derive lens model magnification map for a given source redshift
#   from maps of mass (kappa) and shear (gamma)
#
# Dan Coe <DCoe@STScI.edu>
#
# python magnif.py inkappa.fits inshear.fits lensredshift inredshift outmagnif.fits outredshift

import os, sys
from numpy import sqrt
import astropy.io.fits as pyfits
from scipy.integrate import quad

#################################
# Cosmology

h = 0.7  # Hubble constant / (100 km/s/Mpc)
Om = 0.3  # Omega matter
OL = 1 - Om  # Omega Lambda

H100 = 3.24077674937e-18  # 100 km/s/Mpc in units of 1 / S (inverse seconds)
c = 9.7155e-15  # MPC / S

def H(z):
    """Hubble constant [1/S] at redshift z for a flat universe"""
    return h * H100 * sqrt(Om * (1+z)**3 + OL)

def Hinv(z):
    return 1 / H(z)

def DA(z1, z2):
    """Angular-diameter distance (MPC)
    between two redshifts for a flat universe"""
    #return c / (1.+z2) * integral(lambda z: 1/H(z), z1, z2)
    return c / (1.+z2) * quad(Hinv, z1, z2)[0]

def Dds_Ds(zl, zs):
    """Ratio of angular diameter distances;
    Lensing deflections scale with this quantity"""
    if zs == float('Inf'):
        return 1
    else:
        Dds = DA(zl, zs)
        Ds  = DA(0,  zs)
        return Dds / Ds

#################################
# FITS image handling

def savefits(hdu, data, outfile, overwrite=True):
    if os.path.exists(outfile):
        if overwrite:
            os.remove(outfile)
    hdu[0].data = data
    print 'SAVING %s...' % outfile
    hdu.writeto(outfile, output_verify='fix')

#################################

def savemagnif(inkappa, inshear, zl, zsin, outfile, zsout):
    Dds_Dsin  = Dds_Ds(zl, zsin)   # Factor out old zs
    Dds_Dsout = Dds_Ds(zl, zsout)  # Factor in new zs

    hdu = pyfits.open(inkappa)
    kappa = hdu[0].data

    gamma = pyfits.open(inshear)[0].data

    kappa = kappa * Dds_Dsout / Dds_Dsin
    gamma = gamma * Dds_Dsout / Dds_Dsin

    magnif = 1 / ((1 - kappa)**2 - gamma**2)
    magnif = abs(magnif)

    savefits(hdu, magnif, outfile)


if __name__ == '__main__':
    if not len(sys.argv) > 1:
        print
        print 'USAGE: python hlsp_frontier_model_v1_z-magnif.py inkappa.fits inshear.fits lensredshift inredshift outmagnif.fits outredshift'
        print 'EXAMPLE: python hlsp_frontier_model_v1_z-magnif.py hlsp_frontier_model_abell2744_cats_v1_kappa.fits hlsp_frontier_model_abell2744_cats_v1_gamma.fits 0.308 2.55 magnif_z11.fits 11'
        print 'Requires scipy and pyfits'
        print
        quit()

    inkappa = sys.argv[1]
    inshear = sys.argv[2]
    zl = float(sys.argv[3])     # Lens redshift
    zsin = float(sys.argv[4])   # Source redshift images are calibrated to (input)
    outfile = sys.argv[5]
    zsout = float(sys.argv[6])  # Source redshift you are interested in (output)
    savemagnif(inkappa, inshear, zl, zsin, outfile, zsout)
