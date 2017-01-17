from __future__ import print_function,division
import numpy as np
from astropy.io import fits
from astropy.io import ascii

def read_z(path):
    fitsdata = fits.open(path)
    data = fitsdata[1].data
    if 'DEI' in path:
        zcol = 'Z'
    else:
        zcol = 'Z_1_1'
    gal_z = data[zcol]
    return gal_z

def create_nz(gal_z,nz,nbin,outfile):
    # NOTE: this is written specifically for zrange=[0,0.5]
    z_hist = np.histogram(gal_z,bins=nz,range=(0.,0.5))
    n_of_z = z_hist[0]
    x = z_hist[1]
    z_mids = x[1:]-(x[1]-x[0])/2
    nz_table = np.empty([nz,nbin+1])
    nz_table[:,0] = z_mids
    bin_ranges = np.linspace(0,0.5,num=nbin+1)
    for i in range(nbin):
        nz_table[:,i+1] = np.where(
            (bin_ranges[i]<=z_mids)&(z_mids<bin_ranges[i+1]),n_of_z,0)
    # outfile = outfile+'_nz%snbin%s.asc'%(nz,nbin)
    ascii.write(nz_table,outfile,names=['#z_mid']+['#bin_%s'%(i+1) for i in range(nbin)],delimiter='\t')
    return n_of_z


