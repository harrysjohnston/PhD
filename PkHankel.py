from __future__ import print_function,division
import numpy as np
from astropy.io import fits,ascii

def read_z(shapes,dens):
    shap_z = np.loadtxt(shapes)
    dens_z = np.loadtxt(dens)
    return shap_z,dens_z

def create_nz(gal_z,nz,nbin,outfile):
    # Note: this is written specifically for zrange=[0,0.5]
    # AND would need editing for:  nz != nbin 
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
    if outfile!=None:
        ascii.write(nz_table,outfile,names=['#z_mid']+['#bin_%s'%(i+1) for i in range(nbin)],delimiter='\t')
    return z_mids,n_of_z

def cut_krange(k_h,p_k,kmin=10**-2.2,kmax=10**1.2):
    # limit range in k, prior to FFT (Hankel transfm), to avoid ringing
    kcut = (k_h>=kmin)&(k_h<=kmax)
    #newk_h = k_h[kcut]
    newk_h = k_h
    newp_k = np.empty([len(p_k),len(newk_h)])
    for i in range(len(p_k)):
        #newp_k[i] = p_k[i][kcut]
        newp_k[i] = np.where(kcut==True,p_k[i],0)
    print('k_h range cut to %.4f - %.4f / (h/Mpc)'%(kmin,kmax))
    print('approx. equiv. to %.4f - %.4f / (Mpc/h)'%(1/kmax,1/kmin))
    return (newk_h,newp_k)
