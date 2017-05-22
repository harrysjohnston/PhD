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
    newk_h = k_h[kcut]
    #newk_h = k_h
    newp_k = np.empty([len(p_k),len(newk_h)])
    for i in range(len(p_k)):
        newp_k[i] = p_k[i][kcut]
        #newp_k[i] = np.where(kcut==True,p_k[i],0)
    print('k_h range cut to %.4f - %.4f / (h/Mpc)'%(kmin,kmax))
#    print('approx. equiv. to %.4f - %.4f / (Mpc/h)'%(1/kmax,1/kmin))
    return (newk_h,newp_k)

def zero_pad(k_h,p_k,kmin=1e-5,kmax=1e5):
    # find mean log-step in k_h
    logk = np.log10(k_h)
    mean_diff = np.mean(np.diff(logk))
    lk_zeros = hk_zeros = []

    # extrapolate k-range
    if k_h.min()>kmin:
        ex_logk = np.arange(start=np.log10(kmin),stop=np.log10(k_h.min()),step=mean_diff)[:-1]
        logk = np.array(list(ex_logk)+list(logk))
        lk_zeros = [0.]*len(ex_logk)
    if k_h.max()<kmax:
        ex_logk = np.arange(start=np.log10(k_h.max()),stop=np.log10(kmax),step=mean_diff)[1:]
        logk = np.array(list(logk)+list(ex_logk))
        hk_zeros = [0.]*len(ex_logk)
    newk = 10**logk
    print('PADDING OUTSIDE SCALES OF INTEREST: ')
    my_kmin, my_kmax = np.where(newk<10**-2.2)[0][-1], np.where(newk>10**1.2)[0][0]

    # pad p_k with zeros
    newpk = np.empty([p_k.shape[0],len(lk_zeros)+p_k.shape[1]+len(hk_zeros)])
    for i in range(len(p_k)):
        newpk[i] = np.array(lk_zeros+list(p_k[i])+hk_zeros)
        newpk[i][:my_kmin] = np.zeros_like(newpk[i][:my_kmin])
        newpk[i][my_kmax:] = np.zeros_like(newpk[i][my_kmax:])

    return newk, newpk



