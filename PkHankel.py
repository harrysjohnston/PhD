from __future__ import print_function,division
import numpy as np
from astropy.io import fits,ascii
from scipy.interpolate import PchipInterpolator as interp1d

def read_z(shapes,dens):
    shap_z = np.loadtxt(shapes)
    dens_z = np.loadtxt(dens)
    return shap_z,dens_z

def create_nz(gal_z, nz, nbin, z_range):
    # nz gives n(z) resolution, and nbin gives number of tomographic bins - required for Cl's
    # for my IA analysis nbin doesn't actually do anything, since my weight function only wants
    # a 1d array; just set nz fine enough for nice redshift interpolation of Pks
    zmin, zmax = z_range.min(), z_range.max()
    z_hist = np.histogram(gal_z, bins=nz, range=(zmin, zmax))
    n_of_z = z_hist[0]
    x = z_hist[1]
    z_mids = x[1:]-(x[1]-x[0])/2
    nz_table = np.empty([nz,nbin+1])
    nz_table[:,0] = z_mids
    bin_ranges = np.linspace(zmin, zmax, num=nbin+1)
    for i in range(nbin):
        nz_table[:,i+1] = np.where(
            (bin_ranges[i]<=z_mids)&(z_mids<bin_ranges[i+1]),n_of_z,0 )
    return z_mids, n_of_z

def cut_krange(k_h,p_k,kmin=10**-2.2,kmax=10**1.2):
    # limit range in k, prior to FFT (Hankel transfm), to avoid ringing
    kcut = (k_h>=kmin)&(k_h<=kmax)
    newk_h = k_h[kcut]
    #newk_h = k_h
    newp_k = np.empty([len(p_k),len(newk_h)])
    for i in range(len(p_k)):
        newp_k[i] = p_k[i][kcut]
    print('k_h range cut to %.4f - %.4f / (h/Mpc)'%(kmin,kmax))
    return (newk_h,newp_k)

def zero_pad(k_h, p_k, zerokmin=1e-5, zerokmax=1e5, effkmin=1e-3, effkmax=1e2, linear=0):
    # find mean log-step in k_h
    lk_edge, hk_edge = k_h[:3], k_h[-3:]
    l_kpk_edge, h_kpk_edge = lk_edge*p_k[:,:3], hk_edge*p_k[:,-3:]
    logk = np.log10(k_h)
    mean_diff = np.mean(np.diff(logk))
    lk_zeros = hk_zeros = []

    # extrapolate k-range
    if k_h.min()>zerokmin:
        ex_logk = np.arange(start=np.log10(zerokmin),stop=np.log10(k_h.min()),step=mean_diff)[:-1]
	l_ex_logk = ex_logk.copy()
        logk = np.array(list(ex_logk)+list(logk))
        lk_zeros = [0.]*len(ex_logk)
    if k_h.max()<zerokmax:
        ex_logk = np.arange(start=np.log10(k_h.max()),stop=np.log10(zerokmax),step=mean_diff)[1:]
	h_ex_logk = ex_logk.copy()
        logk = np.array(list(logk)+list(ex_logk))
        hk_zeros = [0.]*len(ex_logk)
    newk = 10**logk
    print('PADDING OUTSIDE SCALES OF INTEREST: ')

    # pad p_k with zeros
    if not linear:
    	my_kmin, my_kmax = np.where(newk<effkmin)[0][-1], np.where(newk>effkmax)[0][0]
    	newpk = np.empty([p_k.shape[0],len(lk_zeros)+p_k.shape[1]+len(hk_zeros)])
    	for i in range(len(p_k)):
    	    newpk[i] = np.array(lk_zeros+list(p_k[i])+hk_zeros)
    	    newpk[i][:my_kmin] = np.zeros_like(newpk[i][:my_kmin])
    	    newpk[i][my_kmax:] = np.zeros_like(newpk[i][my_kmax:])

    else:
    	newpk = np.empty([p_k.shape[0],len(l_ex_logk)+p_k.shape[1]+len(h_ex_logk)])
	for i in range(len(p_k)):
		logm = np.log10(l_kpk_edge[i,-1]/l_kpk_edge[i,0]) / np.log10(lk_edge[-1]/lk_edge[0])
		logc = np.log10(l_kpk_edge[i][0]) - logm*np.log10(lk_edge[0])
		log_interp = logm*(l_ex_logk) + logc
		lk_interp = log_interp.copy()
		logm = np.log10(h_kpk_edge[i,-1]/h_kpk_edge[i,0]) / np.log10(hk_edge[-1]/hk_edge[0])
                logc = np.log10(h_kpk_edge[i][0]) - logm*np.log10(hk_edge[0])
                log_interp = logm*(h_ex_logk) + logc
		hk_interp = log_interp.copy()
		newpk[i] = np.array(list(10**lk_interp) + list(k_h*p_k[i]) + list(10**hk_interp))
		newpk[i] = newpk[i]/newk

    return newk, newpk

def interpolate(ups_ks, hod_ks, hod_pk):
	# given new densely sampled k, interpolate P(k) & return k, newP(k)
	pk_cubic_interp = interp1d(hod_ks, hod_pk, axis=1)#, kind=3)
	pk_cubic = pk_cubic_interp(ups_ks)

	return ups_ks, pk_cubic





