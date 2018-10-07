from __future__ import print_function,division
import numpy as np
from numpy import log10
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
    #print('k_h range cut to %.4f - %.4f / (h/Mpc)'%(kmin,kmax))
    return (newk_h,newp_k)

def zero_pad(k_h, p_k, zerokmin=1e-5, zerokmax=1e5, effkmin=1e-3, effkmax=1e2, linear=0, linzero=0):
	# find mean log-step in k_h
	lk_edge, hk_edge = k_h[:3], k_h[-3:]
	l_kpk_edge, h_kpk_edge = lk_edge*p_k[:,:3], hk_edge*p_k[:,-3:]
	logk = log10(k_h)
	mean_diff = np.mean(np.diff(logk))
	lk_zeros = hk_zeros = []

	# extrapolate k-range
	if k_h.min() > effkmin: # if k doesnt reach desired lower limit
		ex_logk_lower = np.arange(start=log10(zerokmin), stop=log10(effkmin), step=mean_diff)[:-1]
		ex_efflogk_lower = np.arange(start=log10(effkmin), stop=log10(k_h.min()), step=mean_diff)[:-1]
		if linear & (not linzero):
			logk = np.array(list(ex_efflogk_lower) + list(logk))
		else:
			logk = np.array(list(ex_logk_lower) + list(ex_efflogk_lower) + list(logk))
		lk_zeros = np.zeros_like(ex_logk_lower)

	elif linear & (not linzero):
		pass
	else:
		ex_logk_lower = np.arange(start=log10(zerokmin), stop=log10(k_h.min()), step=mean_diff)[:-1]
		logk = np.array(list(ex_logk_lower) + list(logk))
		lk_zeros = np.zeros_like(ex_logk_lower)

	if k_h.max() < effkmax: # if k doesnt reach desired upper limit
		ex_logk_upper = np.arange(start=log10(k_h.max()),stop=log10(effkmax),step=mean_diff)[1:]
		ex_efflogk_upper = np.arange(start=log10(effkmax),stop=log10(zerokmax),step=mean_diff)[1:]
		if linear & (not linzero):
			logk = np.array(list(logk) + list(ex_efflogk_upper))
		else:
			logk = np.array(list(logk) + list(ex_efflogk_upper) + list(ex_logk_upper))
		hk_zeros = np.zeros_like(ex_logk_upper)
	elif linear & (not linzero):
		pass
	else:
		ex_logk_upper = np.arange(start=log10(k_h.max()),stop=log10(zerokmax),step=mean_diff)[1:]
		logk = np.array(list(logk) + list(ex_logk_upper))
		hk_zeros = np.zeros_like(ex_logk_upper)

	newk = 10**logk

	if linzero:
		# extend in a power law to effk limits
		# pad with zeros thereafter to zerok limits
		# gradients denoted 'm'
		newpk = np.empty([p_k.shape[0], len(newk)])
		for i in range(len(p_k)):
			pk = -p_k[i]
			if k_h.min() < effkmin:
				lower_interp = lk_zeros
			else:
				log_m_lower = log10( pk[0]/pk[20] ) / log10( k_h[0]/k_h[20] ) # m = dy/dx
				log_c_lower = log10( pk[0] ) - log_m_lower * log10( k_h[0] ) # c = y - mx
				log_interp_lower = log_m_lower * ex_efflogk_lower + log_c_lower # y = mx + c
				lower_interp = np.concatenate((lk_zeros, 10**log_interp_lower))
				lower_interp = -lower_interp

			if k_h.max() > effkmin:
				upper_interp = hk_zeros
			else:
				log_m_upper = log10( pk[-1]/pk[-21] ) / log10( k_h[-1]/k_h[-21] )
				log_c_upper = log10( pk[-1] ) - log_m_upper * log10( k_h[-1] )
				log_interp_upper = log_m_upper * ex_efflogk_upper + log_c_upper
				upper_interp = np.concatenate((10**log_interp_upper, hk_zeros))
				upper_interp = -upper_interp

			newpk[i] = np.concatenate((lower_interp, p_k[i], upper_interp))

	elif not linear:
		# pad p_k with zeros
		my_kmin, my_kmax = np.where(newk < effkmin)[0][-1], np.where(newk > effkmax)[0][0]
		newpk = np.empty([p_k.shape[0], len(lk_zeros)+p_k.shape[1]+len(hk_zeros)])
		for i in range(len(p_k)):
			newpk[i] = np.array(lk_zeros+list(p_k[i])+hk_zeros)
			newpk[i][:my_kmin] = np.zeros_like(newpk[i][:my_kmin])
			newpk[i][my_kmax:] = np.zeros_like(newpk[i][my_kmax:])
	else:
		# extend in a power law to effk limits
		# without zero-padding
		newpk = np.empty([p_k.shape[0], len(newk)])
		for i in range(len(p_k)):
			pk = -p_k[i]
			if k_h.min() < effkmin:
				newpk[i] = p_k[i]
			else:
				log_m_lower = log10( pk[0]/pk[20] ) / log10( k_h[0]/k_h[20] ) # m = dy/dx
				log_c_lower = log10( pk[0] ) - log_m_lower * log10( k_h[0] ) # c = y - mx
				log_interp_lower = log_m_lower * ex_efflogk_lower + log_c_lower # y = mx + c
				lower_interp = 10**log_interp_lower
				lower_interp = -lower_interp
				newpk[i] = np.concatenate((lower_interp, p_k[i]))

			if k_h.max() > effkmin:
				pass
			else:
				log_m_upper = log10( pk[-1]/pk[-21] ) / log10( k_h[-1]/k_h[-21] )
				log_c_upper = log10( pk[-1] ) - log_m_upper * log10( k_h[-1] )
				log_interp_upper = log_m_upper * ex_efflogk_upper + log_c_upper
				upper_interp = 10**log_interp_upper
				upper_interp = -upper_interp
				newpk[i] = np.concatenate((newpk[i], upper_interp))

	return newk, newpk

def interpolate(ups_ks, hod_ks, hod_pk):
	# given new densely sampled k, interpolate P(k) & return k, newP(k)
	pk_cubic_interp = interp1d(hod_ks, hod_pk, axis=1)#, kind=3)
	pk_cubic = pk_cubic_interp(ups_ks)

	return ups_ks, pk_cubic

def loop(k_h, p_k):
	# extrapolate into small-k to match final amplitude of p_k at large-k
	# i.e. 'close' the loop - might help with ringing
	dk = np.diff(log10(k_h)).mean()
	for i in range(len(p_k)):
		pk = -p_k[i]
		kh = k_h.copy()
		m = log10(pk[5] / pk[0]) / log10(kh[5] / kh[0])
		c = pk[5] - m * log10(kh[5])
		#k_extr = log10(kh[0]) - dk
		#pk_extr = m * k_extr + c
		#pk = np.concatenate((m*(log10(kh[0])-dk)+c, pk))
		#kh = np.concatenate((10**(log10(kh[0])-dk), kh))
		#import pdb ; pdb.set_trace()
		while pk[0] > pk[-1]:
			pk = np.append(pk[::-1], m*(log10(kh[0])-dk)+c)[::-1]
			kh = np.append(kh[::-1], 10**(log10(kh[0])-dk))[::-1]
		kf = np.logspace(log10(kh.min()), log10(kh.max()), len(k_h))
		kh, pk = interpolate(kf, kh, pk)
		p_k[i] = -pk

	return k_h, p_k


















