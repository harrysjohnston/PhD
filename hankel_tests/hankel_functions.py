from __future__ import print_function, division
import numpy as np
from numpy import pi, log10, exp
from scipy.interpolate import interp1d, interp2d
from scipy.integrate import quad
import hankel
from hankel import get_h
from astropy.cosmology import FlatLambdaCDM as FLCDM
cosmo = FLCDM(Om0=0.25, H0=70, Ob0=0.044)

def prep_kpar_integrand(kpar, kperp, Pi):
	karg = np.sqrt(np.add.outer(kperp**2, kpar**2))
	kparm, Pim = np.meshgrid(kpar, Pi)
	cosine = np.cos(kparm * Pim)
	return karg, cosine

def make_kpar_integrand_2d(k, pk, karg, cosine):
	P = interp1d(k, pk, fill_value='extrapolate', bounds_error=0)
	return P(karg)[:, :, np.newaxis] * cosine.T

def make_kpar_integrand_1d(k, pk, kperp):#, pi):
	P = interp1d(k, pk, fill_value='extrapolate', bounds_error=0)
	def make_func(kpar):
		kmode = (kperp**2 + kpar**2)**0.5
		return P(kmode)# * np.cos(kpar*pi)
	return make_func

def k_scale_cuts(kx, low, high):
	return kx[(kx >= low) & (kx <= high)]

def upsample_k(kx, factor, spacing='log', n=None):
	if n is None:
		n = int(factor * len(kx))
	if spacing == 'log':
		return np.logspace(np.log10(kx.min()), np.log10(kx.max()), n)
	if spacing == 'lin':
		return np.linspace(kx.min(), kx.max(), n)

def fivept_stencil(func, x, h):
	# returns f'(x), via 5pt stencil, for grid-spacing h
	return (-func(x+2*h)+8*func(x+h)-8*func(x-h)+func(x-2*h))/(12*h)

def zero_pad(k, pk, nzeros=10):
	kl = 10**(np.log10(k[0]) - np.arange(1, nzeros+1) * np.log10(k[1]/k[0]))[::-1]
	ku = 10**(np.log10(k[-1]) + np.arange(1, nzeros+1) * np.log10(k[-1]/k[-2]))
	p0 = np.zeros_like(ku)
	k = np.concatenate((kl, k, ku))
	if pk.ndim == 2:
		pknew = []
		for i in range(len(pk)):
			pknew.append(np.concatenate((p0, pk[i], p0)))
		pknew = np.array(pknew)
	else:
		pknew = np.concatenate((p0, pk, p0))
	return k, pknew

def choose_h(f, nu, r_min, r_max):
	h_opt, f_at_r, N_opt = get_h(f, nu, np.array([r_min, r_max]))
	if h_opt < 1e-3:
		h_opt = 1e-3
	return h_opt

def choose_kpar_max(Pi):
	#cosine_zero = [pi * (1. / abs(2.*p)) for p in Pi]
	cosine_zero = [np.inf for p in Pi]
	return np.array(cosine_zero)

def compute_Wz(z, nz1, nz2):
	# Wz = [p^2 / X^2*X'] / int[p^2 / X^2*X' dz]
	# compute p(z) = unconditional pdf
	pz1 = nz1 / nz1.sum()
	pz2 = nz2 / nz2.sum()
	assert pz1.shape == z.shape, "p(z) vs. z mismatch"
	assert pz2.shape == z.shape, "p(z) vs. z mismatch"
	# compute X(z) = comoving coordinate
	chiz = lambda x: cosmo.comoving_distance(x) # h (hubble parameter) cancels
	X = chiz(z).value
	Xsq = X**2
	# compute X'(z) = first deriv.
	h = z[1]-z[0]
	Xpr = fivept_stencil(chiz, z, h).value
	# combine & integrate (Riemann sum) over z
	W_nom = (pz1 * pz2) / (Xsq * Xpr)
	W_nom = np.nan_to_num(W_nom)
	W_dom = np.trapz(W_nom, x=z)
	W = W_nom / W_dom
	return W




