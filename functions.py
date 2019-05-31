from __future__ import print_function, division
try:
	from IPython import get_ipython
	ipython = get_ipython()
	ipython.magic("matplotlib")
	#ipython.magic("load_ext autoreload")
	#ipython.magic("autoreload 2")
except:
	print(' -- failed to import ipython --')
import os
import pickle
import numpy as np
import healpy as hp
import scipy
import scipy.interpolate as scint
import scipy.stats as stats
import pandas as pd
from numpy import log10, log
from scipy.interpolate import interp1d, interp2d
from scipy.optimize import curve_fit
norm = stats.norm
from astropy.io import fits, ascii
fitscol = fits.Column
from astropy.table import Table, vstack
from astropy import cosmology
from astropy import units as u
from astropy.cosmology import FlatLambdaCDM as FLCDM
MICEcosmo = FLCDM(Om0=0.25, H0=70, Ob0=0.044)
cosmo = FLCDM(Om0=0.3, H0=70)
cMpc = MICEcosmo.comoving_distance
from os import listdir, mkdir
from os.path import join, isdir, basename, normpath, dirname
def phd():
	os.chdir('/share/splinter/hj/PhD/')
	print('/share/splinter/hj/PhD/')
def euc():
	os.chdir('/share/splinter/hj/PhD/cosmosis/modules/euclid_ias/')
	print('/share/splinter/hj/PhD/cosmosis/modules/euclid_ias/')
def mill():
	os.chdir('/share/data1/hj/millennium/')
	print('/share/data1/hj/millennium/')
def kids():
	os.chdir('/share/splinter/hj/PhD/KiDS_PhotometricClustering')
	print('/share/splinter/hj/PhD/KiDS_PhotometricClustering')
def pau():
	os.chdir('/share/splinter/hj/PhD/PAU')
	print('/share/splinter/hj/PhD/PAU')
import matplotlib as mpl
#mpl.use('Agg')
mpl.rcParams['xtick.labelsize'] = 12
mpl.rcParams['ytick.labelsize'] = 12
mpl.rcParams['axes.labelsize'] = 18
mpl.rcParams['xtick.direction'] = 'in'
mpl.rcParams['xtick.top'] = True
mpl.rcParams['xtick.major.size'] = 6
mpl.rcParams['xtick.minor.size'] = 3.0
mpl.rcParams['xtick.minor.visible'] = True
mpl.rcParams['ytick.direction'] = 'in'
mpl.rcParams['ytick.right'] = True
mpl.rcParams['ytick.major.size'] = 6
mpl.rcParams['ytick.minor.size'] = 3.0
mpl.rcParams['ytick.minor.visible'] = True
mpl.rcParams['axes.linewidth'] = 1.5
import matplotlib.pyplot as plt
import matplotlib.lines as ml
from matplotlib.patches import Ellipse
import matplotlib.colors as colors
import matplotlib.gridspec as gridspec
import sys
sys.path.append("/share/splinter/hj/PhD/")
sys.path.append("/share/splinter/hj/PhD/KiDS_PhotometricClustering")
sys.path.append("/share/splinter/hj/PhD/KiDS_PhotometricClustering/Randoms")
sys.path.append("/share/splinter/hj/PhD/cosmosis/modules/euclid_ias/")
sys.path.append("~/PhD/")
sys.path.append("~/PhD/PAU")
sys.path.append("~/PhD/KiDS_PhotometricClustering")
sys.path.append("~/PhD/KiDS_PhotometricClustering/Randoms")
import skyknife
try:
	import MCplotter as mcp
	import MakeRandoms as MR
	import Clustering_plotter as clp
	from smail_pz_func import *
	import kdestats as kde
	import catalog_sampler as CS
	import jackknife3d as jk3d
	import downsampler as ds
	import compare_coords as cc # plot3dhist() and plot3dims() - each takes 3x 2darray catalogs eg. shapes, dens, randoms.
					# expand to take *args?
	ds_jkfunc = ds.make_jks
except:
	pass
down = lambda x, f: x[np.random.rand(len(x)) <= f]
frac = lambda x: np.sum(x)/np.float32(len(x))
minmax = lambda x: (x.min(), x.max())
midpoints = lambda x: (x[1:] + x[:-1]) / 2.
pdf_baseline = lambda l, s: np.linspace( norm.ppf(0.01, loc=l, scale=s), norm.ppf(0.99, loc=l, scale=s), 100 )
Area = lambda a1,a2,d1,d2: (a2-a1)*(np.sin(d2)-np.sin(d1))
hargs = {'bins':'auto', 'histtype':'step', 'normed':1}
hartlap = lambda N, D: (N - D - 2.) / (N - 1.)
dotdot = lambda v1, M, v2: np.dot(np.dot(v1, M), v2)
chi2pval = lambda chi2, dof: 1. - stats.chi2.cdf(chi2, dof)
p2sigma = lambda p: stats.norm.interval(1. - p, loc=0, scale=1)[1]
new_fitscol = lambda arr, nam, form: fits.ColDefs([fits.Column(array=arr, format=form, name=nam)])
fopen = lambda x: fits.open(x)[1].data
ropen = lambda x: open(x, 'r')
wopen = lambda x: open(x, 'w')
popen = lambda x: pickle.load(ropen(x))

def plotfile(path, **kwargs):
	t, w = np.loadtxt(path, **kwargs).T[:2]
	def plot(**kwargs):
		return plt.plot(**kwargs)
	plot(t, w)

def print_tile_cuts(N=50, max_tile=700):
    n_tiles_per_group = int(max_tile / N)
    tile_bins = np.arange(0, max_tile+1, step=n_tiles_per_group)
    string = ""
    for i in range(len(tile_bins)-1):
        t1, t2 = tile_bins[i:i+2]
        string += "\n(tileID >= %s) & (tileID < %s) //" % (t1, t2)
    print(string)
    print('\n\n\t ~%s deg^2 per grouping'%(t2-1-t1))

def path_to_map(filename, nside=1024, cut_col=None, shift_ra=0, colnames=None):
	if colnames is None:
		colnames = ('RA', 'DEC')
	cat = fits.open(filename)[1].data
	if cut_col is not None:
		if hasattr(cut_col, '__iter__'):
			cut = np.ones(len(cat), dtype=bool)
			for cc in cut_col:
				cut &= cat[cc]
		else:
			cut = cat[cut_col]
		cat = cat[cut]
	ra = cat[colnames[0]]
	dec = cat[colnames[1]]
	if shift_ra:
		ra = ra + 180.
		ra = np.where(ra > 360, ra - 360, ra)
	hpix = hp.ang2pix(nside, ra, dec, lonlat=1)
	hmap = np.bincount(hpix, minlength=hp.nside2npix(nside))
	return hmap
def radec_to_map(ra, dec, nside=1024, shift_ra=0):
	if shift_ra:
		ra = ra + 180.
		ra = np.where(ra > 360, ra - 360, ra)
	hpix = hp.ang2pix(nside, ra, dec, lonlat=1)
	hmap = np.bincount(hpix, minlength=hp.nside2npix(nside))
	return hmap

from matplotlib.ticker import Locator
class MinorSymLogLocator(Locator):
	"""
	Dynamically find minor tick positions based on the positions of
	major ticks for a symlog scaling.
	"""
	def __init__(self, linthresh):
		"""
		Ticks will be placed between the major ticks.
		The placement is linear for x between -linthresh and linthresh,
		otherwise its logarithmically
		"""
		self.linthresh = linthresh

	def __call__(self):
		'Return the locations of the ticks'
		majorlocs = self.axis.get_majorticklocs()
		#if majorlocs[1] == 10*majorlocs[0]:
		majorlocs = np.insert(majorlocs, 0, majorlocs[0]/10.)
		#else:
		#	majorlocs = np.insert(majorlocs, 0, majorlocs[0]-np.diff(majorlocs)[0])
		#if majorlocs[-1] == 10*majorlocs[-2]:
		majorlocs = np.insert(majorlocs, -1, majorlocs[-1]*10.)
		#else:
		#	majorlocs = np.insert(majorlocs, -1, majorlocs[-1]+np.diff(majorlocs)[-1])

		# iterate through minor locs
		minorlocs = []

		# handle the lowest part
		for i in xrange(1, len(majorlocs)):
			majorstep = majorlocs[i] - majorlocs[i-1]
			if abs(majorlocs[i-1] + majorstep/2) < self.linthresh:
				ndivs = 10
			else:
				ndivs = 9
			minorstep = majorstep / ndivs
			locs = np.arange(majorlocs[i-1], majorlocs[i], minorstep)[1:]
			minorlocs.extend(locs)

		return self.raise_if_exceeds(np.array(minorlocs))

	def tick_values(self, vmin, vmax):
		raise NotImplementedError('Cannot get tick locations for a '
								  '%s type.' % type(self))

def read_files(ident='', idnot=None, d=1, dire='.', cols=None, quiet=1, asc=0):
	ident = ident.split('*')
	if idnot is not None:
		idnot = idnot.split('*')
	else:
		idnot = []
	ldir = [i for i in listdir(dire) if all([idt in i for idt in ident]) & all([idn not in i for idn in idnot])]
	ldir.sort()
	file_array = []
	for ld in ldir:
		ld = join(dire, ld)
		if not quiet:
			print(ld)
		if asc:
			file_array.append(ascii.read(ld))
		else:
			try:
				file_array.append(np.loadtxt(ld, usecols=cols))
			except ValueError:
				file_array.append(np.loadtxt(ld, skiprows=1, usecols=cols))
	if d:
		file_dict = dict(zip(ldir, file_array))
		return file_dict
	else:
		return np.asarray(file_array)

def make_randoms(_real_sample, randoms=None, Nrand=10, sky_box=1, plot=0):
	"""
	real_sample: real position catalogue - needs redshifts or comoving dists
	randoms: give random catalogue (2+ column array, where first two are RA/DEC) to draw positions from (optional)
	Nrand: target density of output randoms, in multiples of real_sample size, or (if not (randoms is None))
			in multiples of random sample size
	sky_box: randomize positions in 1x1 RA/DEC box, for scaling up/translating
	plot: make plot of 3D distributions of in/output catalogues

	returns: random catalogue (, axes - if plot=1)
	"""
	if _real_sample.ndim == 1:
		real_z = _real_sample.copy()
		dummy_ra = np.random.rand(len(real_z))
		dummy_dec = np.random.rand(len(real_z))
		real_sample = np.column_stack((dummy_ra, dummy_dec, _real_sample.copy()))
	else:
		real_sample = _real_sample.copy()
		real_z = real_sample.copy().T[2]


	out_size = Nrand * len(real_sample)
	if randoms is not None:
		out_size = Nrand * len(randoms)
		rand_radec = randoms.copy()[:, :2]
		RAND_NUM = np.random.choice(range(len(randoms)), size=out_size)
		new_radec = randoms.copy()[:,:2][RAND_NUM]
	elif sky_box:
		ra = np.random.rand(out_size)
		dec = np.random.rand(out_size)
		new_radec = np.column_stack((ra, dec))

	rand_z = fit_smail(real_z)
	new_rand_z = np.random.choice(rand_z.T[0], p=rand_z.T[1], size=out_size)
	new_rcat = np.column_stack((new_radec, new_rand_z))

	if plot:
		if randoms is not None:
			ax = cc.plot3dhist(real=real_sample, rand=randoms)
		else:
			ax = cc.plot3dhist(real=real_sample)
		h0 = ax[0].hist(new_rcat.T[0], normed=1, histtype='step', color='r', lw=1.5, bins=100)
		h1 = ax[1].hist(new_rcat.T[1], normed=1, histtype='step', color='r', lw=1.5, bins=100)
		h2 = ax[2].hist(new_rcat.T[2], normed=1, histtype='step', color='r', lw=1.5, bins=100)
		plt.show()
		return new_rcat, ax
	else:
		return new_rcat

def myhist(variable, **kwargs):
	if 'bins' not in kwargs.keys():
		kwargs['bins'] = 'auto'
	if 'histtype' not in kwargs.keys():
		kwargs['histtype'] = 'step'
	if 'normed' not in kwargs.keys():
		kwargs['normed'] = 1
	return plt.hist(variable, **kwargs)

def corr_mat(cov, title=None, cmap='RdYlGn', zero_norm=1):
	if type(cov)==str:
		cov = tryload(cov)
	R = pearson_r(cov.copy())
	f, ax = plt.subplots()
	if zero_norm:
		c = ax.imshow(R, norm=MidpointNormalize(midpoint=0), cmap=cmap)
	else:
		c = ax.imshow(R, cmap=cmap)
	plt.figure(f.number)
	plt.colorbar(c)
	if not (title is None):
		ax.set_title(title)
	plt.tight_layout()
	return f, ax

def tryload(filename, **kwargs):
	try:
		file = np.loadtxt(filename, **kwargs)
	except ValueError:
		#print('tryload: skipping row..')
		file = np.loadtxt(filename, skiprows=1, **kwargs)
	return file

def merge_two_dicts(x, y):
	z = x.copy()   # start with x's keys and values
	z.update(y)    # modifies z with y's keys and values & returns None
	return z

def fit_smail(sample_z, quiet=0, args=(3., 1.8, 1.), bins='auto', weights=None):
	"""
	Fit a Smail-type N(z) to a redshift sample
	Returns 2darray with columns: z, P(z)
	"""
	sampleN, sampleZ = np.histogram(sample_z, bins=bins, density=True, weights=weights)

	def smail(z, alpha, beta, height):
		z0 = np.median(sample_z) / np.sqrt(2.)
		pz = z**alpha * np.exp( -(z/z0)**beta )
		norm = height / pz.max()
		return pz * norm

	popt, _ = curve_fit(smail, midpoints(sampleZ), sampleN, p0=[args[0], args[1], args[2]])
	if not quiet:
		print('alpha: %.2f' % popt[0])
		print('beta: %.2f' % popt[1])
		print('height: %.2f' % popt[2])
	z_grid = np.linspace(sampleZ.min(), sampleZ.max(), int(2e7))
	z_dist = smail(z_grid, *popt)
	z_dist /= z_dist.sum()

	return np.column_stack((z_grid, z_dist))

def data_dict():
	dd = {}
	for k in keys:
		try:
			try:
				dd[k] = np.loadtxt('./wgp_%s'%k)
			except ValueError:
				dd[k] = np.loadtxt('./wgp_%s'%k, skiprows=1)
		except IOError:
			print('wgp_%s not found'%k)
	return dd

def plot_dd(dd, rpp=0.7):
	f, ax = plt.subplots()
	splitk = dict(zip(keys, splitN(6)))
	for k in np.sort(dd.keys()):
		d = dd[k]
		rp, w, err = d.T[0], d.T[1], d.T[3]
		try:
			if all(err == 1.):
				ax.errorbar(rp * splitk[k], rp**rpp * w, label=ltags_wcol[k], c=cols[k], marker='o', capsize=4)
			else:
				ax.errorbar(rp * splitk[k], rp**rpp * w, yerr=rp**rpp * err, label=ltags_wcol[k], c=cols[k], marker='o', capsize=4)
		except KeyError:
			if all(err == 1.):
				ax.errorbar(rp, rp**rpp * w, label=k, marker='o', capsize=4)
			else:
				ax.errorbar(rp, rp**rpp * w, yerr=rp**rpp * err, label=k, marker='o', capsize=4)
	ax.set_xscale('log')
	ax.set_ylabel('$r_{p}^{%s}w_{g+}$'%rpp)
	ax.set_xlabel('$r_{p}$')
	ax.axhline(0., c='k', lw=0.7)
	ax.legend(loc='best')
	plt.tight_layout()
	return f, ax

def concatenate_chains(chains1, chains2):
	mc1 = np.loadtxt(chains1)
	mc2 = np.loadtxt(chains2)
	new_chains = np.concatenate((mc1, mc2))
	np.savetxt(chains2, new_chains, header=open(chains2).readline())
	print("%s += %s --> %s"%(chains1, chains2, chains2))

def plot_jk_zbins(z, zmin, zmax, zcut=None):
    f, ax = plt.subplots()
    pdf, zmid = np.histogram(z, bins=np.linspace(zmin, zmax, 1000))
    pdf = pdf / float(pdf.sum())
    cdf = np.cumsum(pdf)
    zmid = midpoints(zmid)
    ax.plot(zmid, cdf, c='g')

    if zcut == None:
        nze = jk3d.slice_jackknife(z, zmin=zmin, zmax=zmax)
        count = None
    else:
        nze1 = jk3d.slice_jackknife(z, zmin=zmin, zmax=zcut)
        nze2 = jk3d.slice_jackknife(z, zmin=zcut, zmax=zmax)
        nze = np.concatenate((nze1[:-1], nze2))
        count = (np.sum(pdf[(zmid >= zmin) & (zmid < zcut)]), np.sum(pdf[(zmid >= zcut) & (zmid <= zmax)]))
    
    if count != None:
        ax.annotate('%.3f'%count[0], xy=((zcut+zmin)/2., 0.1), xycoords='data', fontsize=12)
        ax.annotate('%.3f'%count[1], xy=((zmax+zcut)/2., 0.1), xycoords='data', fontsize=12)
    for i, bin_edge in enumerate(nze):
        if bin_edge != zmax:
            bin_count = np.sum(pdf[(zmid >= nze[i]) & (zmid < nze[i+1])])
            ax.annotate('%.3f'%bin_count, xy=((nze[i+1]+nze[i])/2., 0.2), xycoords='data', fontsize=11)
        if bin_edge in [zmin, zcut, zmax]:
            ax.axvline(bin_edge, c='r')
        else:
            ax.axvline(bin_edge)

    return f, ax

def JKnorm(jkweights):
	Norm_jk_nom = jkweights.sum() - ( (jkweights**2).sum() / jkweights.sum() )
	Norm_jk_denom = Norm_jk_nom + 1  
	Norm_jack = Norm_jk_nom**2 / Norm_jk_denom
	return Norm_jack

class MidpointNormalize(colors.Normalize):
	def __init__(self, vmin=None, vmax=None, midpoint=None, clip=False):
		self.midpoint = midpoint
		colors.Normalize.__init__(self, vmin, vmax, clip)

	def __call__(self, value, clip=None):
		# I'm ignoring masked values and all kinds of edge cases to make a
		# simple example...
		x, y = [self.vmin, self.midpoint, self.vmax], [0, 0.5, 1]
		return np.ma.masked_array(np.interp(value, x, y))

def weighted_quantile(values, quantiles, sample_weight=None, values_sorted=False, old_style=False):
    values = np.array(values)
    quantiles = np.array(quantiles)
    if sample_weight is None:
        sample_weight = np.ones(len(values))
    sample_weight = np.array(sample_weight)
    assert np.all(quantiles >= 0) and np.all(quantiles <= 1), 'quantiles should be in [0, 1]'
    if not values_sorted:
        sorter = np.argsort(values)
        values = values[sorter]
        sample_weight = sample_weight[sorter]
    weighted_quantiles = np.cumsum(sample_weight) - 0.5 * sample_weight
    if old_style:
        # To be convenient with np.percentile
        weighted_quantiles -= weighted_quantiles[0]
        weighted_quantiles /= weighted_quantiles[-1]
    else:
        weighted_quantiles /= np.sum(sample_weight)
    return np.interp(quantiles, weighted_quantiles, values)

def trace_layout():
	plt.tight_layout()
	plt.subplots_adjust(hspace=0)
	return 0

PHD = "/share/splinter/hj/PhD/"
KIDS = PHD + "KiDS_PhotometricClustering"
EUC = PHD + "cosmosis/modules/euclid_ias/"
MILL = "/share/data1/hj/millennium/"

keys = ['z2_b','z2_r','z1_b','z1_r','sdss_b','sdss_r']
names = dict(zip(keys, ['highZ_Blue', 'highZ_Red', 'lowZ_Blue', 'lowZ_Red']))
file_names = dict(zip(keys, ['highZ_Blue', 'highZ_Red', 'lowZ_Blue', 'lowZ_Red', 'highZ_Blue', 'highZ_Red']))
ltags = dict(zip(keys,['$z>0.26$','$z>0.26$','$z<0.26$','$z<0.26$','SDSS Main','SDSS Main']))
ltags_wcol = dict(zip(keys,['$z>0.26$, blue','$z>0.26$, red','$z<0.26$, blue','$z<0.26$, red','SDSS Main, blue','SDSS Main, red']))
cols = dict(zip(keys, ['blue','red','darkcyan','darkorange', 'steelblue', 'saddlebrown']))
cols2 = dict(zip(keys, ['lightblue','lightcoral','paleturquoise','navajowhite','lightsteelblue', 'rosybrown']))
# for plotting in 2 panels; red & blue
aid = dict(zip(keys, np.array([0,1,0,1,0,1])))
splits = dict(zip(keys, [0.95,0.95,1.,1.,1.05,1.05]))

fisher_plt_params = ['$\Omega_{m}$', '$h$', '$\Omega_{b}$', '$\sigma_{8}$', '$n_{s}$', '$w_{0}$', '$w_{a}$', '$A_{IA, red}$', '$A_{IA, blue}$']
NLA_params = ['$A_{IA, red}$', '$A_{IA, blue}$', '$beta_{red}$', '$beta_{blue}$']
fancy_params = {'omega_m' : r'$\Omega_{\rm{m}}$',
			'h0' : '$h$',
			'omega_b' : r'$\Omega_{\rm{b}}$',
			'sigma8_input' : '$\sigma_{8}$',
			'n_s' : '$n_{s}$',
			'w' : '$w_{0}$',
			'wa' : '$w_{a}$',
			'S8' : r'$S_{8} \equiv \sigma_{8}\sqrt{\Omega_{\rm{m}}/0.3}$',
			'S8_SJ' : r'$S_{8}$',
			'a': r'$A_{\rm{IA}}$',
			'a_1_1': r'$A_{\rm{IA}}^{\rm{R}}$',
			'a_2_2': r'$A_{\rm{IA}}^{\rm{B}}$',
			'ab_1_1': r'$A_{\beta}^{\rm{R}}$',
			'ab_2_2': r'$A_{\beta}^{\rm{B}}$',
			'b_g_z2_r': r'$b_{g}^{\rm{Z2R}}$',
			'b_g_z2_b': r'$b_{g}^{\rm{Z2B}}$',
			'b_g_z1_r': r'$b_{g}^{\rm{Z1R}}$',
			'b_g_z1_b': r'$b_{g}^{\rm{Z1B}}$',
			'b_g_sdss_r': r'$b_{g}^{\rm{SR}}$',
			'b_g_sdss_b': r'$b_{g}^{\rm{SB}}$',
			'ic_z2_r': r'$IC^{\rm{Z2R}}$',
			'ic_z2_b': r'$IC^{\rm{Z2B}}$',
			'ic_z1_r': r'$IC^{\rm{Z1R}}$',
			'ic_z1_b': r'$IC^{\rm{Z1B}}$',
			'ic_sdss_r': r'$IC^{\rm{SR}}$',
			'ic_sdss_b': r'$IC^{\rm{SB}}$',
			'a_z2_r': r'$A_{\rm{IA}}^{\rm{Z2R}}$',
			'a_z2_b': r'$A_{\rm{IA}}^{\rm{Z2B}}$',
			'a_z1_r': r'$A_{\rm{IA}}^{\rm{Z1R}}$',
			'a_z1_b': r'$A_{\rm{IA}}^{\rm{Z1B}}$',
			'a_sdss_r' : r'$A_{\rm{IA}}^{\rm{SR}}$',
			'a_sdss_b': r'$A_{\rm{IA}}^{\rm{SB}}$',
			'a_red': r'$A_{\rm{R}}$',
			'a_blue': r'$A_{\rm{B}}$',
			'a_r': r'$A_{\rm{R}}$',
			'a_b': r'$A_{\rm{B}}$',
			'beta_colour': r'$\beta$',
			'beta_z2_r': r'$\beta_{\rm{Z2R}}$',
			'beta_z2_b': r'$\beta_{\rm{Z2B}}$',
			'beta_z1_r': r'$\beta_{\rm{Z1R}}$',
			'beta_z1_b': r'$\beta_{\rm{Z1B}}$',
			'beta_sdss_r' : r'$\beta_{\rm{SR}}$',
			'beta_sdss_b': r'$\beta_{\rm{SB}}$',
			'beta_1_1': r'$\beta_{\rm{R}}$',
			'beta_2_2': r'$\beta_{\rm{B}}$',
			'beta_r': r'$\beta_{\rm{R}}$',
			'beta_b': r'$\beta_{\rm{B}}$',
			'euclid_red_bias_1': r'$a^{\rm{R}}_{1}$', 
			'euclid_red_bias_2': r'$a^{\rm{R}}_{2}$', 
			'euclid_red_bias_3': r'$a^{\rm{R}}_{3}$', 
			'euclid_red_bias_4': r'$a^{\rm{R}}_{4}$', 
			'euclid_red_bias_5': r'$a^{\rm{R}}_{5}$', 
			'euclid_red_bias_6': r'$a^{\rm{R}}_{6}$', 
			'euclid_red_bias_7': r'$a^{\rm{R}}_{7}$', 
			'euclid_red_bias_8': r'$a^{\rm{R}}_{8}$', 
			'euclid_red_bias_9': r'$a^{\rm{R}}_{9}$', 
			'euclid_red_bias_10': r'$a^{\rm{R}}_{10}$',
			'euclid_blue_bias_1': r'$a^{\rm{B}}_{1}$', 
			'euclid_blue_bias_2': r'$a^{\rm{B}}_{2}$', 
			'euclid_blue_bias_3': r'$a^{\rm{B}}_{3}$', 
			'euclid_blue_bias_4': r'$a^{\rm{B}}_{4}$', 
			'euclid_blue_bias_5': r'$a^{\rm{B}}_{5}$', 
			'euclid_blue_bias_6': r'$a^{\rm{B}}_{6}$', 
			'euclid_blue_bias_7': r'$a^{\rm{B}}_{7}$', 
			'euclid_blue_bias_8': r'$a^{\rm{B}}_{8}$', 
			'euclid_blue_bias_9': r'$a^{\rm{B}}_{9}$', 
			'euclid_blue_bias_10': r'$a^{\rm{B}}_{10}$',
			'kids_red_bias_1': r'$a^{\rm{R}}_{1}$', 
			'kids_red_bias_2': r'$a^{\rm{R}}_{2}$', 
			'kids_red_bias_3': r'$a^{\rm{R}}_{3}$', 
			'kids_red_bias_4': r'$a^{\rm{R}}_{4}$', 
			'kids_red_bias_5': r'$a^{\rm{R}}_{5}$', 
			'kids_red_bias_6': r'$a^{\rm{R}}_{6}$', 
			'kids_red_bias_7': r'$a^{\rm{R}}_{7}$', 
			'kids_red_bias_8': r'$a^{\rm{R}}_{8}$', 
			'kids_red_bias_9': r'$a^{\rm{R}}_{9}$', 
			'kids_red_bias_10': r'$a^{\rm{R}}_{10}$',
			'kids_blue_bias_1': r'$a^{\rm{B}}_{1}$', 
			'kids_blue_bias_2': r'$a^{\rm{B}}_{2}$', 
			'kids_blue_bias_3': r'$a^{\rm{B}}_{3}$', 
			'kids_blue_bias_4': r'$a^{\rm{B}}_{4}$', 
			'kids_blue_bias_5': r'$a^{\rm{B}}_{5}$', 
			'kids_blue_bias_6': r'$a^{\rm{B}}_{6}$', 
			'kids_blue_bias_7': r'$a^{\rm{B}}_{7}$', 
			'kids_blue_bias_8': r'$a^{\rm{B}}_{8}$', 
			'kids_blue_bias_9': r'$a^{\rm{B}}_{9}$', 
			'kids_blue_bias_10': r'$a^{\rm{B}}_{10}$'
			}
exp_zero = [i for i in fancy_params.keys() if ('bias' in i) | ('beta' in i) | i.startswith('a_') | i.startswith('ab_') | i.startswith('ic_')]

def mcplotter():
	vanilla_mcp = mcp.Plots(data_dir=mcp.PHD+'/Vanilla+SDSSdatfiles/', LA_dir=mcp.EUC+'/LA_curves/', LAbeta_dir=mcp.EUC+'/LA_curves_wbeta/',
							NLA_dir=mcp.EUC+'/11rpbin_curves/', beta_dir=mcp.EUC+'/11rpbin_curves_wbeta/')
	HH_mcp = mcp.Plots(data_dir=mcp.PHD+'/HH_datfiles/', LA_dir=mcp.EUC+'/HH_LA_curves/', LAbeta_dir=mcp.EUC+'HH_LA_curves_wbeta/',
						NLA_dir=mcp.EUC+'/HH_curves/', beta_dir=mcp.EUC+'/HH_curves_wbeta/', full=1)
	#rp11mcp = mcp.Plots(IAdata_dir=mcp.PHD+'11rpbin_datfiles/', Clusdata_dir=mcp.PHD+'11rpbin_datfiles/')
	return vanilla_mcp, HH_mcp

def kbplotter():
	reload(clp)
	kb = clp.Plot(data_dirs=['K-Br--GAMA_tomobins/kxg_clus_CORRS/', 'GAMA_tomobins/gama_clus_CORRS/'])
	print("\nf, ax = kb.plot_clustering(rpp=0.8, gap=0.05, ylim=(-1, 3))#, pad=6)")
	return kb

def new_handle(**kwargs):
    return ml.Line2D([], [], **kwargs)

def fancy_labels(cosmo_file, include_params=None, NLA=0, fisher=0, return_params=0, **kwargs):
	print('============\t============\t============\t%s'%cosmo_file)

	params = np.array(cosmosis_params(cosmo_file))

	if type(include_params) == type(None):
		include_params = params.copy()
		#include_params = open(EUC+'cosmo_ia_keys', 'r').read().split('\n')
		#include_params = [i for i in include_params if i!='']

	i = np.array([p in include_params for p in params])
	if any(i==False):
		print('params excluded from plots;')
		[print(missing) for missing in params[~i]]
	else:
		print('plotting all parameters!')

	params = params[i]
	fparams = [fancy_params[p] for p in params]

	data = np.loadtxt(cosmo_file)
	if fisher:
		data = data[ i ][ :, i ]
	else:
		#print('removing Posterior column')
		i = np.array(list(i) + [True]) # deprec: remove P column from MCMC
		data = data[:, i]

	data, fparams = filter_params(data, fparams, fisher=fisher, **kwargs)

	if return_params:
		return data, fparams, params
	else:
		return data, fparams

def filter_params(mcmc, params, fisher=0, ignore_params=None, **kwargs):
	params = np.array(params)
	if type(ignore_params) != type(None):
		cut = np.array( [ (i not in ignore_params) for i in range(len(params)) ] )

		print('filtering %s of (marginalising);' % ['MCMC chains', 'Fisher information matrix'][fisher] )
		[print(i) for i in params[~cut]]

		params = params[cut]
		idx = np.argwhere(cut).flatten()

		if fisher:
			mcmc = mcmc[ cut ][ :, cut]
		else:
			mcmc = mcmc[ :, cut ]

	return mcmc, params

def run_startup():
	from IPython import get_ipython
	ipython = get_ipython()
	ipython.magic("matplotlib")
	ipython.magic("%run ~/.ipython/profile_default/startup/00-start.py")

def splitN(N, gap=0.03, **kwargs):
	#gap = 0.03
	span = gap*(N-1)
	halfspan = span/2.
	split = np.linspace(1.-halfspan, 1.+halfspan, N)
	return split

def ellipse_params(cov, ids, CLsigma=1):
	id1, id2 = ids
	alpha = {1:1.52, 2:2.48, 3:3.44} # ellipse axis scalings corresponding to sigma-confidence levels
	alpha = alpha[CLsigma]
	var1, var2, var12 = cov[id1, id1], cov[id2, id2], cov[id1, id2]
	a = np.sqrt( (var1 + var2)/2. + np.sqrt( 0.25*(var1 - var2)**2 + var12**2 ) )
	b = np.sqrt( (var1 + var2)/2. - np.sqrt( 0.25*(var1 - var2)**2 + var12**2 ) )
	a *= 2 * alpha
	b *= 2 * alpha
	tan2th = np.nan_to_num( (2.*var12)/(var1 - var2) )
	theta = np.arctan(tan2th)/2.
	theta *= (180./np.pi)

	# python ellipse plotting rotates theta[deg] anti-clock
	# therefore a +ve correlation & +ve calculated theta
	# means that a(semi-major) is along the x-axis, and so on;
	R = pearson_r(cov.copy())
	if (R[id1, id2] > 0) == (theta > 0):
		a_width = 1
	else:
		a_width = 0

	#print('a = %.3f'%a)
	#print('b = %.3f'%b)
	#print('theta = %.3f'%theta)
	#print('CLsigma = %i'%CLsigma)
	#print('alpha = %.3f'%alpha)

	return a, b, theta, a_width

def ellipses_2sigma(cov, ids, xy=(-1.019, 0.0), nsigma=2, **kwargs):
	# cov = covariance
	# ids = tuple(index_1, index_2) of the parameters
	# xy = ellipse centre, i.e. point of max likelihood, default is that for w_a(y) - w_0(x) plane
	ells = []
	for i in np.arange(nsigma)+1:
		i = int(i)
		exec 'a%s, b%s, t%s, a_width = ellipse_params(cov, ids, CLsigma=%s)' % (i, i, i, i)
		if a_width:
			exec 'ell%s = Ellipse(xy=xy, width=a%s, height=b%s, angle=t%s)' % (i, i, i, i)
		else:
			exec 'ell%s = Ellipse(xy=xy, width=b%s, height=a%s, angle=t%s)' % (i, i, i, i)
		exec "ells.append(ell%s)" % i
	return ells

def make_plot_grid(params, labsz, ticksz, reverse=0):
    mpl.rcParams['xtick.labelsize'] = ticksz
    mpl.rcParams['xtick.direction'] = 'in'
    mpl.rcParams['xtick.top'] = True
    mpl.rcParams['xtick.major.size'] = 6
    mpl.rcParams['xtick.minor.size'] = 3.0
    mpl.rcParams['xtick.minor.visible'] = True
    mpl.rcParams['ytick.labelsize'] = ticksz
    mpl.rcParams['ytick.direction'] = 'in'
    mpl.rcParams['ytick.right'] = True
    mpl.rcParams['ytick.major.size'] = 6
    mpl.rcParams['ytick.minor.size'] = 3.0
    mpl.rcParams['ytick.minor.visible'] = True
    mpl.rcParams['axes.labelsize'] = labsz
    mpl.rcParams['axes.linewidth'] = 1.5

    #if params[-1] == 'P': params = params[:-1]
    nparam = len(params)
    print('nparam = %s'%nparam)

    f = plt.figure(figsize=(18,18))
    gs = gridspec.GridSpec(nparam, nparam)
    axes = []
    axids = []
    if 'beta_1_1' in params:
		params1 = [p.replace('a', 'ab') if p.startswith('a_') else p for p in params]
		fparams = [fancy_params[i] for i in params1]
    else:
		fparams = [fancy_params[i] for i in params]
    for x in range(nparam):
        for y in range(nparam):

			if reverse:
				#axid = 'ax%s%s' % (str(y).zfill(2), str(x).zfill(2))
				#if 'ax%s%s'% (str(x).zfill(2), str(y).zfill(2)) in axids:
				axid = 'ax%s%s' % ('_%s_'%params[y], '_%s_'%params[x])
				if 'ax%s%s' % ('_%s_'%params[x], '_%s_'%params[y]) in axids:
					continue

				exec '%s = plt.subplot(gs[%i, %i])' % (axid, x, y)
				exec 'plt.setp(%s.get_xticklabels(), visible=False)' %(axid)
				exec 'plt.setp(%s.get_yticklabels(), visible=False)' %(axid)
				exec '%s.xaxis.tick_top()'%axid
				exec '%s.yaxis.tick_right()'%axid
				if x==0:
					exec 'plt.setp(%s.get_xticklabels(), visible=True, rotation=45)' %(axid)
					exec '%s.set_xlabel(r"%s")' % (axid, fparams[y])
				if (y==(nparam-1)) & (x!=(nparam-1)):
					exec 'plt.setp(%s.get_yticklabels(), visible=True)' %(axid)
					exec '%s.set_ylabel(r"%s")' % (axid, fparams[x])
				if x==y:
					exec '%s.tick_params(left=0, right=0)' %(axid)

			else:
				#axid = 'ax%s%s' % (str(x).zfill(2), str(y).zfill(2))
				#if 'ax%s%s'% (str(y).zfill(2), str(x).zfill(2)) in axids:
				axid = 'ax%s%s' % ('_%s_'%params[x], '_%s_'%params[y])
				if 'ax%s%s' % ('_%s_'%params[y], '_%s_'%params[x]) in axids:
					continue

				exec '%s = plt.subplot(gs[%i, %i])' % (axid, y, x)
				if y!=(nparam-1) | (y==(nparam-1)&(x==(nparam-1))):
					exec 'plt.setp(%s.get_xticklabels(), visible=False)' %(axid)
				else:
					exec '%s.set_xlabel(r"%s")' % (axid, fparams[x])
					exec 'plt.setp(%s.get_xticklabels(), rotation=45, visible=True)' %(axid)
				if x!=0 | ((x==0)&(y==0)):
					exec 'plt.setp(%s.get_yticklabels(), visible=False)' %(axid)
				else:
					exec '%s.set_ylabel(r"%s")' % (axid, fparams[y])
				if x==y:
					exec '%s.tick_params(which="both", left=0, right=0)' % axid

			#exec '%s.annotate("%s", xy=(0.6,0.6), xycoords="axes fraction", ha="center", va="center")' % (axid, axid)
			exec 'axes.append(%s)' % (axid)
			axids.append(axid)
    #plt.tight_layout()
    #plt.subplots_adjust(wspace=0, hspace=0)
    return dict(zip(axids, axes))

def find_clevels(Z, *fracs):
# normalise the grid
	if Z.sum() != 1: _Z_ = Z/Z.sum()
	else: _Z_ = Z.copy()

# take cumulative sum
	fZ = _Z_.flatten()
	fZ.sort()
	fZ = fZ[::-1]
	cumsum = np.cumsum(fZ)

# interpolate cumulative probability back to grid-value at
# confidence fraction
	interp = scint.PchipInterpolator(cumsum, fZ)
	levels = [float(interp(frac)) for frac in fracs]
	levels.sort()

	return _Z_, levels

def Posterior2contours(P, bins, axes, fill=1, colors=('steelblue', 'navy', 'cornflowerblue'), linestyle='-', linewidth=1, conf_interval=(0.683, 0.954), **kwargs):
	assert P.ndim == len(bins), 'must pass binning of posterior'
	assert P.ndim == len( [i for i in axes.keys() if i[-2:]==i[-4:-2]] ), 'Posterior ndimensions must match axes'
	dbins = np.array([np.diff(i)[0] for i in bins])
	# re-write this to read a dict pointing to axes?
	for key, ax in axes.iteritems():
		p1, p2 = key.replace('ax','')[:2], key.replace('ax','')[2:]
		p1, p2 = int(p1), int(p2)

		sum_axes = range(P.ndim)

		if p1==p2:
			del sum_axes[p1]

			P_1d = P.sum(axis=tuple(sum_axes))
			P_1d /= P_1d.sum() * np.diff(bins[p1])[0]

			ax.plot(midpoints(bins[p1]), P_1d)

		else:
			del sum_axes[p2]
			del sum_axes[p1]

			Z = P.sum(axis=tuple(sum_axes))
			Z /= Z.sum() * dbins[p1] * dbins[p2]
			X, Y = np.meshgrid(midpoints(bins[p1]), midpoints(bins[p2]))

			CIs = conf_interval
			Z, clevels = find_clevels(Z, *CIs)
			Z = Z.T

			# find x-, y-means
			xcs, ycs = np.cumsum( Z.sum(axis=0) / Z.sum(axis=0).sum() ), np.cumsum( Z.sum(axis=1) / Z.sum(axis=1).sum() )
			xint, yint = scint.interp1d( xcs, midpoints(bins[p1]) ), scint.interp1d( ycs, midpoints(bins[p2]) )

			c = ax.contour(X, Y, Z, clevels, colors=colors[0], linewidths=linewidth, linestyles=linestyle, extend='max', zorder=10)
			ax.plot( xint(0.5), yint(0.5), marker='+', ms=8, mfc='none', mec=colors[1], zorder=10)
			if fill:
				cf = ax.contourf(X, Y, Z, clevels, colors=colors[1:], extend='max')

def chains2contours(x, y, ax, grid=(50,50), weights=None, burn=0, fill=1, colors=('k', 'orangered', 'lightsalmon'),
                              linestyle='-', linewidth=1.5, conf_interval=(0.683, 0.954), alpha=0.5, frac=None, **kwargs):
	w = None
	if type(weights) != type(None):
		print('weighting kde..')
		w = burn_mcmc(weights, burn=burn)

	if frac is not None:
		c = np.random.rand(len(x)) < frac
		if w is not None:
			w = w[c]
		x_ = x[c]
		y_ = y[c]
	else:
		x_ = x.copy()
		y_ = y.copy()
	kdehist = kde.kdehist2(x_, y_, grid, weights=w)
	Z = kdehist[0] / kdehist[0].sum()

	CIs = conf_interval
	Z, clevels = find_clevels(Z, *CIs)
	#clevels = kde.confmap(kdehist[0], [(1-np.exp(-1.)), (1-np.exp(-0.5))])

	clabs = ["%.3f"%i for i in CIs]
	clabs.sort()
	clabs = clabs[::-1]
	#fmt = {}
	#for cl, lab in zip(clevels, clabs):
	#	fmt[cl] = lab

	c = ax.contour(kdehist[1], kdehist[2], Z, clevels, colors=colors[0], linewidths=linewidth, linestyles=linestyle, extend='max', zorder=10)#, alpha=alpha)
	if fill:
		if not all(i=='none' for i in colors[1:3]):
			#import pdb ; pdb.set_trace()
			cf = ax.contourf(kdehist[1], kdehist[2], Z, clevels, colors=colors[1:3], extend='max', alpha=alpha)
			cf = ax.contourf(kdehist[1], kdehist[2], Z, [clevels[1], 1], colors=colors[1:3], extend='max', alpha=alpha)
			#if colors[1]==colors[2]:
		#		cf = ax.contourf(kdehist[1], kdehist[2], Z, clevels, colors=(colors[1], 'none'), extend='max', alpha=alpha)

def posterior_grid(mcmc_file, params=None, axes=None, burn=0.3, fill=1, colors=('w', 'orangered', 'lightsalmon'),
					show_percentiles=1, labelsize=24, ticklabelsize=16, NLA=0, reverse=0, weights=None, qlims=(0.03, 0.97), **kwargs):
	# colors=(linecolor, 2sigma color, 1sigma color)
	# NLA requires more detailed handling of individual sample parameters

	if type(mcmc_file) == str:
		# fancy_labels grabs and filters parameters, and NO LONGER removes P column from mcmc
		# kwargs = ignore_params
		mcmc, params, r_params = fancy_labels(mcmc_file, fisher=0, return_params=1, **kwargs)
		nparams = len(params)
		param_ids = range(nparams)
	elif params == None:
		print("must pass list of parameter strings to posterior_grid if not passing path to mcmc file")
		sys.exit()
	elif type(mcmc_file) == np.ndarray:
		mcmc = mcmc_file
		nsamples, nparams = mcmc.shape[0], len(params)
		param_ids = range(nparams)
		params = [fancy_params[i] for i in params]
	else:
		print("unrecognised mcmc type, must be string=path to file, or numpy ndarray")
		sys.exit()

	mcmc = burn_mcmc(mcmc, burn=burn)
	if type(weights)!=type(None):
		weighted = 1
		w_copy = burn_mcmc(weights.copy(), burn=burn)
	else:
		weighted = 0
		#if len(np.unique(mcmc.T[-1])) == 1:
		w_copy = np.ones_like(mcmc.T[-1])
		#else:
		#	w_copy = 10**mcmc.T[-1]

	if axes == None:
		axes = make_plot_grid(r_params, labsz=labelsize, ticksz=ticklabelsize, reverse=reverse)

	keys = []
	maxlikes = {}
	for i in param_ids:
		#stri = str(i).zfill(2)
		stri = '_%s_'%r_params[i]
		for j in param_ids:
			#strj = str(j).zfill(2)
			strj = '_%s_'%r_params[j]

			if not reverse:
				key = '%s%s' % (stri, strj)
				if '%s%s' % (strj, stri) in keys:
					continue
			else:
				key = '%s%s' % (strj, stri)
				if '%s%s' % (stri, strj) in keys:
					continue

			axk = axes['ax'+key]
			p1, p2 = mcmc.T[i].copy(), mcmc.T[j].copy()
			mx, my = p1.mean(), p2.mean()
			sx, sy = p1.std(), p2.std()

			if i!=j:
				# kwargs = ls, lw, grid, weights
				chains2contours(p1, p2, axk, fill=fill, colors=colors, weights=weights, burn=burn, **kwargs)
				p1_mean = weighted_quantile(p1, 0.5, sample_weight=w_copy)
				p2_mean = weighted_quantile(p2, 0.5, sample_weight=w_copy)
				c1 = colors[3]
				if c1 == 'none':
					c1 = colors[0]
				if c1 == colors[2]:
					c1 = colors[1]
				axk.plot(p1_mean, p2_mean, marker='o', mew=1, ms=6, mfc=c1, mec='w', zorder=10)

				#if (stri[1:-1] in exp_zero) & (strj[1:-1] in exp_zero):
					#axk.plot(0., 0., marker='x', ms=14, mfc='none', mec='k', zorder=10)
				if (stri[1:-1] in exp_zero):
					axk.axvline(0, ls='--', lw=1.5, c='gray', zorder=10)
				if (strj[1:-1] in exp_zero):
					axk.axhline(0, ls='--', lw=1.5, c='gray', zorder=10)

				print('contours %s done'%key)
				if not weighted:
					if not ((qlims[0]==0) & (qlims[1]==0)):
						llow, lhigh = weighted_quantile(p2, qlims)
						axk.set_ylim( llow, lhigh )
					y_ticks, y_tlabs = tick_ranger( axk.get_yticks() )
					axk.set_yticks( y_ticks )
					axk.set_yticklabels( y_tlabs )
				#exec 'axk.annotate("%.3f", xy=(0.4,0.4), xycoords="axes fraction", ha="center", va="center")' % (pearson_r(np.cov(p1,p2))[0,1])
			else:
				#axk.hist(p1, bins=100, normed=True, histtype='step', color=colors[1], weights=w_copy, lw=1.5)
				histc, histb = np.histogram(p1, bins=100, density=1, weights=w_copy)
				c1 = colors[3]
				if c1 == 'none':
					c1 = colors[0]
				axk.step(midpoints(histb), histc/histc.max(), color=c1, lw=2.)

				if (stri[1:-1] in exp_zero):
					axk.axvline(0, ls='--', lw=1.5, c='gray', zorder=10)

				ci_low, ci_mean, ci_high = weighted_quantile(p1, [0.16, 0.5, 0.84], sample_weight=w_copy)
				if show_percentiles:
					for ci in (ci_low, ci_mean, ci_high): axk.axvline(ci, ls='--', c=c1)

				axk.tick_params(which='both', left=0, right=0)

				#if reverse: axk.set_xlabel( r"%s$=%.3f_{-%.3f}^{+%.3f}$"%(params[i], ci_mean, abs(ci_mean-ci_low), abs(ci_high-ci_mean)), fontsize=labelsize)#, y=1.05)
				#else: axk.set_title( r"%s$=%.3f_{-%.3f}^{+%.3f}$"%(params[i], ci_mean, abs(ci_mean-ci_low), abs(ci_high-ci_mean)),
				#					fontsize=labelsize, y=1.05, x=0.05, loc='left')

				#exec 'axk.annotate("%.3f", xy=(0.4,0.4), xycoords="axes fraction", ha="center", va="center")' % p1.mean()

			#mx, sdx = p1.mean(), p1.std()
			#my, sdy = p2.mean(), p2.std()
			#xsd = np.linspace(mx-sdx, mx+sdx, 10)
			#ysd = np.linspace(my-sdy, my+sdy, 10)
			#axk.plot(xsd, np.ones_like(xsd)*my, c='darkred', lw=2)
			#axk.plot(np.ones_like(ysd)*mx, ysd, c='darkred', lw=2)

			#if weighted:
			if not ((qlims[0]==0) & (qlims[1]==0)):
				llow, lhigh = weighted_quantile(p1, qlims)
				axk.set_xlim( llow, lhigh )
			x_ticks, x_tlabs = tick_ranger( axk.get_xticks() )
			axk.set_xticks( x_ticks )
			axk.set_xticklabels( x_tlabs )

			if i == j: maxlikes[key] = p1.mean()
			keys.append(key)

	return plt.gcf(), axes, maxlikes

def fisher_grid(fish_file, maxlikes, apply_priors=0, axes=None, params=None, labelsize=24, ticklabelsize=16, xlims=0, **kwargs):
    # priors should be a dict with on-diagonal keys and stdevs, e.g. {'77' : 0.34}
    # maxlikes should be the same, but with central values, e.g. {'77' : 0.44} - generated by posterior_grid

	if type(fish_file) == str:
		fish, params, r_params = fancy_labels(fish_file, fisher=1, return_params=1, **kwargs)
		nparams = len(params)
		param_ids = range(nparams)
	elif type(params) == type(None):
		print("must pass list of parameter strings to fisher_grid if not passing path to fisher file")
		sys.exit()
	elif type(fish_file) == np.ndarray:
		fish = fish_file
		nparams = len(fish)
		param_ids = range(nparams)
	else:
		print("unrecognised fisher type, must be string=path to file, or numpy ndarray")
		sys.exit()

	if axes == None:
		axes = make_plot_grid(r_params, labsz=labelsize, ticksz=ticklabelsize)

	color = 'royalblue'
	ell_ls = '-'
	ell_lw = 1
	if apply_priors:
		prior_matrix = np.loadtxt(kwargs['priors'])
		assert prior_matrix.shape == fish.shape, "prior - fisher matrix shape mismatch!"
		color = 'darkorange'
		ell_ls = '-.'
		ell_lw = 1.5
		#prior_matrix = np.zeros_like(fish)
		#for i in param_ids:
		#	key = '%s%s' % (i, i)
		#	if key in priors.keys():
		#		prior_matrix[i, i] += 1. / (priors[key]**2.)
		fish += prior_matrix

	cov = np.linalg.inv(fish)
 
	keys = []
	ells = {}
	for i in param_ids:
		#stri = str(i).zfill(2)
		stri = '_%s_'%r_params[i]
		for j in param_ids:
			#strj = str(j).zfill(2)
			strj = '_%s_'%r_params[j]

			key = '%s%s' % (stri, strj)
			if '%s%s' % (strj, stri) in keys:
				continue

			axk = axes['ax'+key]
			x_mean = maxlikes["%s%s" % (stri, stri)]
			y_mean = maxlikes["%s%s" % (strj, strj)]
			x_var = cov[i, i]
			y_var = cov[j, j]

			#xsd = np.linspace(x_mean-x_var**0.5, x_mean+x_var**0.5, 10)
			#ysd = np.linspace(y_mean-y_var**0.5, y_mean+y_var**0.5, 10)
			#axk.plot(xsd, np.ones_like(xsd)*y_mean, c='c')
			#axk.plot(np.ones_like(ysd)*x_mean, ysd, c='c')

			# define 0.01 - 0.99 baseline for pdf
			x_baseline = pdf_baseline(x_mean, x_var**0.5)
			y_baseline = pdf_baseline(y_mean, y_var**0.5)
			pdf = norm.pdf(x_baseline, loc=x_mean, scale=x_var**0.5)

			# normalise pdf
			pdfnorm = np.sum(np.diff(x_baseline)[0] * pdf)
			pdf /= np.trapz(pdf, x_baseline)

			ells[key] = ellipses_2sigma(cov, (i,j), xy=(x_mean, y_mean), **kwargs)

			if i!=j:
				# grid_ells(list_of_ellipse_artists, axis_for_plot)
				grid_ells(ells[key], axk, c=color, a=1, lw=ell_lw, ls=ell_ls)
				axk.set_ylim(minmax(y_baseline))
				y_ticks, y_tlabs = tick_ranger(y_baseline)
				axk.set_yticks( y_ticks )
				axk.set_yticklabels( y_tlabs )
			else:
				axk.plot(x_baseline, pdf, c=color)
				axk.set_title( r"%s$=%.3f\pm{}%.3f$"%(params[i], x_mean, x_var**0.5),
                               fontsize=labelsize, y=1.05, x=0.05, loc='left')#, rotation=20, rotation_mode='anchor')

			if xlims:
				axk.set_xlim(minmax(x_baseline))
				x_ticks, x_tlabs = tick_ranger(x_baseline)
				axk.set_xticks( x_ticks )
				axk.set_xticklabels( x_tlabs )
			#exec 'axk.annotate("%.3f", xy=(0.6,0.2), xycoords="axes fraction", ha="center", va="center")' % (pearson_r(cov.copy())[i, j])
			keys.append(key)

	return plt.gcf(), axes

def plot_trace(mcmc, params=None, nwalkers=50, burn=0.3, idx=None, ax=None, colors=('maroon', 'darkgreen', 'k'), ms=0.2):
	if type(mcmc)==str:
		mcmc, params = cosmosis_mcmc(mcmc, burn=0)
	nparams = len(params)
	for i in range(len(params)):
		if params[i] in fancy_params.keys():
			params[i] = fancy_params[params[i]]

	# calculate Gelman Rubin statistic
	GR = gelman_rubin(mcmc, nwalkers, burn=burn)

	# remove posterior column if nec.
	#if nparams == mcmc.shape[1]-1:
	#	mcmc = mcmc[:, :-1]

	# isolate parameters of interest
	if type(idx) != type(None):
		params, mcmc = np.array(params)[np.array(idx)], mcmc[:, np.array(idx)]
		nparams = len(params)

	# pull chains from flattened array
	nsamples = int(len(mcmc) / nwalkers)
	start = int(np.round(len(mcmc)*burn))
	chains = np.empty([nwalkers, nsamples, len(params)+1])
	for i in range(nwalkers):
		chains[i] = mcmc[i::nwalkers]

	if type(ax) == type(None):
		if nparams%2 == 0: x = 0
		else: x = 1
		f, axe = plt.subplots(nparams//2 + x, 2, sharex=True, figsize=(14, nparams))
		ax = axe.flatten()

	Niter = np.arange(nsamples)
	wb = int(start / nwalkers) # burn per walker
	c1, c2, c3 = colors
	for i in range(nparams):
		for chain in chains:
			ax[i].plot(Niter[:wb], chain[:wb, i], ls='-', marker='', lw=0.2, c=c1)
			ax[i].plot(Niter[wb:], chain[wb:, i], ls='-', marker='', lw=0.2)
		if c3!='w':
			ci_ll, ci_low, ci_mean, ci_high, ci_hh = weighted_quantile(chains[:, wb:, i].flatten(), [.025, .16, .5, .84, .975])#, sample_weight=10**chains[:, wb:, -1].flatten())
			ax[i].axhline(ci_low, ls='--', c='k')
			ax[i].axhline(ci_mean, ls='--', c='k')
			ax[i].axhline(ci_high, ls='--', c='k')
		constraint = r"$%.3f^{+%.3f}_{-%.3f}$" % (ci_mean, abs(ci_high-ci_mean), abs(ci_mean-ci_low))
		#exfrac = frac((chains[:, wb:, i].flatten() > ci_hh)|(chains[:, wb:, i].flatten() < ci_ll))
		#ax[i].annotate(r'%s$\quad{}(ex2sigma=%.3f ,\quad GR=%.3f)$'%(constraint, exfrac, GR[i]), xy=(0.07, 0.07), xycoords='axes fraction', fontsize=15)
		ax[i].annotate(r'%s$\quad{}(GR=%.3f)$'%(constraint, GR[i]), xy=(0.07, 0.07), xycoords='axes fraction', fontsize=15)
		ax[i].set_ylabel(params[i], fontsize=18)
		y_ticks = ax[i].get_yticks()
		y_ticks, y_labs = tick_ranger(y_ticks, nticks=7)
		ax[i].set_yticks( y_ticks )
		ax[i].set_yticklabels( y_labs )
	ax[-1].set_xlabel('N iterations', fontsize=18)
	#f.set_size_inches(20, 10)
	try:
		plt.tight_layout()
	except:
		print('tight_layout() failed')
	plt.subplots_adjust(hspace=0)
	return f, ax
pt = plot_trace

def tick_ranger(data_range, nticks=5, short=0, **kwargs):
	#def myround(x, prec=2, base=.05, **kwargs):
		#return round(base * round(float(x)/base),prec)
	def short_range(mid=0, step=0.05, **kwargs):
		mid_up = np.arange(mid, data_range.max()+step, step=step)
		mid_down = np.arange(data_range.min()-step, mid+step, step=step)
		return np.concatenate((mid_up, mid_down))
	if short:
		ticks = short_range(**kwargs)
		#ticks = [myround(t) for t in ticks]
		tlabs = ["%.2f"%t for t in ticks]
	else:
		ticks = np.linspace(data_range.min(), data_range.max(), nticks)
		tlabs = ["%.3f"%t for t in ticks]
	tlabs = [tl.rstrip('0') for tl in tlabs]
	tlabs[0] = tlabs[-1] = ''
	return ticks, tlabs

def grid_ells(ells, ax, c='c', a=0.3, lw=1.5, ls='-', fc='none'):
	# feed this function with lists/tuples of ellipse artists from
	# ellipses_2sigma(), and an Axes instance on which to plot
	for ell in ells:
		ax.add_artist(ell)
		ell.set_clip_box(ax.bbox)
		ell.set_alpha(a)
		ell.set_facecolor(fc)
		ell.set_edgecolor(c)
		ell.set_linewidth(lw)
		ell.set_linestyle(ls)

def pearson_r(covar_matrix):
	c = covar_matrix
	d = np.diag(c)
	stddev = np.sqrt(d.real)
	c /= stddev[:, None]
	c /= stddev[None, :]
	# Clip real and imaginary parts to [-1, 1]
	np.clip(c.real, -1, 1, out=c.real)
	if np.iscomplexobj(c):
	    np.clip(c.imag, -1, 1, out=c.imag)
	return c

def gelman_rubin(chains, nwalkers, burn=0.3):
	if len(chains.shape) != 3:
		chain = []
		for i in range(nwalkers):
			chain.append(chains[i::nwalkers])
		chain = np.array(chain)
		assert len(chain.shape) == 3, "error reassembling walker chains"
	else:
		chain = chains

	print('chain shape: (%s, %s, %s)' % chain.shape)
	print('burning %s..'%burn)
	remove_burn = int(np.round(chain.shape[1]*burn))
	chain = chain[:, remove_burn:, :]
	print('chain shape: (%s, %s, %s)' % chain.shape)

	ssq = np.var(chain, axis=1, ddof=0)
	W = np.mean(ssq, axis=0)
	tb = np.mean(chain, axis=1)
	tbb = np.mean(tb, axis=0)
	m = float(chain.shape[0])
	n = float(chain.shape[1])
	B = (n / (m - 1)) * np.sum((tbb - tb)**2, axis=0)
	var_t = ((n - 1) / n) * W + (1 / n) * B
	R = np.sqrt(var_t / W)
	return R

#def find_CI(sample, nsig=1):#						IDIOT
#	ci = {1:0.683, 2:0.954, 3:0.997}[nsig]
#
#	# bin likelihood
#	pdf, xbins = np.histogram(sample, bins='auto', density=True)
#	xmid, dx = xbins[:-1]+np.diff(xbins), np.diff(xbins)[0]
#
#	# max-likelihood - better with more samples
#	pdfmax = float(xmid[ np.where(pdf == pdf.max()) ])
#
#	# find y-value of CI integral
#	clipped_pdf_integral = lambda P: pdf[pdf > P].sum()*dx - ci
#	P0 = scipy.optimize.brentq(clipped_pdf_integral, pdf.max(), pdf.min())
#	print('P0 = %.1f'% P0)
#
#	# find corresponding x-values
#	interp = scint.PchipInterpolator(xmid, pdf)
#	findP0 = lambda x: interp(x) - P0
#	xlow = scipy.optimize.brentq(findP0, xmid.min(), pdfmax)
#	xhigh = scipy.optimize.brentq(findP0, pdfmax, xmid.max())

#	print('%ssigma: (-%.3f) %.3f <= %.3f <= %.3f (+%.3f)' % (nsig, (pdfmax-xlow), xlow, pdfmax, xhigh, (xhigh-pdfmax)) )
#	return xlow, pdfmax,  xhigh

def pdfmax(Px_sample):
	pdf, xbins = np.histogram(Px_sample, bins='auto', density=True)
	xmid, dx = xbins[:-1]+np.diff(xbins), np.diff(xbins)[0]
	max = xmid[np.where(pdf == pdf.max())]
	return max

def pdfmean(Px_sample):
	pdf, xbins = np.histogram(Px_sample, bins='auto', density=True)
	xmid, dx = xbins[:-1]+np.diff(xbins), np.diff(xbins)[0]
	mean = np.sum(xmid*pdf*dx)
	return mean

def pdfvariance(Px_sample):
	pdf, xbins = np.histogram(Px_sample, bins='auto', density=True)
	xmid, dx = xbins[:-1]+np.diff(xbins), np.diff(xbins)[0]
	mean = np.sum(xmid*pdf*dx)
	var = np.sum((xmid**2)*pdf*dx) - mean**2
	return var

def cosmosis_params(outfile_path):
	header = open(outfile_path).readline()
	params = header.split('\t')
	params = [param.replace('\n', '') for param in params]
	params = [param for param in params if '--' in param]
	params_ = [param.split('--')[1] for param in params]
	for i in range(len(params)):
		if sum([params_[i]==j for j in params_]) > 1:
			params[i] = 'dummy--'+'_'.join(params[i].split('--'))
	params = [param.split('--')[1] for param in params]
	return params

def cosmosis_mcmc(mcmc_file, burn=0.3):
	print(mcmc_file)
	mcmc = np.loadtxt(mcmc_file)
	mcmc = burn_mcmc(mcmc, burn=burn)
	params = cosmosis_params(mcmc_file)
	return mcmc, params
cmc = cosmosis_mcmc

def cosmosis_mc2dict(mcmc_file, return_P=0, **kwargs):
	mcmc, p = cmc(mcmc_file, **kwargs)
	if return_P:
		mcd = dict(zip(p + ['Post'], mcmc.T))
	else:
		mcd = dict(zip(p, mcmc.T[:-1]))
	keys = np.sort(mcd.keys())
	print(keys)
	return mcd
cmcd = cosmosis_mc2dict

def print_minmax(*cats):
	cats = cats[0]
	for cat in cats:
		ra,dec,z = np.loadtxt(cat).T[:3]
		print(cat)
		print('ra: %.3f - %.3f'%(ra.min(), ra.max()))
		print('dec: %.3f - %.3f'%(dec.min(), dec.max()))
		print('z/chi: %.3f - %.3f'%(z.min(), z.max()))

def burn_mcmc(mcmc, burn=0.3):
	remove_burn = int(np.round(len(mcmc)*burn))
	print('burn = %s ; removing first %s samples..!!' % (burn, remove_burn) )
	mcmc_new = mcmc[remove_burn:]
	return mcmc_new

#print("keys = ['z2_b','z2_r','z1_b','z1_r','sdss_b','sdss_r']\n",
#	"names = dict(zip(keys, ['highZ_Blue', 'highZ_Red', 'lowZ_Blue', 'lowZ_Red']))\n",
#	"file_names = dict(zip(keys, ['highZ_Blue', 'highZ_Red', 'lowZ_Blue', 'lowZ_Red', 'highZ_Blue', 'highZ_Red']))\n",
#	"ltags = dict(zip(keys,['$z>0.26$','$z>0.26$','$z<0.26$','$z<0.26$','SDSS Main','SDSS Main']))\n",
#	"ltags_wcol = dict(zip(keys,['$z>0.26$, blue','$z>0.26$, red','$z<0.26$, blue','$z<0.26$, red','SDSS Main, blue','SDSS Main, red']))\n",
#	"cols = dict(zip(keys, ['blue','red','darkcyan','darkorange','steelblue','maroon']))\n",
#	"cols2 = dict(zip(keys, ['lightblue','lightcoral','paleturquoise','navajowhite','lightsteelblue','rosybrown']))\n",
#	"# for plotting in 2 panels; red & blue\n",
#	"aid = dict(zip(keys, np.array([0,1,0,1,0,1])))\n",
#	"splits = dict(zip(keys, [0.95,0.95,1.,1.,1.05,1.05]))\n")
