from tqdm import tqdm
from numpy import log10
from scipy import stats
from astropy.io import fits
from astropy.table import Table
from astropy.cosmology import FlatLambdaCDM
from scipy.interpolate import interp1d, interp2d, InterpolatedUnivariateSpline
import h5py
import argparse
import warnings
import numpy as np
import pandas as pd
import astropy.units as u
import scipy.optimize as so
fitscol = fits.Column
pchoice = np.random.choice
cosmo = FlatLambdaCDM(H0=70, Om0=0.3)
cmd = lambda x: cosmo.comoving_distance(x).value * cosmo.h
zg = np.linspace(0., 10., 10000)
dg = cmd(zg)
get_d = lambda z: interp1d(zg, dg, bounds_error=0, fill_value=-99.)(z)
get_z = lambda d: interp1d(dg, zg, bounds_error=0, fill_value=-99.)(d)
interp1d_ = lambda x, xx, xy: interp1d(xx, xy, bounds_error=0, fill_value=-99.)(x)
midpoints = lambda x: (x[1:] + x[:-1]) / 2.
sqdeg2ster = lambda a: a / (180./np.pi)**2.
tnorm = stats.truncnorm(-2, 2, loc=0, scale=1)
rtnorm_prep = tnorm.rvs(size=10000)
minmax = lambda x: [x.min(), x.max()]

def smooth(y, box_pts):
	box = np.ones(box_pts)/box_pts
	y_smooth = np.convolve(y, box, mode='same')
	return y_smooth

def find_nearest(array, values, give_indices=False):
	array = np.asarray(array)
	# the last dim must be 1 to broadcast in (array - values) below.
	values = np.expand_dims(values, axis=-1)
	indices = np.abs(array - values).argmin(axis=-1)
	if give_indices:
		return indices
	else:
		return array[indices]

def get_k_z(kcorrs):
	"""
	Read maggy ratios from file, where each can be converted
	to a k-correction from z (given by column) to z=0.
	Return:
		1. k-correction matrix: shape = (N galaxies x len(z_grid))
		2. z_grid: at which the k(z->0) is to computed
		3. mask: indicates bad photometry
	"""
	#df_kcorrs = pd.read_csv(kcorrs, delimiter=' ')
	df_kcorrs = h5py.File(kcorrs, 'r')
	ID = df_kcorrs['ID'][:]
	z = df_kcorrs['z'][:]
	maggyratio = df_kcorrs['maggy_ratio'][:]
	idsort = np.argsort(ID[:len(set(ID))])

	zero_ratio = maggyratio[np.where(z == 0)[0]] # mgy(z=0) / mgy(z=0)
	mask = zero_ratio != 1 # ratio != 1 indicates error in computation of kcorrections -- should trace back to bad photometry

	z_grid = np.unique(z)
	uID = np.unique(ID)
	maggyratio_matrix = np.array(maggyratio).reshape(len(z_grid), len(uID)).T
	maggyratio_matrix[mask] = 1e-1
	k_z = -2.5 * log10(maggyratio_matrix)
	return k_z[idsort], z_grid, mask[idsort]

def interp_k(z_grid, k_z, z_out):
	#k_out = []
	#for i in range(len(k_z)):
	#	k_out.append( InterpolatedUnivariateSpline(z_grid, k_z[i], bbox=[0., z_grid.max()], ext='const')(z_out[i]) )
	#return np.array(k_out)
	print '\tinterpolating k-corrections..'
	x = np.array([interp1d(z_grid, k_z[i], #axis=1,
				assume_sorted=True, copy=False,
				bounds_error=False, fill_value=z_grid.max())(z_out) for i in range(len(k_z))])
	return x

def interp_dm(z_grid, d_grid, z_out):
	#print '\tinterpolating distance moduli..'
	x = interp1d(z_grid, d_grid,
				assume_sorted=True, copy=False,
				bounds_error=False, fill_value=z_grid.max())(z_out)
	return x

def fit_zmax(fluxlim, m_obs, z_obs, z_grid, k_z, min=0):
	"""
	Fit a maximum redshift to each galaxy
	given the observed redshift and magnitude,
	the flux-limit FL, and the k(z).

	Arguments:
		fluxlim: limiting magnitude of survey (float)
		m_obs: observed magnitudes (1d-array of floats, length = # of galaxies)
		z_obs: observed redshifts (as m_obs)
		z_grid: redshifts at which k-corrections are provided (1d-array of floats)
		k_z: k-corrections as functions of redshift (2d-array of floats, shape = N_galaxies x len(z_grid))

	Returns:
		z_max: maximum observable redshift per object (as m_obs & z_obs)

	"""
	# get observed distance moduli
	z_grid_fine = np.linspace(z_grid.min(), z_grid.max(), 1000)
	d_grid_fine = cosmo.distmod(z_grid_fine)
	dm_obs = interp_dm(z_grid_fine, d_grid_fine, z_obs)

	# build k-correction interpolator per object
	print '\tbuilding k-correction interpolators..'
	i1d_args = dict(assume_sorted=True, copy=False,
                	bounds_error=False, fill_value=0.)
	k_interpolators = np.array([interp1d(z_grid, k_z[i], **i1d_args) for i in range(len(k_z))])

	# get k-corrections from z_obs -> z=0
	k_obs = np.array([ki(z) for ki, z in np.column_stack((k_interpolators, z_obs))])

	# define the function to minimise
	LHS = fluxlim - m_obs + dm_obs + k_obs
	def distmod_relation(z_max, idx, get_kmax=1):
		# LHS = mu(z_max) + k(z_max)
		# the array z_max has length = # of galaxies
		# pass idx to specify indices of galaxies within main array

		dm_max = interp_dm(z_grid_fine, d_grid_fine, z_max)
		if get_kmax:
			k_max = np.array([ki(z) for ki, z in np.column_stack((k_interpolators[idx], z_max))])
		else:
			k_max = np.zeros_like(z_max)

		loss_fn = LHS[idx] - (dm_max + k_max)
		loss_fn[np.isnan(loss_fn) | np.isinf(loss_fn)] = 0.
		return np.abs(loss_fn)

	if not min:
		print '\tfitting z_max per galaxy..'
	else:
		print '\tfitting z_min per galaxy..'
	x = minimize(distmod_relation, z_obs, bounds=[0., z_grid.max()])
	return x

def minimize(fun, x0, tol=1e-2, bounds=[-np.inf, np.inf], quiet=True):
	# find variables x that bring fun(x) down to < tol[mags]
	# starting with guess x0
	# perturb x0 by +/- 10% -> x1, and evaluate fun(x1)
	# if fun(x1) < fun(x0), perturb x1 by 10% again,
	# but now with a smaller chance to flip the sign
	# if fun(x1) > fun(x0), perturb x1 by 10% again,
	# but now with a larger chance to flip the sign
	nflip = np.array([-1.,-1.,-1.,-1.,1.])
	pflip = np.array([1.,1.,1.,1.,-1.])
	def perturb(y, sign=None, boost=None, bounds=None):
		# return y1 = y +/- 10%, the sign of the shift and the perturbation
		perturbation = (np.random.rand(len(y)) - 0.5) / 5.
		if sign is not None:
			# impose sign, if given
			perturbation = abs(perturbation) * sign
		if boost is not None:
			# in/deflate perturbation, if given
			perturbation *= boost
		sign = np.sign(perturbation)
		y1 = y * (1. + perturbation)

		if bounds is not None:
			y1[y1 < bounds[0]] = bounds[0]
			y1[y1 > 10.] = bounds[1]

		return y1, sign, perturbation

	# initially shift all outward, since starting with observed z
	x, s, p = perturb(x0, sign=np.ones_like(x0), bounds=bounds)
	x_out = np.zeros_like(x0)
	idx = np.arange(len(x0), dtype=int)
	idx1 = idx.copy()
	get_kmax = True # True=slower but strictly more accurate, will assume kmax=0 otherwise
	fx = fun(x, idx, get_kmax=get_kmax)
	fx0 = fun(x0, idx, get_kmax=get_kmax)

	fins = []
	while any(x_out == 0):
		fin = abs(fx) < tol
		fins.append(fin.sum())

		# catch brightest galaxies with z_max >> survey limit
		if all(np.array(fins[-100:]) == 0):
			print '\trounding-off bright galaxies..'
			fracs = np.linspace(0.01, 1., 200)
			roundoff_array = np.array([fun(x*frac, idx1) for frac in fracs])
			minima = np.argmin(roundoff_array, axis=0)
			x_out[idx1] = x * fracs[minima]
			break

		# remove successfully shifted galaxies
		x_out[idx1[fin]] = x[fin]
		x = x[~fin]
		idx1 = idx1[~fin]
		fx = fx[~fin]
		fx0 = fx0[~fin]
		s = s[~fin]

		gain = (fx < fx0) # closer to z_max == unlikely to change shift direction
		loss = (fx > fx0) # further from z_max == very likely to change shift direction
		s[gain] *= np.random.choice(pflip, size=gain.sum())
		s[loss] *= np.random.choice(nflip, size=loss.sum())

		boost = np.ones_like(fx)
		boost[gain] -= np.random.rand(gain.sum()) / 2. # if closer, decrease the shift by up to 50%
		boost[loss] += np.random.rand(loss.sum()) / 2. # if further, increase the shift by up to 50%

		fx0 = fx.copy()
		x, s, p = perturb(x, sign=s, boost=boost, bounds=bounds)
		# restart frozen entries
		if all(np.array(fins)[-5:] == 0):
			x[x < x0[idx1]] = x0[idx1][x < x0[idx1]]
			s[x < x0[idx1]] = 1.
		fx = fun(x, idx1, get_kmax=get_kmax)

		if not quiet:
			try:
				for i in np.random.choice(range(len(x)), size=10):
					print 'galaxy %s: delta=%.3f[mag], last z_max=%.3f, z of object=%.3f' % (i+1, fx[i], x[i], x0[idx1][i])
				print 'N outside tolerance =', (~fin).sum()
			except:
				pass

	return x_out

def get_volume_limits(d, area=180., volume=3.5e6):
	# get asymmetric LoS limits for a given
	# volume centred on each galaxy with distance d
	Om = sqdeg2ster(area)

	def invvol(r1, vol):
		r2 = ((3. * vol / Om) + r1**3.) ** (1./3.)
		# calculate asymmetric limits of volume -- see notes
		with warnings.catch_warnings():
			warnings.simplefilter("ignore")
			new_r1 = ((3.*r1**3. - r2**3.) / 2.) ** (1./3.)
		new_r2 = ((r1**3. + r2**3.) / 2.) ** (1./3.)
		# fix NaNs in new_r1 -- make positive for the cube-root, and then negative again
		badr1 = np.isnan(new_r1)
		try:
			new_r1[badr1] = (-(3.*r1[badr1]**3. - r2[badr1]**3.) / 2.) ** (1./3.)
		except IndexError:
			new_r1[badr1] = (-(3.*r1**3. - r2[badr1]**3.) / 2.) ** (1./3.)
		new_r1[badr1] *= -1.
		# for symmetric limits (will override the above):
		#new_r1 = 2.*r1 - r2
		#new_r2 = r2
		return [new_r1, r1, new_r2]

	widths = np.asarray(invvol(d, volume)).squeeze()

	return widths

def get_volume_depth(volume, area=180.):
	Om = sqdeg2ster(area)
	return (3. * volume / Om)**(1./3.)

def get_volume(limits, area=180.):
	# get volume defined by limits, for a given area
	Om = sqdeg2ster(area)
	volume = (Om / 3.) * (limits[1]**3. - limits[0]**3.)
	return volume

def get_tgauss_window(d, hilim, lolim, area=180., volume=3.5e6, d_base=None, zmax=0.6):
	Om = sqdeg2ster(area)
	dmax = get_d(zmax)
	mid_bins = midpoints(d_base)
	window_fn = np.zeros_like(mid_bins)

	if d > dmax or hilim == -99. or lolim == -99. or lolim == hilim:
		# throw away those with poor k-corrections->limits
		window_fn *= np.nan
	elif d < 10:
		# flat window for nearby galaxies
		window_fn[(mid_bins >= lolim) & (mid_bins <= hilim)] = 1.
	else:
		# draw from a truncated [-2, 2] Gaussian (pre-drawn); convert to volumes
		rvols = np.abs(rtnorm_prep) * volume

		# convert draws to comoving distances
		limits = get_volume_limits(d, volume=rvols, area=area)[[0, 2]]
		limits = np.concatenate(limits)
		# reflect in boundaries
		while any(limits > hilim) or any(limits < lolim):
			limits[limits < lolim] += 2.*(lolim - limits[limits < lolim])
			limits[limits > hilim] += 2.*(hilim - limits[limits > hilim])

		# bin, smooth and normalise for the pdf
		if d_base is not None: bins = d_base
		else: bins = 'auto'
		window_fn, bin_edges = np.histogram(limits, bins=bins)
		mid_bins = midpoints(bin_edges)
		window_fn = smooth(window_fn, 5)
		window_fn[mid_bins < lolim] = 0.
		window_fn[mid_bins > hilim] = 0.

	# normalisation must be such that int[W(V) dV] = 1
	N = Om * np.trapz(window_fn * mid_bins**2., x=mid_bins)
	window_fn = np.asarray(window_fn, dtype=np.float32) / N
	#window_fn *= Om * mid_bins**2.

	window_fn = np.stack((mid_bins, window_fn))
	return window_fn

def clone_galaxies(idcol, maxcol, Nrand=10, zlims=None, window_vol=None, area=180.,
				   dobs=None, dres=1., Niter=15, load_windows=True, runid='', mincol=None):
	# take arrays of IDs and max distances/redshifts
	# draw Nrand redshifts from quadratic/fitted distribution [0, {d/z}max] per galaxy
	Om = sqdeg2ster(area)
	Ngal = len(idcol)
	zmin, zmax = zlims or (0., 0.6)
	dmin, dmax = cmd([zmin, zmax])
	print '\t\t(', (dobs > maxcol).sum(), 'objects observed beyond calculated maximum -- REMOVING )'
	#maxcol[dobs > maxcol] = dobs[dobs > maxcol]
	maxcol[dobs > maxcol] = -99.
	maxcol[maxcol < dmin] = -99.
	if mincol is not None:
		print '\t\t(', (dobs < mincol).sum(), 'objects observed closer than bright limit -- REMOVING )'
		#mincol[dobs < mincol] = dobs[dobs < mincol]
		mincol[dobs < mincol] = -99.
	else:
		mincol = dmin * np.ones_like(maxcol)

	# get limits; zmax/zmin or survey edges
	upps = np.stack((maxcol, np.ones_like(maxcol)*dmax))
	lows = np.stack((np.zeros_like(maxcol), mincol))
	hilim = np.min(upps, axis=0)
	lolim = np.max(lows, axis=0)

	# ready comoving distance baseline
	dres_fine = dres / 5.
	d_base = np.arange(dmin, dmax+dres, step=dres)
	d_mid = midpoints(d_base)
	d_mid_fine = np.arange(dmin, dmax+dres, step=dres_fine)
	dd_mid = np.diff(d_mid)[0]
	dd_mid_fine = np.diff(d_mid_fine)[0]

	# take galaxy n(chi)
	n_g = np.asarray(np.histogram(dobs, bins=d_base)[0], dtype=np.float32)

	if window_vol is not None:
		print '\t\t( 1sigma window volume given as %.1e [Mpc/h]^3 )'%window_vol
		if load_windows:
			print '\t\tloading Gaussian windows..'
			try:
				f = h5py.File('windows%s.h5'%runid, 'r')
				windows = f['windows'][:]
				f.close()
				load_failed = False
			except IOError:
				print '\t\twindow-loading failed;'
				load_failed = True
		if (not load_windows) or load_failed:
			# get truncated Gaussian window per object
			windows = [get_tgauss_window(dobs[i], hilim[i], lolim[i], volume=window_vol, d_base=d_base, zmax=zmax, area=area)
						for i in tqdm(range(len(dobs)), desc='\t\twindowing', ncols=100)]
			windows = np.array(windows)[:, 1, :]
			assert (not all(np.isnan(windows.flatten()))), "Windowing failed! All windows are NaN; d > dmax or bad limits"

		with h5py.File('windows%s.h5'%runid, 'w') as f:
			f.create_dataset('windows', data=windows)
			f.close()

		# interpolate windows onto finer comoving grid
		windows_fine = np.empty([len(windows), len(d_mid_fine)])
		for i, w in enumerate(tqdm(windows, desc='\t\tinterpolating windows', ncols=100)):
			window_int = np.interp(d_mid_fine, d_mid, w)
			windows_fine[i] = window_int

		# combine window with pdf
		print '\t\tnormalising windows/building pdfs..',
		windows_fine = (windows_fine.T / np.sum(windows_fine * dres_fine, axis=1)).T
		#windows = (windows.T / np.sum(windows * dres, axis=1)).T
#		quadratic_weight = np.array([d_mid_fine**2. / np.sum(d_mid_fine[w!=0]**2. * dres_fine)
#										for w in tqdm(windows_fine, desc='\t\tgetting quad-weights', ncols=100)])
#		pdfs = np.array([windows_fine[i] * quadratic_weight[i]
#							for i in tqdm(range(len(windows_fine)), desc='\t\tapplying weights', ncols=100)])
		pdfs = (windows_fine.T / np.sum(windows_fine, axis=1)).T
		print 'done.'

	else:
		windows = np.ones([len(maxcol), len(d_mid)], dtype=np.float32)
		for i in tqdm(range(len(windows)), desc='\t\tlimiting integrals', ncols=100):
			windows[i][(d_mid < lolim[i]) | (d_mid > hilim[i])] = 0.
		norm = Om * np.trapz(windows * d_mid**2., x=d_mid, axis=1)
		windows = (windows.T / norm).T

	# get Vmax
	Vmax = Om * np.trapz(d_mid**2. * windows, x=d_mid)
	if window_vol is not None:
		rawV = get_volume([lolim, hilim], area=area)
		small_Vmax = rawV < window_vol
		print '\t\trestricting small-Vmax (< window) objects: %s'%small_Vmax.sum()
		small_windows = [get_tgauss_window(dobs[small_Vmax][i], hilim[small_Vmax][i], lolim[small_Vmax][i],
											 volume=rawV[small_Vmax], d_base=d_base, zmax=zmax, area=area)
								for i in tqdm(range(len(dobs[small_Vmax])), desc='\t\tre-windowing', ncols=100)]
		windows[small_Vmax] = np.array(small_windows)[:, 1, :]

	# setup diagnostic save-outs
	Vmax_dc_list = []
	n_clones_list = []
	Delta_d_list = []
	pdf_list = []
	n_r_list = []
	
	this_iter = 1
	print '\t\titerating n_clones calculation..'
	while this_iter < Niter+1:
		print '\t\t======= ITERATION #%s =======' % this_iter
		# get Vmax,dc ; density-weighted Vmax for n_clone calculation
		if this_iter == 1:
			Delta_d = np.ones_like(d_mid, dtype=np.float64)
		else:
			n_r = np.asarray(np.histogram(ddraw, bins=d_base)[0], dtype=np.float32)
			Delta_d = np.nan_to_num(Nrand * n_g / n_r)
		delta_mask = ~np.isnan(Delta_d) & ~np.isinf(Delta_d)
		#assert all(delta_mask), "Delta going undefined!"
		Vmax_dc = Om * np.trapz(d_mid[delta_mask]**2. \
								* Delta_d[delta_mask] \
								* windows.T[delta_mask].T, \
								x=d_mid[delta_mask], axis=1)
		Vmax_dc = Vmax / Vmax_dc

		Vmax_dc_list.append(Vmax_dc)

		bad_Vmaxdc = np.isnan(Vmax_dc) | np.isinf(Vmax_dc) | (Vmax_dc <= 0)
		if this_iter > 1:
			prev_N = n_clones.sum()
		n_clones = np.asarray(np.round(Nrand * Vmax / Vmax_dc), dtype=int)
		n_clones[bad_Vmaxdc] = 0
		n_clones[n_clones < 0] = 0
		n_clones_list.append(n_clones)
		print '\t\ttotal N clones = %s'%n_clones.sum()
		if this_iter > 1:
			print '\t\t\t\ti.e. %+d'%(n_clones.sum() - prev_N)

		Nrand_idcol = np.repeat(idcol, n_clones)
		Nrand_dobs = np.repeat(dobs, n_clones)
		Nrand_maxcol = np.repeat(maxcol, n_clones)
		Nrand_Vmax = np.repeat(Vmax, n_clones)
		Nrand_lolim = np.repeat(lolim, n_clones)
		Nrand_hilim = np.repeat(hilim, n_clones)

		Nrand_maxcol[Nrand_dobs > Nrand_maxcol] = Nrand_dobs[Nrand_dobs > Nrand_maxcol]
		badmax = ((Nrand_maxcol == -99.) |
				  (Nrand_hilim == -99.) |
				  (Nrand_lolim == -99.) |
				  (Nrand_dobs < 30.) | (Nrand_maxcol < 30.) |
				  (Nrand_dobs > dmax) |
				  (Nrand_lolim == Nrand_hilim))

		if window_vol is not None:
			# draw from truncated [-2, 2] Gaussian
			#tnorm = stats.truncnorm(-2, 2, loc=0, scale=1) # this draw needs to be weighted as chi^2?
			#rtnorm = tnorm.rvs(size=len(Nrand_dobs))
			# convert draws to window-volumes
			#rvols = np.abs(rtnorm) * window_vol
			# calculate limits for the volumes and take as the distance-draws
			#lims = get_volume_limits(Nrand_dobs, volume=rvols)
			#ddraw = np.where(rtnorm < 0, lims[0], lims[2])

			# draw n_clones from each comoving distance probdens function
			ddraw = []
			for pdf, nc in tqdm(zip(pdfs, n_clones), desc='\t\tfilling windows', ncols=100):
				if nc != 0:
					draw = np.random.choice(d_mid_fine, p=pdf, size=nc) \
							+ (np.random.rand(nc) - 0.5) * dres
					ddraw.append(draw)
				else:
					continue
			ddraw = np.concatenate(ddraw)

			ddraw[badmax] = -99. # set unwanted to -99

			# reflect out-of-bounds draws in the boundaries
			print '\t\treflecting at boundaries..'
			oob = ~badmax & ((ddraw < Nrand_lolim) | (ddraw > Nrand_hilim))
			oob1 = len(oob)
			while any(oob):
				ddraw[ddraw < Nrand_lolim] += 2.*(Nrand_lolim[ddraw < Nrand_lolim] - ddraw[ddraw < Nrand_lolim])
				ddraw[ddraw > Nrand_hilim] += 2.*(Nrand_hilim[ddraw > Nrand_hilim] - ddraw[ddraw > Nrand_hilim])
				oob = ~badmax & ((ddraw < Nrand_lolim) | (ddraw > Nrand_hilim))
				if oob.sum() < oob1 and oob.sum() != 0:
					print oob.sum()
					oob1 = oob.sum()
		else:
			# draw from d^2 (== growth of volume element) between [0, upperlimit]
			ddraw = np.random.power(3, n_clones.sum()) * Nrand_hilim
			while any(ddraw[~badmax] < Nrand_lolim[~badmax]):
				ddraw[ddraw < Nrand_lolim] = np.random.power(3, (ddraw < Nrand_lolim).sum()) * Nrand_hilim[ddraw < Nrand_lolim]
			#vdraw = np.random.rand(n_clones.sum()) * Nrand_Vmax
			#ddraw = get_volume_depth(vdraw, area=area)
			ddraw[badmax] = -99.

		Delta_d_list.append(Delta_d)
		try: pdf_list.append(pdf)
		except: pass
		try: n_r_list.append(n_r)
		except: pass
		this_iter += 1

	with h5py.File('diagnostics_clonerandoms%s.h5'%runid, 'w') as h5_diag:
		h5_diag.create_dataset('Vmax', data=Vmax)
		h5_diag.create_dataset('d_mid', data=d_mid)
		h5_diag.create_dataset('Vmax_dc', data=Vmax_dc_list)
		h5_diag.create_dataset('n_clones', data=n_clones_list)
		h5_diag.create_dataset('Delta_d', data=Delta_d_list)
		try: h5_diag.create_dataset('pdf', data=pdf_list)
		except: pass
		h5_diag.create_dataset('n_r', data=n_r_list)
		h5_diag.close()

	mask = ddraw > 0.
	zdraw = np.ones_like(ddraw) * -99.
	zdraw[mask] = get_z(ddraw[mask])

	randoms_id_z = np.column_stack((Nrand_idcol, zdraw))
	return randoms_id_z

#def mask_randoms(args):
	# expand randoms from 1x1 sky-box into catalogue footprint
	# create mutliple boxes for disjoint footprints?
	# use jackknife code to identify these?
	# apply a mask to resulting window
# OR
	# write code to sample from fitted 2D ra/dec distributions?

def main(args):
	print '\n'
	for cat_path in args.catalogues:
		print 'reading %s..' % cat_path
		cat = fits.open(cat_path)[1].data
		t = Table(cat)

		if args.refresh_zmax:
			# establish zmax per detection band
			for mcol, mlim, kcorr, blim in zip(args.magcols, args.maglims, args.kcorrs, args.brightlim):
				print '\t"%s" at limit: %s' % (mcol, mlim)
				print '\tbright limit: %s' % blim
				print '\treading "%s" k-corrections..' % kcorr

				k_z, z_grid, mask = get_k_z(kcorr)
				max_redshift = fit_zmax(mlim, cat[mcol], cat[args.zcol], z_grid, k_z)
				max_redshift[mask] = -99.

				zmax_col = mcol+'_fl%.1f_zmax'%mlim
				t[zmax_col] = max_redshift

				if blim is not None:
					zmin_col = mcol+'_fl%.1f_zmin'%blim
					min_redshift = fit_zmax(blim, cat[mcol], cat[args.zcol], z_grid, k_z, min=1)
					t[zmin_col] = min_redshift

			t.write(cat_path, format='fits', overwrite=1)
		cat = t.copy()

		if args.randoms:
			print '\tbuilding clone randoms..'
			# drop points at random in 1x1 square
			Nrand = len(cat) * args.Nrand
			idcol = cat[args.idcol]
			zcol = np.asarray(cat[args.zcol])
			mask = (zcol > 0)

			# apply any weighting
			if args.zweight is None:
				wcol = None
				bins = 'auto'
			else:
				wcol = cat[args.zweight]
				bins = np.linspace(zcol[mask].min(), zcol[mask].max(), 40)

			# construct redshift distribution per band
			#print '\t\tcloning..'
			random_cols = []
			for mcol, mlim, blim in zip(args.magcols, args.maglims, args.brightlim):
				zmax_col = mcol+'_fl%.1f_zmax'%mlim
				#zmax_col = 'zmax_19p8'

				maxcol = get_d(cat[zmax_col])
				dobs = get_d(zcol)
				if blim is not None:
					zmin_col = mcol+'_fl%.1f_zmin'%blim
					mincol = get_d(cat[zmin_col])
				else:
					mincol = None
				randoms_id_z = clone_galaxies(idcol, maxcol, args.Nrand, zlims=args.zlims, window_vol=args.window, area=args.area, dres=args.dres,
											  dobs=dobs, Niter=args.niter, load_windows=args.load_windows, runid=args.id, mincol=mincol)
				randoms_comoving = get_d(randoms_id_z[:, 1])

				random_cols.append(fitscol(array=np.asarray(randoms_id_z[:, 0], dtype=np.int32), name=mcol+'_cloneID', format='K'))
				random_cols.append(fitscol(array=randoms_id_z[:, 1], name=mcol+'_cloneZ', format='D'))
				random_cols.append(fitscol(array=randoms_comoving, name=mcol+'_cloneComovingDist', format='D'))

			ralims = minmax(cat[args.colnames[0]])
			declims = minmax(cat[args.colnames[1]])
			dec = np.arange(declims[0], declims[1]+0.1/60., step=0.1/60.)
			cos = np.cos(np.deg2rad(dec))
			cos /= cos.sum()

			ra = np.random.rand(len(randoms_id_z)) * np.diff(ralims) + ralims[0]
			dec = np.random.choice(dec, p=cos, size=len(randoms_id_z))

			random_cols.append(fitscol(array=ra, name='ra', format='D'))
			random_cols.append(fitscol(array=dec, name='dec', format='D'))

			# save randoms catalogue with galaxy IDs and redshifts for selection
			random_cols = fits.ColDefs(random_cols)
			rand_hdu = fits.BinTableHDU.from_columns(random_cols)
			if args.o is None:
				out = cat_path.replace('.fits', '_CloneZIDRandoms.fits')
			else:
				out = args.o
			if args.id is not None:
				out = out.replace('.fits', '%s.fits'%args.id)
			rand_hdu.writeto(out, overwrite=1)

	print 'done!'

if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument(
		'catalogues',
		nargs='*',
		default='give paths to fits catalogue(s) with redshifts and absolute rest-frame magnitudes')
	parser.add_argument(
		'-id',
		type=str,
		default='',
		help='give run-identifier for outputs, default is a blank string')
	parser.add_argument(
		'-magcols',
		nargs='*',
		type=str,
		help='specify observer-frame magnitude column name(s) from which to compute zmax(s)')
	parser.add_argument(
		'-maglims',
		nargs='*',
		type=float,
		help='specify apparent magnitude detection limit per -magcols column')
	parser.add_argument(
		'-brightlim',
		nargs='*',
		type=float,
		help='specify apparent magnitude *minimum* (bright limit) per -magcols column')
	parser.add_argument(
		'-kcorrs',
		nargs='*',
		type=str,
		help='specify k-corrections filename(s), for computation of zmax(s). File must have columns ID, z, maggie_ratio')
	parser.add_argument(
		'-randoms',
		type=int,
		help='1 = create a randoms catalogue with 1 redshift distribution per -magcol')
	parser.add_argument(
		'-dres',
		type=float,
		default=1.,
		help='optionally, specify grid resolution (float) for window calculation, default= 1 Mpc/h')
	parser.add_argument(
		'-window',
		type=float,
		help='optionally, give 1s window volume for windowed scattering of randoms (e.g. GAMA is 3.5e6)')
	parser.add_argument(
		'-load_windows',
		type=int,
		default=1,
		help='1 = load per-object truncatedw Gaussian window functions (default), or 0 = re-compute them')
	parser.add_argument(
		'-refresh_zmax',
		type=int,
		default=0,
		help='1 = re-compute z_max for each galaxy; do so if k-corrs/mag-lims/bright-lims have changed')
	parser.add_argument(
		'-niter',
		type=int,
		default=15,
		help='optionally, give number of times to iterate the calculation of n_clone per object, default=15')
	parser.add_argument(
		'-idcol',
		default='ID',
		type=str,
		help='specify galaxy identifier column-name -- default=ID')
	parser.add_argument(
		'-zcol',
		default='z',
		type=str,
		help='specify galaxy redshift column-name -- default=z')
	parser.add_argument(
		'-colnames',
		type=str,
		default=['RA','DEC'],
		nargs=2,
		help='specify galaxy ra/dec column names, default= RA DEC')
	parser.add_argument(
		'-area',
		type=float,
		default=180.,
		help='specify survey area for volume calculations, default=180 i.e. GAMA')
	parser.add_argument(
		'-zlims',
		nargs=2,
		default=[0., 2.],
		type=float,
		help='specify redshift window, s/t no randoms will be dropped outside of this range -- default=[0, 2]')
	parser.add_argument(
		'-zweight',
		type=str,
		help='specify galaxy redshift weighting for n(z) fitting')
	parser.add_argument(
		'-Nrand',
		default=10,
		type=int,
		help='specify integer multiple of catalogue galaxies for cloning -- default=10')
	parser.add_argument(
		'-o',
		type=str,
		help='specify output filename for randoms -- default is <catalogue>_CloneRandoms.fits')
	args = parser.parse_args()

	#if args.kcorrs is None:
	#	args.kcorrs = [None]*len(args.magcols)
	if args.id != '':
		args.id = '_'+args.id
	if args.brightlim is None:
		args.brightlim = [None]*len(args.maglims)
	main(args)





#def dmax(lim, M_, ceiling=np.inf, Mcorr=None):
#	# maximum distance of absolute M galaxy in Mpc, given apparent mag lim
#	if Mcorr is not None:
#		M = M_ + Mcorr
#	else:
#		M = M_.copy()
#	dist = 10. ** (0.2 * (lim - M + 5.)) # parsecs
#	dist /= 1e6
#	bad_M = (M < -30.) | (M > -10.) # conservative removal of non-galaxy objects
#	dist = np.where(bad_M, -99., dist)
#	try:
#		dist = np.where(dist > ceiling, ceiling, dist)
#	except:
#		dist = np.where(dist > ceiling.value, ceiling.value, dist)
#	return dist
#
#def get_zmax_hdu(cat, mcol, mlim, zcol, ceiling=None, kcorrs=None):
#	print '\t"%s" at limit: %s' % (mcol, mlim)
#
#	# calculate maximum distance to each galaxy
#	mag = cat[mcol]
#	if kcorrs is None:
#		kcorr = np.zeros_like(mag)
#	else:
#		kcorr = cat[kcorrs]
#		print '\t"%s" correction: mean mag %.2f -> %.2f' % (kcorrs, mag.mean(), (mag-kcorr).mean())
#	max_distance = dmax(mlim, mag, Mcorr=-kcorr) #, Mcorr=-5*log10(0.7)) for a cosmolgy-correction
#	max_redshift = np.zeros_like(max_distance)
#	mask = max_distance > 0
#
#	# interpolate these onto the maximum redshift at which the object would be observed
#	zmin, zmax = 0., 6.
#	z_grid = np.linspace(zmin, zmax, 100)
#	distance_grid = cmd(z_grid)
#	max_redshift[mask] = interp1d_(max_distance[mask], distance_grid, z_grid)
#	max_redshift[~mask] = -99.
#	max_distance[~mask] = -99.
#	zmax_colname = mcol+'_fl%.1f_zmax'%mlim
#	dmax_colname = mcol+'_fl%.1f_dmax'%mlim
#	max_cols = [fitscol(array=max_redshift, format='D', name=zmax_colname),
#				fitscol(array=max_distance, format='D', name=dmax_colname)]
#	# create fits table from new+old columns
#	for col in max_cols:
#		if col.name in cat.columns.names:
#			cat.columns.del_col(col.name)
#	max_cols = fits.ColDefs(max_cols)
#	hdu = fits.BinTableHDU.from_columns(max_cols + cat.columns)
#	return hdu

	# draw from fit to redshift distribution of real galaxies
	#if zg is not None and Pz is not None:
		#pass
		#for i in tqdm(range(Nrand), desc='\t\tclone batches', ascii=True, ncols=100):
		#	zdraw = pchoice(zg, p=Pz, size=Ngal)
		#	mask = maxcol > 0
		#	nbad = (zdraw[mask] > maxcol[mask]).sum() # number of long-draws
		#	dnbad = nbad*1
		#	while nbad > 0:
		#		#print nbad
		#		badz = np.where((zdraw > maxcol) & mask)[0]
		#		if dnbad <= 10: # if re-draws are slowing down, loop over remaining draws with truncated P(z)
		#			for bz in badz:
		#				zg1, Pz1 = zg[zg < maxcol[bz]], Pz[zg < maxcol[bz]]
		#				Pz1 /= Pz1.sum()
		#				zdraw[bz] = pchoice(zg1, p=Pz1, size=1)[0]
		#		else: # otherwise keep re-drawing
		#			zdraw[badz] = pchoice(zg, p=Pz, size=len(badz))
		#		dnbad = nbad - (zdraw[mask] > maxcol[mask]).sum()
		#		nbad = (zdraw[mask] > maxcol[mask]).sum()
		#	zdraw[~mask] = -99.
		#	id_zdraw = np.column_stack((idcol.copy(), zdraw))
		#	randoms_id_z[i*Ngal:(i+1)*Ngal] = id_zdraw
		#mask = maxcol > 0
		#mask = np.array(list(mask) * Nrand)
		#randoms_id_z[:, 1] = np.where(mask, randoms_id_z[:, 1], -99.)
		#return randoms_id_z

#def get_zmax(z, obsmag, kcorr, fluxlimit):
#	distmod = cosmo.distmod(z).value
#
#	absmag = obsmag - distmod - kcorr
#	max_distmod = fluxlimit - absmag
#	max_distance = 10. ** (max_distmod/5. - 5.) # Mpc
#
#	zgrid = np.linspace(z.min(), z.max(), 200)
#	dgrid = cmd(zgrid).value
#	max_redshift = interp1d_(max_distance, dgrid, zgrid)
#	max_redshift = np.where(np.isnan(max_redshift) | np.isinf(max_redshift), -99., max_redshift)
#	return max_redshift
#
#def find_M_crossing(M_, fluxlimit, kcorrs, Mgrid=0):
#	print '\tk-correcting for M(z)..'
#	# define analytical absolute flux-limit
#	Mlim_z = lambda z: fluxlimit - 5.*log10(cmd(z).value * 1e6 / 10.)
#
#	# apply k-corrections to get M(z) per galaxy
#	zgrid = np.unique(kcorrs['z'])
#	kmatrix = np.zeros([len(zgrid), len(M_)])
#	for i, zpoint in enumerate(zgrid):
#		maggy_ratio = kcorrs['maggy_ratio'][np.where(kcorrs['z'] == zpoint)[0]].values
#		kcorrection = -2.5 * log10(maggy_ratio)
#		kmatrix[i] += kcorrection
#	Mmatrix = kmatrix + M_
#
#	# interpolate onto a finer redshift grid and find the
#	# redshift intersection of M(z) and the flux-limit
#	print '\tfinding flux-limit intersection..'
#	fine_zgrid = np.linspace(0.001, zgrid[-1], 10*len(zgrid))
#	fine_Mmatrix = interp1d(zgrid, Mmatrix, axis=0)(fine_zgrid)
#	Mlim = Mlim_z(fine_zgrid)
#	Mz_minus_Mlim = (fine_Mmatrix.T - Mlim).T
#	minima = fine_zgrid[np.argmin(abs(Mz_minus_Mlim), axis=0)]
#	mask = (M_ < -26) | (M_ > -15) # remove intrinsically very faint/bright objects
#	minima[mask] = -99.
#
#	print '\tdone.'
#	if Mgrid:
#		return fine_Mmatrix, fine_zgrid
#	else:
#		return minima
