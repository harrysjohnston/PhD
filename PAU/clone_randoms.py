from tqdm import tqdm
from numpy import log10
from scipy import stats
from astropy.io import fits
from astropy.table import Table
from astropy.cosmology import FlatLambdaCDM as FLCDM
from scipy.interpolate import interp1d, interp2d, InterpolatedUnivariateSpline
import os
import h5py
import argparse
import warnings
import numpy as np
import pandas as pd
import astropy.units as u
import scipy.optimize as so
fitscol = fits.Column
pchoice = np.random.choice
cosmo = FLCDM(Om0=0.25, H0=70, Ob0=0.044)
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
norm = stats.norm(loc=0, scale=1)
rnorm_prep = norm.rvs(size=10000)
minmax = lambda x: [x.min(), x.max()]
frac = lambda x: np.sum(x)/np.float32(len(x))
nn = lambda x: x[~np.isnan(x) & ~np.isinf(x)]

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

def get_k_z(kcorrs, cat_ids):
	"""
	Read maggy ratios from file, where each can be converted
	to a k-correction from z (given by column) to z=0.
	Return:
		1. k-correction matrix: shape = (N galaxies x len(z_grid))
		2. z_grid: at which the k(z->0) is to computed
		3. mask: indicates bad photometry
	"""
	df_kcorrs = h5py.File(kcorrs, 'r')
	ID = df_kcorrs['ID'][:]
	z = df_kcorrs['z'][:]
	maggyratio = df_kcorrs['maggy_ratio'][:]
	print '\tcutting IDs..'
	idcut = np.isin(ID, cat_ids)
	ID = ID[idcut]
	z = z[idcut]
	maggyratio = maggyratio[idcut]
	idsort = np.argsort(ID[:len(np.unique(ID))])

	zero_ratio = maggyratio[np.where(z == 0)[0]] # mgy(z=0) / mgy(z=0)
	mask = zero_ratio != 1 # ratio != 1 indicates error in computation of kcorrections -- should trace back to bad photometry

	z_grid = np.unique(z)
	uID = np.unique(ID)
	maggyratio_matrix = np.array(maggyratio).reshape(len(z_grid), len(uID)).T
	maggyratio_matrix[mask] = 1e-1
	k_z = -2.5 * log10(maggyratio_matrix)
	return k_z[idsort], z_grid, mask[idsort]

def interp_k(z_grid, k_z, z_out):
	print '\tinterpolating k-corrections..'
	x = np.array([interp1d(z_grid, k_z[i], #axis=1,
				assume_sorted=True, copy=False,
				bounds_error=False, fill_value=z_grid.max())(z_out) for i in range(len(k_z))])
	return x

def interp_dm(z_grid, d_grid, z_out):
	x = interp1d(z_grid, d_grid,
				assume_sorted=True, copy=False,
				bounds_error=False, fill_value=z_grid.max())(z_out)
	return x

def fit_zmax(fluxlim, m_obs, z_obs, z_grid, k_z, min=0, Q=0):
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

		loss_fn = LHS[idx] - (dm_max + k_max + Q * z_max)
		loss_fn[np.isnan(loss_fn) | np.isinf(loss_fn)] = 0.
		return np.abs(loss_fn)

	if not min:
		print '\tfitting z_max per galaxy..'
	else:
		print '\tfitting z_min per galaxy..'
	x = minimize(distmod_relation, z_obs, bounds=[0., z_grid.max()])
	return x

def minimize(fun, x0, tol=1e-1, bounds=[-np.inf, np.inf], quiet=False):
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
		#if all(np.array(fins[-50:]) == 0):
		if np.array(fins[-50:]).mean() < 10:
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

def get_tgauss_window(d, hilim, lolim, area=180., volume=3.5e6, V_base=None, zmax=0.6):
	Om = sqdeg2ster(area)
	dmax = get_d(zmax)
	V_mid = midpoints(V_base)
	window_fn = np.zeros_like(V_mid)
	# volume coordinates
	hilim_vol = (Om/3.) * hilim**3.
	lolim_vol = (Om/3.) * lolim**3.
	vol_at_gal = (Om/3.) * d**3.
	Vmax = (Om/3.) * (hilim**3. - lolim**3.)

	if d > dmax or hilim == -99. or lolim == -99. or lolim == hilim or hilim < 10.:
		# throw away those with poor k-corrections->limits
		window_fn *= np.nan
	else:
	#	if Vmax < 4.*volume:
	#		# reduce window size for galaxies with small Vmax
	#		volume = Vmax / 4.
		# draw from a truncated [-2, 2] Gaussian (pre-drawn); convert to volumes
		rvols = rtnorm_prep * volume # units of sigma
		vol_draw = vol_at_gal + rvols # volume at galaxy + some deviation
		vol_draw = np.abs(vol_draw)
		# reflect in boundaries, staying in volume units
		x = 0
		while any(vol_draw > hilim_vol) or any(vol_draw < lolim_vol):
			vol_draw[vol_draw < lolim_vol] += 2.*(lolim_vol - vol_draw[vol_draw < lolim_vol])
			vol_draw[vol_draw > hilim_vol] += 2.*(hilim_vol - vol_draw[vol_draw > hilim_vol])
			vol_draw = np.abs(vol_draw)
			x += 1
			if x > 200:
				#print 'flattening window for zmax = %.4f galaxy' % get_z(hilim)
				window_fn = np.ones_like(V_mid, dtype=float)
				window_fn[(V_mid < lolim_vol) | (V_mid > hilim_vol)] = 0.
				return np.stack((V_mid, window_fn))
				#break

		# bin, smooth and normalise for the pdf
		vol_minmax = minmax(vol_draw)
		window_fn, bin_edges = np.histogram(vol_draw, bins='auto')
		window_fn = np.interp(V_mid, midpoints(bin_edges), window_fn)
		window_fn[(V_mid < vol_minmax[0]) | (V_mid > vol_minmax[1])] = 0.

	# normalisation must be such that int[W(V) dV] = 1
	#window_fn = np.asarray(window_fn, dtype=np.float32) / np.trapz(window_fn, x=V_mid)
	window_fn = np.asarray(window_fn, dtype=np.float32)
	window_fn /= 1.*window_fn.sum()
	window_fn = np.stack((V_mid, window_fn))
	return window_fn

def clone_galaxies(idcol, maxcol, Nrand=10, zlims=None, window_vol=None, area=180., dspec=None, save_diag=False,
				   dobs=None, drawn_dobs=None, dres=1., Niter=15, load_windows=True, runid='', mincol=None):
	# take arrays of IDs and max distances/redshifts
	# draw Nrand redshifts from quadratic/fitted distribution [0, {d/z}max] per galaxy
	windowed = window_vol is not None
	Om = sqdeg2ster(area)
	Ngal = len(idcol)
	zmin, zmax = zlims or (0., 0.6)
	dmin, dmax = cmd([zmin, zmax])
	Vminimum = (Om/3.) * dmin**3.
	Vmaximum = (Om/3.) * dmax**3.
	maxcol[maxcol < dmin] = -99.
	if dspec:
		zph_bins = np.arange(zmin, zmax+0.1, step=0.1)
		zspec = get_z(dspec)
	if drawn_dobs is None:
		drawn_dobs = dobs
	if dspec is not None:
		print '\t\t(', (drawn_dobs > maxcol).sum(), 'objects observed beyond calculated maximum -- REMOVING )'
		maxcol[drawn_dobs > maxcol] = -99.#drawn_dobs[drawn_dobs > maxcol]
	if mincol is not None:
		print '\t\t(', (drawn_dobs < mincol).sum(), 'objects observed closer than bright limit -- REMOVING )'
		mincol[drawn_dobs < mincol] = -99.#drawn_dobs[drawn_dobs < mincol]
	else:
		mincol = dmin * np.ones_like(maxcol)

	# get limits; zmax/zmin or survey edges
	upps = np.stack((maxcol, np.ones_like(maxcol)*dmax))
	lows = np.stack((np.zeros_like(maxcol), mincol))
	hilim = np.min(upps, axis=0)
	lolim = np.max(lows, axis=0)
	hilim_vol = (Om/3.) * hilim**3.
	lolim_vol = (Om/3.) * lolim**3.

	# ready comoving distance & volume baselines
	d_base = np.arange(0, dmax, step=dres)
	V_base = (Om/3.) * d_base**3.
	V_mid = midpoints(V_base)
	d_mid = (3.*V_mid/Om)**(1./3.)
	Vobs = (Om/3.) * dobs**3.

	# take galaxy n(chi)
	n_g = np.asarray(np.histogram(Vobs, bins=V_base)[0], dtype=np.float32)

	if windowed:
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
			windows = [get_tgauss_window(drawn_dobs[i], hilim[i], lolim[i], volume=1.*window_vol, V_base=V_base, zmax=zmax, area=area)
						for i in tqdm(range(len(drawn_dobs)), desc='\t\twindowing', ncols=100)]
			windows = np.array(windows)[:, 1, :]
			assert (not all(np.isnan(windows.flatten()))), "Windowing failed! All windows are NaN; d > dmax or bad limits"

		#with h5py.File('windows%s.h5'%runid, 'w') as f:
		#	f.create_dataset('windows', data=windows)
		#	f.close()

	# windows in volume coords
	flat_windows = np.ones([len(maxcol), len(V_mid)], dtype=np.float32)
	for i in tqdm(range(len(flat_windows)), desc='\t\tlimiting flat windows', ncols=100):
		flat_windows[i][(V_mid < lolim_vol[i]) | (V_mid > hilim_vol[i])] = 0.

	if not windowed:
		windows = flat_windows

	# get Vmax
	Vmax = Om * np.trapz(d_mid**2. * windows, x=d_mid, axis=1)
	flat_Vmax = np.trapz(d_mid**2. * flat_windows, x=d_mid, axis=1)

	# setup diagnostic save-outs
	Vmax_dc_list = []
	n_clones_list = []
	Delta_d_list = []
	n_r_list = []
	
	this_iter = 1
	print '\t\titerating n_clones calculation..'
	while this_iter < Niter+1:
		print '\t\t======= ITERATION #%s =======' % this_iter
		# get Vmax,dc ; density-weighted Vmax for n_clone calculation
		if this_iter == 1:
			Delta_d = np.ones_like(V_mid, dtype=np.float64)
		else:
			n_r = np.asarray(np.histogram((Om/3.)*ddraw**3., bins=V_base)[0], dtype=np.float32)
			Delta_d = Nrand * n_g / n_r
		delta_mask = ~np.isnan(Delta_d) & ~np.isinf(Delta_d)
		print '\t\tfrac(delta==nan/inf) = %.2f'%frac(~delta_mask)
		Vmax_dc = Om * np.trapz( d_mid[delta_mask]**2. \
								* Delta_d[delta_mask] \
								* windows.T[delta_mask].T, \
								x=d_mid[delta_mask], axis=1)

		bad_Vmaxdc = np.isnan(Vmax_dc) | np.isinf(Vmax_dc) | (Vmax_dc <= 0)
		if this_iter > 1:
			prev_N = n_clones.sum()
			n_clones = np.asarray(np.round(Nrand * Vmax / Vmax_dc), dtype=int)
		else:
			n_clones = np.repeat(Nrand, len(Vmax))
		n_clones[bad_Vmaxdc] = 0
		n_clones[n_clones < 0] = Nrand
		print '\t\ttotal N clones = %s'%n_clones.sum()
		if this_iter > 1:
			print '\t\t\t\ti.e. %+d'%(n_clones.sum() - prev_N)

		ddraw, clone_ids = [], []
		for i in tqdm(range(len(n_clones)), desc='\t\tcloning', ncols=100):
			# skip bad galaxies
			if n_clones[i] == 0:
				continue

			if not windowed or this_iter == 1:
				# draw uniformly from the allowed volume
				Vdraw = np.random.rand(n_clones[i]) * (hilim_vol[i] - lolim_vol[i]) + lolim_vol[i]
			else:
				# over-draw uniformly and reject excess clones according to window function
				Vdraw1 = np.random.rand(n_clones[i]*100) * (hilim_vol[i] - lolim_vol[i]) + lolim_vol[i]
				window_at_Vdraw = np.interp(Vdraw1, V_mid, windows[i])
				if window_at_Vdraw.sum() == 0:
					# skip bad windows
					print 'galaxy', idcol[i], 'at z=%.3f'%get_z(dobs[i]), ' -- window is too small?'
					continue
				window_at_Vdraw /= window_at_Vdraw.sum()
				Vdraw = np.random.choice(Vdraw1, p=window_at_Vdraw, size=n_clones[i])

			ddraw_i = (3.*Vdraw/Om)**(1./3.)
			assert all(ddraw_i <= hilim[i]) and all(ddraw_i >= lolim[i]), "z-limiting of clones is broken!"
			for ddi in ddraw_i:
				ddraw.append(ddi)
				clone_ids.append(idcol[i])

		ddraw = np.array(ddraw).flatten()
		clone_ids = np.array(clone_ids).flatten()
		assert len(ddraw) == len(clone_ids), "cloning going wrong!"

		Delta_d_list.append(Delta_d)
		Vmax_dc_list.append(Vmax_dc)
		n_clones_list.append(n_clones)	
		try: n_r_list.append(n_r)
		except: pass
		this_iter += 1

	if save_diag:
		with h5py.File('diagnostics_clonerandoms%s.h5'%runid, 'w') as h5_diag:
			h5_diag.create_dataset('Vmax', data=Vmax)
			h5_diag.create_dataset('Vmax_dc', data=Vmax_dc_list)
			h5_diag.create_dataset('V_mid', data=V_mid)
			h5_diag.create_dataset('d_mid', data=d_mid)
			h5_diag.create_dataset('Delta_d', data=Delta_d_list)
			h5_diag.create_dataset('n_clones', data=n_clones_list)
			h5_diag.create_dataset('n_r', data=n_r_list)
			h5_diag.create_dataset('Om', data=Om)
			h5_diag.close()

	mask = ddraw != -99.
	zdraw = np.ones_like(ddraw) * -99.
	zdraw[mask] = get_z(np.abs(ddraw[mask]))

	randoms_id_z = np.column_stack((clone_ids, zdraw))
	return randoms_id_z

def main(args):
	print '\n'
	for cat_path in args.catalogues:
		print 'reading %s..' % cat_path
		cat = fits.open(cat_path)[1].data
		cat = cat[np.argsort(cat[args.idcol])]
		t = Table(cat)

		if args.refresh_zmax and not args.zmax_col and not args.zph_max:
			# establish zmax per detection band
			for mcol, mlim, kcorr, blim in zip(args.magcols, args.maglims, args.kcorrs, args.brightlim):
				print '\t"%s" at limit: %s' % (mcol, mlim)
				print '\tbright limit: %s' % blim
				print '\treading "%s" k-corrections..' % kcorr

				k_z, z_grid, mask = get_k_z(kcorr, cat[args.idcol])
				max_redshift = fit_zmax(mlim, cat[mcol], cat[args.zcol], z_grid, k_z, Q=args.Q)
				max_redshift[mask] = -99.

				zmax_col = mcol+'_fl%.1f_zmax'%mlim
				t[zmax_col] = max_redshift

				if blim is not None:
					zmin_col = mcol+'_fl%.1f_zmin'%blim
					min_redshift = fit_zmax(blim, cat[mcol], cat[args.zcol], z_grid, k_z, min=1, Q=args.Q)
					t[zmin_col] = min_redshift
				else:
					min_redshift = np.zeros(len(t))

				Vmax = (sqdeg2ster(args.area)/3.) * (get_d(max_redshift)**3. - get_d(min_redshift)**3.)
				Vmax[mask] = -99.
				t['Vmax'] = Vmax

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
			#if args.zweight is None:
			#	wcol = None
			#	bins = 'auto'
			#else:
			#	wcol = cat[args.zweight]
			#	bins = np.linspace(zcol[mask].min(), zcol[mask].max(), 40)

			# construct redshift distribution per band
			#print '\t\tcloning..'
			random_cols = []
			for mcol, mlim, blim in zip(args.magcols, args.maglims, args.brightlim):
				if not args.zmax_col:
					zmax_col = mcol+'_fl%.1f_zmax'%mlim
				else:
					zmax_col = args.zmax_col
				if blim is not None:
					zmin_col = mcol+'_fl%.1f_zmin'%blim
					mincol = get_d(cat[zmin_col])
				else:
					mincol = None

				dobs = get_d(zcol)
				if not args.zph_max:
					maxcol = get_d(cat[zmax_col])
					drawn_dobs = None
				else:
					maxcol = get_d(args.zmax_draw)
					drawn_dobs = dobs = get_d(args.zspec_draw)

				if len(maxcol) != len(dobs):
					idcut1 = np.isin(idcol, zmax_id)
					idcut2 = np.isin(zmax_id, idcol)
					idcol, dobs = idcol[idcut1], dobs[idcut1]
					if mincol: mincol = mincol[idcut1]
					maxcol = maxcol[idcut2]

				randoms_id_z = clone_galaxies(idcol, maxcol, args.Nrand, zlims=args.zlims, window_vol=args.window, area=args.area, dres=args.dres, save_diag=args.save_diag,
											  dobs=dobs, drawn_dobs=drawn_dobs, Niter=args.niter, load_windows=args.load_windows, runid=args.id, mincol=mincol)
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
			dec = np.random.choice(dec, p=cos, size=len(randoms_id_z)) + \
					(np.random.rand(len(randoms_id_z))-0.5) * (1./60.)

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
	return out

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
	parser.add_argument(
		'-zmax_col',
		type=str,
		help='column name for zmax -- overrides fitting of zmax given k-corrections')
	parser.add_argument(
		'-Q',
		type=float,
		default=0.,
		help='include a scalar evolution term (float), default=0')
	parser.add_argument(
		'-zph_max',
		type=str,
		nargs=2,
		help='2 args: (1) path to .zmaxtable file, (2) number of zmax vectors to sample')
	parser.add_argument(
		'-save_diag',
		type=int,
		default=0,
		help='(1) save diagnostics or (0) not (default)')
	args = parser.parse_args()

	#if args.kcorrs is None:
	#	args.kcorrs = [None]*len(args.magcols)
	if args.id != '':
		args.id = '_'+args.id
	if args.brightlim is None:
		args.brightlim = [None]*len(args.maglims)

	if args.zph_max:
		import os
		import multiprocessing as mp
		from glob import glob
		zmaxtable_h5 = h5py.File(args.zph_max[0], 'r')
		zmaxtable = zmaxtable_h5['zmax'][:]
		zspectable = zmaxtable_h5['zspec'][:]
		if zmaxtable.shape[0] > zmaxtable.shape[1]:
			zmaxtable = zmaxtable.T
		if zspectable.shape[0] > zspectable.shape[1]:
			zspectable = zspectable.T
		zmax_id = zmaxtable_h5['ID'][:]
		zmaxtable = zmaxtable[:, np.argsort(zmax_id)]
		zspectable = zspectable[:, np.argsort(zmax_id)]
		Ndraws = int(args.zph_max[1])
		indices = np.random.choice(len(zmaxtable), replace=False, size=Ndraws)
		orig_id = args.id

		def ensemble_randoms(x):
			setattr(args, 'zmax_draw', zmaxtable[x])
			setattr(args, 'zspec_draw', zspectable[x])
			args.id = orig_id + '_%s'%str(x).zfill(2)
			outfile_x = main(args)
			return outfile_x

		if mp.cpu_count() > 1:
			pool = mp.Pool(mp.cpu_count())
			draws = pool.map(ensemble_randoms, indices)
			pool.close()
		else:
			draws = []
			for x in indices:
				draws.append(ensemble_randoms(x))

		cats = []
		for fil in draws:
			cats.append(fits.open(fil)[1].data)

		mcat = Table(np.concatenate(cats))
		if args.o is None:
			out = args.catalogues[0].replace('.fits', '_CloneZIDRandoms.fits')
		else:
			out = args.o
		if args.id is not None:
			out = out.replace('.fits', '%s.fits'%args.id)

		mcat.write(out, format='fits', overwrite=1)
		for fil in draws:
			os.system('rm '+fil)
	else:
		main(args)



