from tqdm import tqdm
from numpy import log10
from scipy import stats
from astropy.io import fits
from astropy.table import Table
from astropy.cosmology import FlatLambdaCDM
from scipy.interpolate import interp1d, interp2d, InterpolatedUnivariateSpline
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

def get_k_z(kcorrs):
	"""
	Read maggy ratios from file, where each can be converted
	to a k-correction from z (given by column) to z=0.
	Return:
		1. k-correction matrix: shape = (N galaxies x len(z_grid))
		2. z_grid: at which the k(z->0) is to computed
		3. mask: indicates bad photometry
	"""
	df_kcorrs = pd.read_csv(kcorrs, delimiter=' ')
	ID = df_kcorrs['ID']
	idsort = np.argsort(ID[:len(set(ID))])
	z = df_kcorrs['z']
	maggyratio = df_kcorrs['maggy_ratio']
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

def fit_zmax(fluxlim, m_obs, z_obs, z_grid, k_z):
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
	def distmod_relation(z_max, idx, kcorrect=1):
		# LHS = mu(z_max) + k(z_max)
		# the array z_max has length = # of galaxies
		# pass idx to specify indices of galaxies within main array

		dm_max = interp_dm(z_grid_fine, d_grid_fine, z_max)
		if kcorrect:
			k_max = np.array([ki(z) for ki, z in np.column_stack((k_interpolators[idx], z_max))])
		else:
			k_max = np.zeros_like(z_max)

		loss_fn = LHS[idx] - (dm_max + k_max)
		loss_fn[np.isnan(loss_fn) | np.isinf(loss_fn)] = 0.
		return abs(loss_fn)

	print '\tfitting z_max per galaxy..'
	x = minimize(distmod_relation, z_obs, bounds=[0., z_grid.max()])
	return x

def minimize(fun, x0, tol=1e-2, bounds=[-np.inf, np.inf]):
	# find variables x that bring fun(x) down to < tol[mags]
	# starting with guess x0
	# perturb x0 by +/- 10% -> x1, and evaluate fun(x1)
	# if fun(x1) < fun(x0), perturb x1 by 10% again,
	# but now with a smaller chance to flip the sign
	# if fun(x1) > fun(x0), perturb x1 by 10% again,
	# but now with a larger chance to flip the sign
	nflip = np.array([1., -1., -1., -1. , -1.])
	pflip = np.array([1., 1., 1., 1., -1.]) # 20% chances to go the 'wrong' way
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
	kcorrect = False
	fx = fun(x, idx, kcorrect=kcorrect)
	fx0 = fun(x0, idx, kcorrect=kcorrect)

	fins = []
	while any(x_out == 0):
		fin = abs(fx) < tol
		fins.append(fin.sum())
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
		if all(np.array(fins)[-5:] == 0):
			x[x < x0[idx1]] = x0[idx1][x < x0[idx1]]
			s[x < x0[idx1]] = 1.
		fx = fun(x, idx1, kcorrect=kcorrect)

		try:
			for i in range(8):
				print 'galaxy %s: delta=%.3f[mag], last z_max=%.3f, z of object=%.3f' % (i+1, fx[i], x[i], x0[idx1][i])
			print 'N outside tolerance =', (~fin).sum()
		except:
			pass
		#if all(np.array(fins)[-200:] == 0):
		#	np.savetxt('badgalaxies.txt', np.column_stack((idx1, x)))
		#	print 'saving failures to ./badgalaxies.txt !'
		#	x_out[idx1] = x
		#	break

	return x_out

def get_windows(z, area=180., volume=3.5e6):
	area = area / (180. / np.pi)**2. # convert sqdeg to steradians
	def invvol(r1, vol):
		r2 = ((3. * vol / area) + r1**3.) ** (1./3.)
		# calculate asymmetric limits of volume -- see notes
		with warnings.catch_warnings():
			warnings.simplefilter("ignore")
			new_r1 = ((3.*r1**3. - r2**3.) / 2.) ** (1./3.)
		new_r2 = ((r1**3. + r2**3.) / 2.) ** (1./3.) # 1s limits
		# fix NaNs in new_r1 -- make positive for the cube-root, and then negative again
		badr1 = np.isnan(new_r1)
		new_r1[badr1] = np.abs(((3.*r1[badr1]**3. - r2[badr1]**3.) / 2.)) ** (1./3.)
		new_r1[badr1] *= -1.
#		new_r1 = r1 - 2*(r1 - new_r1)
#		new_r2 = r1 + 2*(new_r2 - r1) # 2s limits
		return [new_r1, r1, new_r2]

	d = get_d(z)
	widths = np.asarray(invvol(d, volume)).squeeze()
	np.savetxt('./asymmlimits_2s_vol.txt', widths)
	return widths

def clone_galaxies(idcol, maxcol, Nrand=10, zlims=None, windows=None, dobs=None):
	# take arrays of IDs and max distances/redshifts
	# draw Nrand redshifts from quadratic/fitted distribution [0, {d/z}max] per galaxy
	Ngal = len(idcol)
	zmin, zmax = zlims or (0., 6.)
	dmin, dmax = cmd([zmin, zmax])

	Nrand_windows = np.array(list(windows.T) * Nrand)
	Nrand_maxcol = np.array(list(maxcol) * Nrand)
	Nrand_idcol = np.array(list(idcol) * Nrand)
	Nrand_dobs = np.array(list(dobs) * Nrand)
	print '\t\t(', (dobs > maxcol).sum(), 'objects observed beyond calculated maximum -- setting observed=max )'
	Nrand_maxcol[Nrand_dobs > Nrand_maxcol] = Nrand_dobs[Nrand_dobs > Nrand_maxcol]
	badmax = (Nrand_maxcol == -99.) | (Nrand_maxcol < 60.)

	if windows is not None: # if windowing, 'maxcol' must be in comoving Mpc,
						   # and must also pass observed distances 'dobs'
						   # 'windows' should be the 2sigma comoving Gaussian width PER OBJECT
		# galaxy-density weighting - may not be needed
		#ddraw = np.ones(Ngal * Nrand) * Nrand_dobs
		#dbins = np.linspace(dmin, dmax, 51)
		#dmids = (dbins[1:] + dbins[:-1]) / 2.
		#dVdz = dmids ** 2.
		#rhist = np.histogram(ddraw, bins=dbins)
		#ghist = np.histogram(dobs, bins=dbins)
		#ovd_z = np.ones(len(dmids))
		#limits = np.array(list(dbins[1:]) * len(maxcol))
		#limits = np.where(# this needs to apply a top-hat to the integral, with a zmax per object
		#ovd_max = np.trapz(ovd_z * dVdz * (dbins[1]-dbins[0]))

		#draw1 = (np.random.rand(len(scales)) - 0.5) * 4. * scales
		# windows gives -1s, d_obs, +1s limits
		lohi_sigmas = np.column_stack((Nrand_windows[:, 1] - Nrand_windows[:, 0],
									   Nrand_windows[:, 2] - Nrand_windows[:, 1]))
		tnorm = stats.truncnorm(-2, 2, loc=0, scale=1)
		rtnorm = tnorm.rvs(size=len(lohi_sigmas))
		draw1 = np.where(rtnorm < 0, -rtnorm*lohi_sigmas[:, 0], rtnorm*lohi_sigmas[:, 1])

		ddraw = Nrand_dobs + draw1
		lowbound = ~badmax & (ddraw < 0)
		highbound = ~badmax & (ddraw > dmax)
		ddraw[lowbound] = -ddraw[lowbound] # reflection in z=0
		ddraw[highbound] = 2.*dmax - ddraw[highbound] # reflection in z-limit
		ddraw[badmax] = -99.

		oob = (ddraw > Nrand_maxcol) & ~badmax
		# re-draw if outside of bounds
		while oob.sum() != 0:
			print 1.*oob.sum()/len(oob)
			#draw1 = (np.random.rand(len(scales[oob])) - 0.5) * 4. * scales[oob]

			rtnorm = tnorm.rvs(size=len(lohi_sigmas[oob]))
			draw1 = np.where(rtnorm < 0, -rtnorm*lohi_sigmas[oob][:, 0], rtnorm*lohi_sigmas[oob][:, 1])
			ddraw[oob] = Nrand_dobs[oob] + draw1
			lowbound = ~badmax & (ddraw < 0)
			highbound = ~badmax & (ddraw > dmax)
			ddraw[lowbound] = -ddraw[lowbound] # reflection in z=0
			ddraw[highbound] = 2.*dmax - ddraw[highbound] # reflection in z-limit
			oob = (ddraw > Nrand_maxcol) & ~badmax
	else:
		# draw from d^2 (== growth of volume element) between [0, d_max]
		ddraw = np.random.power(3, Ngal * Nrand) * Nrand_maxcol * u.Mpc
		ddraw[badmax] = -99.
		ddraw = ddraw.value

	mask = ddraw > 0
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

		if args.window is None:
			# establish zmax per detection band
			for mcol, mlim, kcorr in zip(args.magcols, args.maglims, args.kcorrs):
				print '\t"%s" at limit: %s' % (mcol, mlim)
				print '\treading "%s" k-corrections..' % kcorr

				k_z, z_grid, mask = get_k_z(kcorr)
				max_redshift = fit_zmax(mlim, cat[mcol], cat[args.zcol], z_grid, k_z)
				max_redshift[mask] = -99.

				zmax_col = mcol+'_fl%.1f_zmax'%mlim
				t[zmax_col] = max_redshift

			t.write(cat_path, format='fits', overwrite=1)
		cat = t.copy()

		if args.randoms:
			print '\tbuilding clone randoms..'
			# drop points at random in 1x1 square
			Nrand = len(cat) * args.Nrand
			ra = np.random.rand(Nrand)
			dec = np.random.rand(Nrand)
			random_cols = [fitscol(array=ra, name='ra', format='D'),
						   fitscol(array=dec, name='dec', format='D')]
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
			print '\t\tcloning..'
			for mcol, mlim in zip(args.magcols, args.maglims):
				zmax_col = mcol+'_fl%.1f_zmax'%mlim
				#zmax_col = 'zmax_19p8'

				if args.window:
					windows = get_windows(zcol)
				else:
					windows = None					

				maxcol = get_d(cat[zmax_col])
				dobs = get_d(zcol)
				randoms_id_z = clone_galaxies(idcol, maxcol, args.Nrand, zlims=args.zlims, windows=windows, dobs=dobs)
				randoms_comoving = get_d(randoms_id_z[:, 1])

				random_cols.append(fitscol(array=randoms_id_z[:, 0], name=mcol+'_cloneID', format='K'))
				random_cols.append(fitscol(array=randoms_id_z[:, 1], name=mcol+'_cloneZ', format='D'))
				random_cols.append(fitscol(array=randoms_comoving, name=mcol+'_cloneComovingDist', format='D'))

			# save randoms catalogue with galaxy IDs and redshifts for selection
			random_cols = fits.ColDefs(random_cols)
			rand_hdu = fits.BinTableHDU.from_columns(random_cols)
			if args.o is None:
				out = cat_path.replace('.fits', '_CloneZIDRandoms.fits')
			else:
				out = args.o
			rand_hdu.writeto(out, overwrite=1)

	print 'done!'

if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument(
		'catalogues',
		nargs='*',
		default='give paths to fits catalogue(s) with redshifts and absolute rest-frame magnitudes')
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
		'-kcorrs',
		nargs='*',
		type=str,
		help='specify k-corrections filename(s), for computation of zmax(s). File must have columns ID, z, maggie_ratio')
	parser.add_argument(
		'-randoms',
		type=int,
		help='1 = create a randoms catalogue with 1 redshift distribution per -magcol')
	parser.add_argument(
		'-window',
		help='1 = scatter clones within a truncated Gaussian window, per-object')
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
