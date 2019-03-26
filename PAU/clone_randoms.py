from astropy.io import fits
from astropy.table import Table
from astropy.cosmology import FlatLambdaCDM, z_at_value
from functions import fit_smail
from tqdm import tqdm
from numpy import log10
from scipy.interpolate import interp1d
import astropy.units as u
import numpy as np
import argparse
fitscol = fits.Column
pchoice = np.random.choice
cosmo = FlatLambdaCDM(H0=70, Om0=0.3)
interp1d_ = lambda x, xx, xy: interp1d(xx, xy, bounds_error=0, fill_value=-99.)(x)

def dmax(lim, M_, ceiling=np.inf, hcorr=None):
	# maximum distance of absolute M galaxy in Mpc, given apparent mag lim
	if hcorr is not None:
		M = M_ + hcorr
	else:
		M = M_.copy()
	dist = 10. ** (0.2 * (lim - M + 5.)) # parsecs
	dist /= 1e6
	bad_M = (M < -30.) | (M > -10.) # conservative removal of non-galaxy objects
	dist = np.where(bad_M, -99., dist)
	try:
		dist = np.where(dist > ceiling, ceiling, dist)
	except:
		dist = np.where(dist > ceiling.value, ceiling.value, dist)
	return dist

def get_zmax_hdu(cat, mcol, mlim, zcol, ceiling=None):
	print '\t"%s" at limit: %s' % (mcol, mlim)
	# calculate maximum distance to each galaxy
	mag = cat[mcol]
	max_distance = dmax(mlim, mag)#, hcorr=-5*log10(0.7))
	max_redshift = np.zeros_like(max_distance)
	mask = max_distance > 0
	# interpolate these onto the maximum redshift at which the object would be observed
	zmin, zmax = 0., 6.
	z_grid = np.linspace(zmin, zmax, 100)
	distance_grid = cosmo.comoving_distance(z_grid)
	max_redshift[mask] = interp1d_(max_distance[mask], distance_grid, z_grid)
	max_redshift[~mask] = -99.
	max_distance[~mask] = -99.
	zmax_colname = mcol+'_fl%.1f_zmax'%mlim
	dmax_colname = mcol+'_fl%.1f_dmax'%mlim
	max_cols = [fitscol(array=max_redshift, format='D', name=zmax_colname),
				fitscol(array=max_distance, format='D', name=dmax_colname)]
	# create fits table from new+old columns
	for col in max_cols:
		if col.name in cat.columns.names:
			cat.columns.del_col(col.name)
	max_cols = fits.ColDefs(max_cols)
	hdu = fits.BinTableHDU.from_columns(max_cols + cat.columns)
	return hdu

def clone_galaxies(idcol, maxcol, Nrand=10, zg=None, Pz=None, zlims=None):
	# take arrays of IDs and max distances/redshifts
	# draw Nrand redshifts from quadratic/fitted distribution [0, {d/z}max] per galaxy
	Ngal = len(idcol)
	zmin, zmax = zlims or (0., 6.)
	dceil = cosmo.comoving_distance(zmax)
	randoms_id_z = np.empty([Ngal * Nrand, 2])

	# draw from fit to redshift distribution of real galaxies
	if zg is not None and Pz is not None:
		for i in tqdm(range(Nrand), desc='\t\tclone batches', ascii=True, ncols=100):
			zdraw = pchoice(zg, p=Pz, size=Ngal)
			mask = maxcol > 0
			nbad = (zdraw[mask] > maxcol[mask]).sum() # number of long-draws
			dnbad = nbad*1
			while nbad > 0:
				#print nbad
				badz = np.where((zdraw > maxcol) & mask)[0]
				if dnbad <= 10: # if re-draws are slowing down, loop over remaining draws with truncated P(z)
					for bz in badz:
						zg1, Pz1 = zg[zg < maxcol[bz]], Pz[zg < maxcol[bz]]
						Pz1 /= Pz1.sum()
						zdraw[bz] = pchoice(zg1, p=Pz1, size=1)[0]
				else: # otherwise keep re-drawing
					zdraw[badz] = pchoice(zg, p=Pz, size=len(badz))
				dnbad = nbad - (zdraw[mask] > maxcol[mask]).sum()
				nbad = (zdraw[mask] > maxcol[mask]).sum()
			zdraw[~mask] = -99.
			id_zdraw = np.column_stack((idcol.copy(), zdraw))
			randoms_id_z[i*Ngal:(i+1)*Ngal] = id_zdraw
		mask = maxcol > 0
		mask = np.array(list(mask) * Nrand)
		randoms_id_z[:, 1] = np.where(mask, randoms_id_z[:, 1], -99.)
		return randoms_id_z

	# draw from d^2 (== growth of volume element) between [0, d_max]
	else:
		Nrand_maxcol = np.array(list(maxcol) * Nrand)
		Nrand_idcol = np.array(list(idcol) * Nrand)
		ddraw = np.random.power(3, Ngal * Nrand) * Nrand_maxcol * u.Mpc
		#ddraw = np.where(ddraw > dceil, -99., ddraw)
		ddraw = ddraw.value
		mask = ddraw > 0
		#zmin = z_at_value(cosmo.comoving_distance, ddraw[mask].min())
		#zmax = z_at_value(cosmo.comoving_distance, ddraw[mask].max())
		#zgrid = np.logspace(log10(zmin), log10(zmax), 100)
		zgrid = np.linspace(zmin, zmax, 1000)
		dgrid = cosmo.comoving_distance(zgrid).value
		zdraw = np.ones_like(ddraw) * -99.
		zdraw[mask] = interp1d_(ddraw[mask], dgrid, zgrid)
#		print 'zmin =', zmin
#		print 'zmax =', zmax
#		print 'dceil =', dceil
#		print 'ddraw =', ddraw.min(), ddraw[mask].min(), ddraw[mask].mean(), ddraw[mask].max()

		randoms_id_z = np.column_stack((Nrand_idcol, zdraw))
		return randoms_id_z

#def mask_randoms(args):
	# expand randoms from 1x1 sky-box into catalogue footprint
	# create mutliple boxes for disjoint footprints?
	# use jackknife code to identify these?
	# apply a mask to resulting window

def main(args):
	print '\n'
	for cat_path in args.catalogues:
		print 'reading %s..' % cat_path
		cat = fits.open(cat_path)[1].data
		# establish zmax per detection band
		for mcol, mlim in zip(args.magcols, args.maglims):
			#ceil = cosmo.comoving_distance(args.zlims[1]).value
			hdu = get_zmax_hdu(cat, mcol, mlim, args.zcol)
			cat = hdu.data
		# save new columns
		hdu.writeto(cat_path, overwrite=1)

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

			# fit to n(z)
			#print '\t\tfitting redshift distribution..'
			#zg, Pz = fit_smail(zcol[mask], quiet=1, weights=wcol, bins=bins).T
			zg = Pz = None

			# construct redshift distribution per band
			print '\t\tcloning..'
			for mcol, mlim in zip(args.magcols, args.maglims):
				maxcol = cat[mcol+'_fl%.1f_dmax'%mlim]
				randoms_id_z = clone_galaxies(idcol, maxcol, args.Nrand, zg=zg, Pz=Pz, zlims=args.zlims)
				zmin, zmax = 0., randoms_id_z[:, 1].max()
				zgrid = np.linspace(zmin, zmax, 100)
				dgrid = cosmo.comoving_distance(zgrid)
				randoms_comoving = interp1d_(randoms_id_z[:, 1], zgrid, dgrid)
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
		default='give paths to fits catalogue(s) with redshifts')
	parser.add_argument(
		'-magcols',
		nargs='*',
		type=str,
		help='specify rest-frame magnitude column name(s) from which to compute Vmax(s)')
	parser.add_argument(
		'-maglims',
		nargs='*',
		type=float,
		help='specify apparent magnitude detection limit per -magcols column')
	parser.add_argument(
		'-randoms',
		type=int,
		help='1 = create a randoms catalogue with 1 redshift distribution per -magcol')
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

	main(args)





