from __future__ import print_function, division
import numpy as np
from astropy.io import fits, ascii
import sys
from os import listdir, mkdir
from os.path import basename, normpath, join, dirname, isdir
import argparse
from joblib import Parallel, delayed
import gc
import multiprocessing
from astropy.cosmology import FlatLambdaCDM as FLCDM
MICEcosmo = FLCDM(Om0=0.25, H0=70, Ob0=0.044)

def betwixt((re1,re2,de1,de2,ze1,ze2)):
        def make_cut(ra,dec,z):
                return (ra>=re1)&(ra<=re2)&(dec>=de1)&(dec<=de2)&(z>=ze1)&(z<=ze2)
        return make_cut

def read_randoms(path, cols=None, weights=None):
	print('reading randoms...')
	if '.fits' in path:
		columns = ('RA','DEC','Z')
		if cols!=None:
			columns = cols
		fd = fits.open(path)
		data = fd[1].data
	else:
		columns = (0,1,2)
		if cols!=None:
                        columns = cols
		data = np.loadtxt(path).T
	col1,col2,col3 = data[columns[0]],data[columns[1]],data[columns[2]]
	w = np.ones_like(col1)
	randoms = np.column_stack((col1,col2,col3,w))
	
	print('done.')
	return randoms

def read_reals(path, cols=None, weights=None):
	print('reading reals...')
	if '.fits' in path:
		columns = ('RA','DEC','Z')
		if cols!=None:
			columns = cols
		fd = fits.open(path)
                data = fd[1].data	
	else:
		columns = (0,1,2)
		if cols!=None:
			columns = cols
		data = np.loadtxt(path).T
	col1,col2,col3 = data[columns[0]],data[columns[1]],data[columns[2]]
	reals = np.column_stack((col1,col2,col3))
	if weights!=None:
		col4 = data[weights]
		reals = np.column_stack((reals,col4))
	
	print('done.')
	return reals

def save_reals(reals, path, zrange=None, rad=0):
	fpath, zrang, radfac = (0, 0, 0)
	if path.endswith('fits'):
		asc_out = path[:-4]+'asc'
		fpath = True
	elif path.endswith('asc'):
		asc_out = path
	if zrange!=None:
		zstring = '_z%s-%s'%(zrange[0],zrange[1])
		asc_out = asc_out[:-4]+zstring+'.asc'
		zrang = True
	colnames = ['# RA','DEC','z']
	if rad & any((reals.T[:2]>2*np.pi).flatten()):
		reals.T[:2] = reals.T[:2] * (np.pi/180)
		radfac = True
	if fpath|zrang|radfac:
		try:
			print('saving reals.asc..')
			ascii.write(reals, asc_out, names=colnames, delimiter='\t')
		except ValueError:
			print('..with weights..')
			colnames += ['weights']
			ascii.write(reals, asc_out, names=colnames, delimiter='\t')	
	else:
		print('nothing changed, not re-saving reals')

def cut_zrange(dataset,zmin,zmax):
	if dataset.ndim==1:
		new_dat = dataset[ (zmin<dataset) & (dataset<=zmax) ]
	else:
		zdat = dataset.T[2]
		new_dat = dataset[ (zmin<zdat) & (zdat<=zmax) ]
	return new_dat

def save_newrandoms(new_randoms, path, zrange=None, outfile=None, rad=0):
	if 'fits' in path:
		path = path[:-4]+'asc'
	rand_out = join(dirname(path),'downsampled_'+basename(normpath(path)))
	if outfile!=None:
		rand_out = outfile
	if zrange!=None:
		zstring = '_z%s-%s'%(zrange[0],zrange[1])
		rand_out = rand_out[:-4]+zstring+'.asc'
	if rad:
		new_randoms.T[:2] = new_randoms.T[:2] * (np.pi/180)
	print('saving randoms...')
	try:
        	ascii.write(new_randoms,rand_out,delimiter='\t',names=['# RA','DEC','z','weight'])
	except ValueError:
		new_randoms = np.column_stack((new_randoms, np.ones_like(new_randoms.T[0])))
        	ascii.write(new_randoms,rand_out,delimiter='\t',names=['# RA','DEC','z','weight'])

def find_limits(reals):
	ramin, decmin, zmin = (reals.T[i].min() for i in range(3))
        ramax, decmax, zmax = (reals.T[i].max() for i in range(3))
        xyzrange = ( ramax-ramin, decmax-decmin, zmax-zmin )
        minima = ( ramin, decmin, zmin )
	return xyzrange, minima

def make_randoms(reals, starting_factor=25):
	xyzrange, minima = find_limits(reals)
	real_size = reals.shape[0]
        # randomly, densely populate the survey cuboid
        print('laying down random points...')
        xmax,ymax,zmax = xyzrange
        X,Y,z = map(lambda x: np.random.random(starting_factor * real_size) * x, [xmax,ymax,zmax])
	X += minima[0]
	Y += minima[1]
	z += minima[2]
        rand_cat = np.column_stack((X,Y,z))
        return rand_cat

def trim_randoms(randoms, reals, max_sep=0.1):
	# idea is;
	# identify random points near/beyond real survey edges (somehow) = rpoints
	# for each rpoint, identify nearest real galaxy, and its distance from the rpoint
	# if >max_sep, remove rpoint ---> CPU heavvvyyyy??
	pass

def unit_check(cat, give_back='degrees', tag=''):
	assert give_back in ['degrees', 'radians'], 'give_back kwarg must == "degrees" | "radians"'
	coords = cat.T[:2]
	if all(coords.flatten() <= 2*np.pi):
		current = 'radians'
	else:
		current = 'degrees'
	if current != give_back:
		cat_ = cat.copy()
		#print('converting %s from %s to %s..'%(tag, current, give_back))
		if current == 'radians':
			cat_.T[:2] *= (180/np.pi)
		else:
			cat_.T[:2] *= (np.pi/180)
		return cat_
	else:
		return cat

def downsample(randoms, sample_z, nbin=12, target_nz=11):
	# downsample artificial randoms to match mocks' N(z)
	if all(sample_z < 1.):
		print('reals in REDSHIFT..')
		assert all(randoms.T[2] < 1.), 'RANDOMS vs. REALS z/chi mismatched! Patches in REDSHIFT'
		hist_range = (0., 0.6)
	else:
		print('reals in COMOVING DIST..')
		assert (not all(randoms.T[2] < 1.)), 'RANDOMS vs. REALS z/chi mismatched! Patches in COMOVING'
		hist_range = (0., 1580)
	hist_range = (sample_z.min(), sample_z.max())
	real_nz = np.histogram(sample_z, bins=nbin, range=hist_range)[0]
	real_Nz = real_nz/len(sample_z)
	print('real z.shape: ', sample_z.shape, '\nreal n(z): ', real_nz)

	# trim edges of z-distn
	random_z = randoms[:,2]
	zmin, zmax = sample_z.min(), sample_z.max() 
	ztrim = (random_z >= zmin) & (random_z <= zmax)
	print('trimming random-z edges to reals..')
	new_randoms, random_z = randoms[ztrim], random_z[ztrim]
	print('sample z minmax: %s, %s'%(sample_z.min(), sample_z.max()))
	print('random z minmax: %s, %s'%(random_z.min(), random_z.max()))

	rand_num = np.random.random(size=len(random_z))
	rand_nz = np.histogram(random_z, bins=nbin, range=hist_range)
	zbins = rand_nz[1]
	rand_nz = rand_nz[0]
	#print('rand z.shape: ', random_z.shape, '\nrand n(z): ', rand_nz)
	reduce_factor = target_nz*(real_nz/rand_nz)
	reduce_factor = np.nan_to_num(reduce_factor)

	print('downsampling random points...')
	nztune = np.zeros_like(random_z, dtype=bool)
	for i in range(len(rand_nz)):
		bincut = (zbins[i] < random_z) & (random_z <= zbins[i+1])
		if (reduce_factor[i] <= 0.95):
			nzcut = (rand_num >= 0.05) & (rand_num <= reduce_factor[i] + 0.05)
			nzbincut = bincut & nzcut
		else:
			print('random-z bin %s / %s dN/dz too low..!'%(i, len(reduce_factor)))
			nzbincut = bincut
		nztune = np.where(nzbincut, True, nztune)

	new_randoms = new_randoms[nztune]
	new_rand_z = new_randoms[:,2]
	rand_Nz = np.histogram(new_rand_z, bins=nbin, range=hist_range)[0]/len(new_rand_z)
	print('real/random N(z) (should be ~1):\n',real_Nz/rand_Nz)

#	else:
#		print('reduce_factor: ', reduce_factor)
#		new_rand_z = new_randoms[:,2]
#		print('no random downsampling required.')

	print('done.')
	return new_randoms

def para_jk_save(jkdir, patches, save_jks, jk_randoms, randoms, units, names, rand_names, j):
	# define JK samples, real and/or random, 1-per-patch
	jksample_str = join(jkdir, 'JKsample%s.asc'%(str(j).zfill(3)))
	jk_rand_str = join(jkdir, 'rand_JKsample%s.asc'%(str(j).zfill(3)))
		
	if save_jks:
		del_one_patches = np.delete(patches, j, axis=0)
		jksample = np.concatenate(del_one_patches)
		jksample = unit_check(jksample, give_back=units, tag='JKsample%s'%(str(j).zfill(3)))

		ascii.write(jksample, jksample_str, names=names, delimiter='\t')
		del jksample
		gc.collect()

	if jk_randoms:
		# find limits of patch
		patch_lims = np.array([( x.min(), x.max() ) for x in patches[j].T[:3]]).flatten()

		# match units with randoms
		if all(patches[j].T[:2].flatten() <= 2*np.pi) & ( not all(randoms.T[:2].flatten() <= 2*np.pi) ):
			patch_lims[:4] *= (180./np.pi)
		elif ( not all(patches[j].T[:2].flatten() <= 2*np.pi) )  &  all(randoms.T[:2].flatten() <= 2*np.pi):
			patch_lims[:4] *= (np.pi/180.)

		# reverse patch-cut for delete-1 jackknife
		patch_cut = betwixt( tuple(patch_lims) )
		jkrand_cut = ~patch_cut(randoms.T[0], randoms.T[1], randoms.T[2])

		jkrands = randoms[jkrand_cut]
		jkrands = unit_check(jkrands, give_back=units, tag='rand_JKsample%s'%(str(j).zfill(3)))

		# match N columns with reals
		Nmissing_cols = len(rand_names) - jkrands.shape[1]
		if Nmissing_cols>0:
			dummy_columns = np.ones([jkrands.shape[0], Nmissing_cols])
			jkrands = np.column_stack((jkrands, dummy_columns))
		elif Nmissing_cols<0:
			jkrands = jkrands[:, :len(rand_names)]

		ascii.write(jkrands, jk_rand_str, names=rand_names, delimiter='\t')
		del jkrands
		gc.collect()

def make_jks(wdir, randoms=None, random_cutter=None, empty_patches=None, radians=0, save_jks=0, jk_randoms=1, patch_str='patch', paths='all', largePi=0, sdss=0):
	# empty_patches is boolean array of length uncut-Npatches
	# True where all skinny-patch cuts are met
	# index 0-3: shapes, 4-7: densities, 8-9: all-colour densities
	patch_str = patch_str.split('*')
	units = ['degrees', 'radians'][radians]
	pathdict = {'all': ['highZ_Red_UnMasked', 'highZ_Blue_UnMasked', 'lowZ_Red_UnMasked', 'lowZ_Blue_UnMasked'],
		    'swot-all': ['swot_highZ_Red_UnMasked', 'swot_highZ_Blue_UnMasked', 'swot_lowZ_Red_UnMasked', 'swot_lowZ_Blue_UnMasked'],
		     'shapes': ['highZ_Red', 'highZ_Blue', 'lowZ_Red', 'lowZ_Blue']}

	samples = pathdict[paths]
	if sdss: samples = samples[:2]
	for s, sample in enumerate(samples):
		print('reading from: %s' % (['', '(____largePi)'][largePi]) )
		print(sample)

		pdir = join(wdir, sample)
		if largePi:
			pdir += '_largePi'
		if not isdir(pdir):
			print('no directory %s - skipping..'%pdir)
			continue

		copy_randoms = randoms.copy()
		if (paths!='swot-all') & (all(copy_randoms.T[2] <= 1.)) & ('swot' not in sample):
			print('converting random redshifts to comoving distances [Mpc/h]..')
			copy_randoms[:,2] = MICEcosmo.comoving_distance(copy_randoms.T[2]) * 0.7 # * h
			
		jkdir = join(pdir, 'JKsamples')
		if not isdir(jkdir):
			mkdir(jkdir)

		# read patches
		ldir = [i for i in listdir(pdir) if all([ps in i for ps in patch_str])]
		if len(ldir)==0:
			print('invalid patch_str for locating patch.asc files - give str patterns punctuated by *s')
			sys.exit()
		patches = np.array([np.loadtxt(join(pdir, li)) for li in ldir])
		print('N unmasked patches: %s'%patches.shape[0])

		# if SDSS, cut randoms by colour
		if sdss:
			redcut = copy_randoms.T[-1] > 0.63
			colour = np.array( [~redcut, redcut][int('Red' in sample)], dtype=bool )
			random_cutter_ = np.array( [np.array(rc, dtype=bool) & colour for rc in random_cutter], dtype=bool )
		else:
			random_cutter_ = np.array(random_cutter, dtype=bool)

		if paths=='all':
			ep_cut = np.array(empty_patches[s])
			print('N shapes patches: %s'%np.sum(ep_cut))
		elif paths=='swot-all':
			ep_cut = np.array(empty_patches[s+4])
			print('N swot patches: %s'%np.sum(ep_cut))

		# apply skinny-patch cut to random-cubes
		random_cutter_ = random_cutter_[ep_cut]

		# apply filtered cube-cuts to randoms & patches (cubes)
		combined_random_cutter = np.zeros(len(copy_randoms), dtype=bool)
		for random_cube in random_cutter_:
			combined_random_cutter |= np.array(random_cube, dtype=bool)
		copy_randoms = copy_randoms[combined_random_cutter]
		print('combined_random_cutter: ', sum(combined_random_cutter), '/', len(combined_random_cutter))
		#patches = patches[pop_patch_cut]

		# downsample randoms to patches
		patches_z = np.concatenate(patches, axis=0).T[2]
		print('downsampling patched randoms before JK sampling..')
		ds_randoms = downsample(copy_randoms, patches_z, nbin=(3, 1)[sdss], target_nz=10)
		print('ds_randoms.shape: ', ds_randoms.shape)

		names = ascii.read(join(pdir, ldir[0])).keys()
		names[:2] = [['# ra[deg]', 'dec[deg]'], ['# ra[rad]', 'dec[rad]']] [radians]
		rand_names = names # redundant??

		gc.collect()
		if save_jks | jk_randoms:
			num_cores = multiprocessing.cpu_count()
			Parallel(n_jobs=num_cores)(delayed(para_jk_save)(jkdir, patches, save_jks, jk_randoms, ds_randoms, units, names, rand_names, j) for j in range(len(patches)) )
		else:
			print('no JK function called.')

if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument(
	'reals_path',
	help='path of real galaxies for z-distribution, fits or ascii/csv')
	parser.add_argument(
	'-rand_path',
	default=None,
	help='path of randoms file to be downsampled, fits or ascii/csv, default=None; randoms will be created instead')
	parser.add_argument(
	'-random_cols',
	nargs=3,
	default=None,
	help='column identifiers for coords and redshifts in randoms. give strings as eg. "RA", or give column numbers. 3x args')
	parser.add_argument(
	'-real_cols',
	nargs=3,
	default=None,
	help='column identifiers for coords and redshifts in real data, for saving to ascii. give strings as eg. "Z", or give column numbers. ignore if input is ascii and redshift is 3rd column')
	parser.add_argument(
	'-weights',
	default=None,
	type=str,
	help='column identifier (key) for optional weights, if saving .fits file to ascii')
	parser.add_argument(
        '-zrange',
        default=None,
        nargs=2,
        type=np.float32,
        help='optionally specify limits to redshift range for saving of real/random ascii files, will be appended to filenames')
	parser.add_argument(
	'-outfile',
	type=str,
	default=None,
	help='FULL PATH & output filename.asc for downsampled randoms')
	parser.add_argument(
	'-nbin',
	default=9,
	type=int,
	help='number bins in which to downsample, default=9')
	parser.add_argument(
	'-target_nz',
	type=int,
	default=10,
	help='random/real density factor target, default=10')
	parser.add_argument(
	'-ds',
	type=int,
	default=1,
	help='perform downsampling (1), or just cut/save reals (0), default=1')
	parser.add_argument(
	'-readrand',
	help='read in existing randoms for downsampling (1) -> must give -path arg, or lay down random points in a wedge defined by real coord limits (0), default=1',
	type=int,
	default=1)
	parser.add_argument(
	'-radians',
	type=int,
	default=0,
	help='save .asc catalogs with coordinates in radians (1), or degrees (0), default=0')
	parser.add_argument(
	'-jkfunc',
	help='(1): jackknife randoms & save in JKdir, (2): jk reals only & save, (3): jk reals AND randoms, & save, or (0) none of the above! Running jkfunc bypasses all other downsampler utilities (must still pass some/anything to reals_path!!. Default=0',
	type=int,
	default=0)
	parser.add_argument(
	'-patch_str',
	help='(for jkfunc) patch string identifier, in shell format eg. patch*asc for all file strings containing "patch" and "asc". Default=patch',
	type=str,
	default='patch')
	parser.add_argument(
	'-patch_dirs',
	help='"all" for all z/colour unmasked samples (comoving dists), or "swot-all" for equivalent swot samples (redshifts)',
	default='all',
	type=str,
	choices=['all', 'swot-all'])
	parser.add_argument(
	'-wdir',
	type=str,
	default=None,
	help='Wcorr working directory containing samples/patches/jksamples/, required for save jackknife samples')
	args = parser.parse_args()

	if args.jkfunc!=0:
		if args.wdir==None:
			print('-wdir arg must specify Wcorr working directory for patch directories')
			sys.exit()
		jkf = args.jkfunc
		assert jkf in [1,2,3], 'jkfunc arg must be one of 0,1,2,3'
		jkfunc_dict = {1 : (0, 1), 2 : (1, 0), 3: (1, 1)}
		jkf = jkfunc_dict[jkf]
		make_jks(args.wdir, radians=args.radians, save_jks=jkf[0], jk_randoms=jkf[1], patch_str=args.patch_str, paths=args.patch_dirs)
		sys.exit()
		
	# read catalogues
	reals = read_reals(args.reals_path, args.real_cols, args.weights)
	if args.ds:
		if args.readrand:
			randoms = read_randoms(args.rand_path, args.random_cols)
		else:
			print('MAKING randoms..')
			print('NOTE: should only make_randoms for box-like real footprint - no functionality for curved edges/masks YET(?!)')
			randoms = make_randoms(reals)
			#randoms = trim_randoms(randoms, reals, maxsep=0.1)
	
	# cut redshift range
	if args.zrange!=None:
		zmin, zmax = args.zrange
		reals = cut_zrange(reals, zmin, zmax)	
	
	# resave real.fits as ascii
	save_reals(reals, args.reals_path, args.zrange, args.radians)

	# downsample randoms
	if args.ds:
		sample_z = reals[:,2]
		new_randoms = downsample(randoms, sample_z, args.nbin, args.target_nz)
		if args.zrange!=None:
			zmin,zmax = args.zrange
			new_randoms = cut_zrange(new_randoms, zmin, zmax)
		randpath = args.rand_path
		if args.rand_path==None:
			randpath = dirname(args.reals_path) + '/rand_' + basename(normpath(args.reals_path))
		save_newrandoms(new_randoms, randpath, args.zrange, args.outfile, args.radians)
	print('DONE.')






